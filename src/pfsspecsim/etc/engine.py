"""Orchestrator: ports `main()` from `src/gsetc.c` (gsetc.c:1690-2178).

`run_etc(params) -> EtcResults` builds the spectrograph configuration,
computes (or reloads) the per-arm noise vectors, and then -- for whichever
of the [OII]-curve/single-line/continuum/[OII]-catalog outputs `params`
requests (`None` skips that output entirely, mirroring the C engine's
``-``-named-file convention) -- assembles the corresponding
`astropy.table.Table`. `run_etc_files(params) -> EtcResults` additionally
writes each requested table to `params.outdir`-relative ECSV files via
`io.write_table`.

Index/arm-gating conventions (see the task brief's "実行者への注意"):
`ia` is always the internal 0-based arm loop index; only the *output*
``arm`` column uses `config.spectro_arm`'s id (which differs from `ia`
only for the red arm in MR mode). `noise.py`/`snr.py` deliberately leave
arm-range gating (whether a redshifted line's wavelength actually falls
inside a given arm's `[lmin, lmax)`) and the ``sqrt(n_exp)`` "coadd
n_exp exposures" rescaling to their caller -- this module.
"""

from __future__ import annotations

import dataclasses
import math
from pathlib import Path

import numpy as np
from astropy import units as u
from astropy.table import Table, vstack

from . import io, psf, snr
from .config import (
    MAXARM,
    Spectrograph,
    find_config_file,
    load_spectrograph_config,
    spectro_arm,
)
from .constants import DZ_OII, NZ_OII, OII_LAMBDA, Z_CHUNK, ZMIN_OII
from .noise import NoiseResult, compute_noise
from .params import EtcParams, MagSpec, resolve_degrade

#: Rest-frame midpoint of the [OII] 3727 doublet (gsetc.c's literal
#: ``372.845`` -- exactly ``(372.71 + 372.98) / 2``; derived from
#: `constants.OII_LAMBDA` rather than repeating the literal, per the
#: "no hand-typed duplicate constants" discipline).
_OII_MID = 0.5 * (OII_LAMBDA[0] + OII_LAMBDA[1])

#: Rest-frame wavelength of the single line curve's fiducial line
#: (gsetc.c:2067-2084's literal ``345.5``; not derivable from any other
#: constant -- an independent algorithm-defining literal of gsetc.c, kept
#: verbatim here since it is the SNL curve's *own* fiducial line, unrelated
#: to `EtcParams.line_flux`/`line_width`'s wavelength-independent meaning).
_SNL_LAMBDA_REST = 345.5

#: [OII]-doublet arm-range gate half-widths (gsetc.c:2036/2134/2153):
#: ``lmin[ia] < 373.8*(1+z)`` and ``371.8*(1+z) < lmax[ia]`` -- offset from
#: the doublet's actual rest wavelengths (372.71/372.98) by enough margin
#: that the *whole* NP_WIN-pixel deposit window (not just the line center)
#: stays inside the arm; not derived from `constants.OII_LAMBDA`, since the
#: margin itself is an independent, hand-tuned gsetc.c literal.
_OII_GATE_LO_REST = 371.8
_OII_GATE_HI_REST = 373.8


@dataclasses.dataclass
class EtcResults:
    """Result of one `run_etc` call; mirrors `main`'s in-memory state just
    before it would have written each output file.

    Attributes
    ----------
    noise : Table
        Always populated (computed fresh, or reloaded when
        `EtcParams.noise_reused`): concatenation of every arm's per-pixel
        noise variance/sky vectors (gsetc.c:1979-2017).
    snc, snl, oii_curve, oii_catalog : Table or None
        `None` exactly when the corresponding `EtcParams.outfile_*` (or,
        for `oii_catalog`, `EtcParams.oii_cat_in`) is `None` -- the
        ``-``-named-file "skip this output" convention (gsetc.c:1967-1978).
        `oii_curve` is the legacy "sno2" [OII]-doublet SNR-vs-redshift
        curve (gsetc.c:2019-2058); `snl` is the single-emission-line SNR
        curve (gsetc.c:2060-2093); `snc` is the continuum SNR curve
        (gsetc.c:2096-2116).
    aperture_factor_800_target, aperture_factor_800_point : float
        The two encircled-energy-at-800nm diagnostics `main` prints before
        doing any real work (gsetc.c:1959-1960): for the fiducial target
        (`EtcParams.reff`) and for a point source, respectively.
    oii_histogram : ndarray of shape (constants.NZ_OII,), or None
        Redshift-recovery histogram counts (gsetc.c:2171-2176), `None`
        unless `EtcParams.oii_cat_in` is set. Bin `j` covers
        ``[ZMIN_OII + j*DZ_OII, ZMIN_OII + (j+1)*DZ_OII)``.
    oii_n_targets : int or None
        gsetc.c's ``ngtot``: the number of detected catalog objects whose
        redshift also falls inside the histogram range (equivalently,
        ``oii_histogram.sum()``) -- *not* necessarily ``len(oii_catalog)``,
        since gsetc.c:2160-2165 counts ``ngtot`` only inside the
        ``j>=0 && j<NZ_OII`` histogram-bin check while still writing every
        detected object to the catalog (preserved quirk). `None` unless
        `EtcParams.oii_cat_in` is set.
    """

    noise: Table
    snc: Table | None
    snl: Table | None
    oii_curve: Table | None
    oii_catalog: Table | None
    aperture_factor_800_target: float
    aperture_factor_800_point: float
    oii_histogram: np.ndarray | None = None
    oii_n_targets: int | None = None


# --- z-grid "private hooks" -------------------------------------------------
#
# Each of these is exactly one `np.arange` call reproducing a `main` `for`
# loop's grid (not vectorized-for-performance -- just factored out so tests
# can `monkeypatch` them to a handful of points for a fast smoke run that
# still exercises the full column-assembly code path).


def _oii_curve_z_grid() -> np.ndarray:
    """z grid for the [OII] SNR curve (gsetc.c:2031: ``z=0.1; z<2.5001;
    z+=0.0002``)."""
    return np.arange(0.1, 2.5001, 0.0002)


def _snl_z_grid() -> np.ndarray:
    """z grid for the single-line SNR curve (gsetc.c:2067: ``z=0.1;
    z<2.7627; z+=0.0002``)."""
    return np.arange(0.1, 2.7627, 0.0002)


def _oii_prescan_z_grid() -> np.ndarray:
    """z grid for the [OII]-catalog MDLF prescan (gsetc.c:2130: ``z=0.1;
    z<2.5001; z+=0.00025``)."""
    return np.arange(0.1, 2.5001, 0.00025)


def _resolve_config_path(params: EtcParams) -> Path:
    """The spectrograph config file this run uses: `EtcParams.instr_config`
    if given explicitly, else the packaged file `config.find_config_file`
    resolves from `throughput_model`/`spectrograph`/`mr_mode`.
    """
    if params.instr_config is not None:
        return Path(params.instr_config)
    return find_config_file(
        params.throughput_model, params.spectrograph, params.mr_mode
    )


def _map_masked(
    n: int, mask: np.ndarray, chunk_size: int, fn, *arrays: np.ndarray
) -> np.ndarray:
    """Evaluate `fn(*arrays_at_masked_positions)` only where `mask` is
    True, chunked in groups of at most `chunk_size`; positions outside
    `mask` are left at 0.

    Mirrors gsetc.c's per-arm range gate (e.g. gsetc.c:2036, 2134, 2151),
    which skips the `gsGetSNR_OII`/`gsGetSNR_Single` call entirely rather
    than computing and discarding -- important not just for performance but
    for correctness on caller-supplied catalog rows, where an ungated
    `r_eff`/`FOII` might not even be physically valid input (gsetc.c's
    ``r_eff>=0.`` catalog eligibility check, gsetc.c:2151). Chunking bounds
    peak memory for large z/row sweeps (gsetc.c:2031/2067/2129 -- thousands
    of points): each chunk's `snr.snr_oii`/`snr.snr_single` call internally
    materializes a `(chunk_size, ~1000)`-shaped MTF array (see
    `constants.Z_CHUNK`).
    """
    out = np.zeros(n, dtype=np.float64)
    idx = np.flatnonzero(mask)
    for start in range(0, idx.size, chunk_size):
        sl = idx[start : start + chunk_size]
        out[sl] = fn(*(arr[sl] for arr in arrays))
    return out


# --- Noise: compute or reload (gsetc.c:1979-2017) --------------------------


def _build_noise_table(spectro: Spectrograph, noise_result: NoiseResult) -> Table:
    """Assemble the "noise" ECSV table (all arms concatenated) from a
    freshly computed `NoiseResult` (gsetc.c:2002-2017).
    """
    tables = []
    for ia in range(spectro.N_arms):
        arm = noise_result[ia]
        npix = arm.noise.size
        tables.append(
            Table(
                {
                    "arm": np.full(npix, spectro_arm(ia, spectro.MR), dtype=np.int64),
                    "pixel": np.arange(npix, dtype=np.int64),
                    "wavelength": arm.lam * u.nm,
                    "variance": arm.noise * u.electron**2,
                    "sky": arm.sky * u.electron,
                }
            )
        )
    return vstack(tables)


def _noise_arrays_from_table(
    spectro: Spectrograph, noise_table: Table
) -> tuple[list[np.ndarray], list[np.ndarray]]:
    """Regroup a (reloaded) noise Table into per-internal-arm noise/sky
    arrays, ordered by ascending `pixel`.

    Ports the *intent* of gsetc.c:1980-2000's noise-reload branch --
    "assign each row to its arm's array" -- but NOT its implementation:
    the C code assumes arm 0 and arm 1 always have equal pixel counts,
    using the fixed offsets ``spectro.npix[0]``/``spectro.npix[0]+
    spectro.npix[1]`` against a flat row counter `k` rather than grouping
    by the file's own `arm` column -- a latent bug for any spectrograph
    config where the three arms don't all share one `npix`. Grouping
    explicitly by the `arm` column (mapped back to the internal arm index
    via `config.spectro_arm`) sidesteps that bug entirely, as the task
    brief instructs.
    """
    arm_col = np.asarray(noise_table["arm"])
    pixel_col = np.asarray(noise_table["pixel"])
    variance_col = np.asarray(noise_table["variance"], dtype=np.float64)
    sky_col = np.asarray(noise_table["sky"], dtype=np.float64)

    noise_arrays = []
    sky_arrays = []
    for ia in range(spectro.N_arms):
        arm_id = spectro_arm(ia, spectro.MR)
        sel = arm_col == arm_id
        npix = int(spectro.npix[ia])
        n_sel = int(np.count_nonzero(sel))
        if n_sel != npix:
            raise ValueError(
                f"noise_reused: reloaded noise table has {n_sel} row(s) for "
                f"arm id {arm_id} (internal ia={ia}), expected {npix} from "
                "the spectrograph config"
            )
        order = np.argsort(pixel_col[sel])
        noise_arrays.append(variance_col[sel][order])
        sky_arrays.append(sky_col[sel][order])
    return noise_arrays, sky_arrays


# --- [OII] SNR curve, "sno2" (gsetc.c:2019-2058) ---------------------------


def _compute_oii_curve(
    params: EtcParams, spectro: Spectrograph, noise_arrays: list[np.ndarray]
) -> Table:
    z = _oii_curve_z_grid()
    n = z.size
    sqrt_nexp = math.sqrt(params.exp_num)

    wavelength0 = OII_LAMBDA[0] * (1.0 + z)
    wavelength1 = OII_LAMBDA[1] * (1.0 + z)
    lam_mid = _OII_MID * (1.0 + z)

    # QUIRK, preserved verbatim (gsetc.c:2049): the diagnostic aperture
    # factor column uses `fieldang=0`, unlike every other use of
    # `geometric_throughput` in this module (which use the real
    # `params.field_ang`) -- see the task brief and `snr.py`'s module
    # docstring.
    fiber_aperture_factor = psf.geometric_throughput(
        spectro, lam_mid, params.reff, params.fiber_offset, 0.0, params.seeing
    )

    aeff = np.zeros(n, dtype=np.float64)
    snr_arms = np.zeros((spectro.N_arms, n), dtype=np.float64)
    for ia in range(spectro.N_arms):
        # Aeff gate (gsetc.c:2045): centered on the doublet midpoint,
        # distinct from the SNR gate below.
        aeff_gate = (spectro.lmin[ia] < lam_mid) & (lam_mid < spectro.lmax[ia])
        aeff += np.where(
            aeff_gate, psf.effective_area(spectro, ia, lam_mid, params.field_ang), 0.0
        )

        snr_gate = (spectro.lmin[ia] < _OII_GATE_HI_REST * (1.0 + z)) & (
            _OII_GATE_LO_REST * (1.0 + z) < spectro.lmax[ia]
        )
        raw = _map_masked(
            n,
            snr_gate,
            Z_CHUNK,
            lambda zc, ia=ia: snr.snr_oii(
                params,
                spectro,
                ia,
                zc,
                params.line_flux,
                params.line_width,
                params.reff,
                0.0,
                1.0,
                noise_arrays[ia],
                snr_type=2,
            ),
            z,
        )
        snr_arms[ia] = raw * sqrt_nexp

    snr_tot = np.sqrt(np.sum(snr_arms**2, axis=0))

    return Table(
        {
            "z": z,
            "wavelength0": wavelength0 * u.nm,
            "wavelength1": wavelength1 * u.nm,
            "fiber_aperture_factor": fiber_aperture_factor,
            "effective_area": aeff * u.m**2,
            "snr_b": snr_arms[0],
            "snr_r": snr_arms[1],
            "snr_n": snr_arms[2],
            "snr_tot": snr_tot,
        }
    )


# --- Single-line SNR curve, "snl" (gsetc.c:2060-2093) ----------------------


def _compute_snl(
    params: EtcParams,
    spectro: Spectrograph,
    noise_arrays: list[np.ndarray],
    magspec: MagSpec,
) -> Table:
    z = _snl_z_grid()
    n = z.size
    sqrt_nexp = math.sqrt(params.exp_num)

    lam = _SNL_LAMBDA_REST * (1.0 + z)
    mag_at_lam = magspec(lam)

    # Unlike the [OII] curve's fiducial aperture factor, this one uses the
    # real `params.field_ang` (gsetc.c:2084 vs 2049 -- see task brief).
    fiber_aperture_factor = psf.geometric_throughput(
        spectro, lam, params.reff, params.fiber_offset, params.field_ang, params.seeing
    )

    aeff = np.zeros(n, dtype=np.float64)
    snr_arms = np.zeros((spectro.N_arms, n), dtype=np.float64)
    for ia in range(spectro.N_arms):
        gate = (spectro.lmin[ia] < lam) & (lam < spectro.lmax[ia])
        aeff += np.where(
            gate, psf.effective_area(spectro, ia, lam, params.field_ang), 0.0
        )

        raw = _map_masked(
            n,
            gate,
            Z_CHUNK,
            lambda zc, mc, ia=ia: snr.snr_single(
                params,
                spectro,
                ia,
                mc,
                _SNL_LAMBDA_REST * (1.0 + zc),
                params.line_flux,
                params.line_width,
                params.reff,
                noise_arrays[ia],
                snr_type=0,
            ),
            z,
            mag_at_lam,
        )
        snr_arms[ia] = raw * sqrt_nexp

    snr_tot = np.sqrt(np.sum(snr_arms**2, axis=0))

    mid_name = "snr_m" if spectro.MR else "snr_r"
    return Table(
        {
            "wavelength": lam * u.nm,
            "fiber_aperture_factor": fiber_aperture_factor,
            "effective_area": aeff * u.m**2,
            "snr_b": snr_arms[0],
            mid_name: snr_arms[1],
            "snr_n": snr_arms[2],
            "snr_tot": snr_tot,
        }
    )


# --- Continuum SNR curve, "snc" (gsetc.c:2096-2116) -------------------------


def _compute_snc(
    params: EtcParams,
    spectro: Spectrograph,
    noise_arrays: list[np.ndarray],
    sky_arrays: list[np.ndarray],
    magspec: MagSpec,
) -> Table:
    sqrt_nexp = math.sqrt(params.exp_num)
    tables = []
    for ia in range(spectro.N_arms):
        npix = int(spectro.npix[ia])
        # QUIRK, preserved verbatim (gsetc.c:2109, snr.snr_continuum's own
        # convention): pixel *left edge*, unlike the noise table's pixel
        # *center* -- see snr.py's module docstring.
        lam = spectro.lmin[ia] + spectro.dl[ia] * np.arange(npix, dtype=np.float64)
        mag_at_lam = magspec(lam)

        result = snr.snr_continuum(
            params, spectro, ia, mag_at_lam, params.reff, noise_arrays[ia]
        )

        tables.append(
            Table(
                {
                    "arm": np.full(npix, spectro_arm(ia, spectro.MR), dtype=np.int64),
                    "pixel": np.arange(npix, dtype=np.int64),
                    "wavelength": lam * u.nm,
                    "snr": result.snr * sqrt_nexp,
                    "signal": result.counts * u.electron,
                    "noise_variance": noise_arrays[ia] * u.electron**2,
                    "noise_variance_tot": result.noise * u.electron**2,
                    "input_mag": result.mag,
                    "conversion_factor": result.trans,
                    "sampling_factor": result.sample_factor,
                    "sky": sky_arrays[ia] * u.electron,
                }
            )
        )
    return vstack(tables)


# --- [OII] catalog + MDLF prescan + redshift histogram (gsetc.c:2119-2177) -


def _compute_mdlf(
    params: EtcParams, spectro: Spectrograph, noise_arrays: list[np.ndarray]
) -> float:
    """Minimum detectable line flux for point sources (gsetc.c:2124-2142):
    the fiducial (1:1 line ratio) [OII] flux at which the best-case (over
    the whole prescan z grid) SNR equals `min_snr`, divided by 1.6 (neither
    doublet line can exceed 80% of the total flux) and by a 0.9 safety
    factor.
    """
    z = _oii_prescan_z_grid()
    n = z.size
    sqrt_nexp = math.sqrt(params.exp_num)

    snr_arms = np.zeros((spectro.N_arms, n), dtype=np.float64)
    for ia in range(spectro.N_arms):
        gate = (spectro.lmin[ia] < _OII_GATE_HI_REST * (1.0 + z)) & (
            _OII_GATE_LO_REST * (1.0 + z) < spectro.lmax[ia]
        )
        raw = _map_masked(
            n,
            gate,
            Z_CHUNK,
            lambda zc, ia=ia: snr.snr_oii(
                params,
                spectro,
                ia,
                zc,
                1e-16,
                70.0,
                0.0,
                0.0,
                1.0,
                noise_arrays[ia],
                snr_type=2,
            ),
            z,
        )
        snr_arms[ia] = raw * sqrt_nexp / 1.6

    snrtot = np.sqrt(np.sum(snr_arms**2, axis=0))
    snrmax16 = max(1.0, float(np.max(snrtot)) if snrtot.size else 1.0)
    return 1e-16 / snrmax16 * params.min_snr * 0.9


def _read_oii_catalog(
    path: str | Path,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Read an [OII]-emitter input catalog (gsetc.c:2147's fscanf column
    order): ``id z r_eff ROII FOII contOII sigma``. `sigma` is read but
    never used (gsetc.c always passes a hardcoded ``sigma_v=70`` instead,
    gsetc.c:2154 -- see `snr.py`'s module docstring); it is dropped here
    rather than threaded through unused.
    """
    data = np.atleast_2d(np.loadtxt(path))
    obj_id = data[:, 0].astype(np.int64)
    z = data[:, 1]
    r_eff = data[:, 2]
    roii = data[:, 3]
    foii = data[:, 4]
    cont_oii = data[:, 5]
    return obj_id, z, r_eff, roii, foii, cont_oii


def _compute_oii_catalog(
    params: EtcParams, spectro: Spectrograph, noise_arrays: list[np.ndarray]
) -> tuple[Table, np.ndarray, int]:
    obj_id, z, r_eff, roii, foii, cont_oii = _read_oii_catalog(params.oii_cat_in)
    n = z.size
    sqrt_nexp = math.sqrt(params.exp_num)

    mdlf = _compute_mdlf(params, spectro, noise_arrays)
    # gsetc.c:2151 eligibility gate: objects below the point-source MDLF or
    # with an invalid (negative) r_eff never even get a snrtot computed
    # (implicitly 0, so they fail the `min_snr` cut below regardless).
    eligible = (foii >= mdlf) & (r_eff >= 0.0)

    snr_arms = np.zeros((spectro.N_arms, n), dtype=np.float64)
    for ia in range(spectro.N_arms):
        gate = (
            eligible
            & (spectro.lmin[ia] < _OII_GATE_HI_REST * (1.0 + z))
            & (_OII_GATE_LO_REST * (1.0 + z) < spectro.lmax[ia])
        )
        # `snr.snr_oii`'s `r_eff` is a *scalar* per call, exactly like the
        # C original's (gsetc.c:2154 calls gsGetSNR_OII once per catalog
        # row with that row's r_eff; `psf.geometric_throughput` likewise
        # takes one r_eff per call) -- so gated rows are grouped by unique
        # r_eff value and each group vectorized over its rows (chunked by
        # Z_CHUNK, as elsewhere). Real catalogs quote r_eff at coarse
        # precision (the C output format is %5.3lf), so groups are
        # typically large; the worst case degenerates to the C code's own
        # one-call-per-row behavior.
        gated = np.flatnonzero(gate)
        for r_val in np.unique(r_eff[gated]):
            rows = gated[r_eff[gated] == r_val]
            for start in range(0, rows.size, Z_CHUNK):
                sl = rows[start : start + Z_CHUNK]
                snr_arms[ia][sl] = snr.snr_oii(
                    params,
                    spectro,
                    ia,
                    z[sl],
                    foii[sl],
                    70.0,
                    float(r_val),
                    cont_oii[sl],
                    roii[sl],
                    noise_arrays[ia],
                    snr_type=2,
                )
        snr_arms[ia] *= sqrt_nexp

    snrtot = np.sqrt(np.sum(snr_arms**2, axis=0))
    detected = snrtot >= params.min_snr

    histogram = np.zeros(NZ_OII, dtype=np.int64)
    bin_idx = np.floor((z - ZMIN_OII) / DZ_OII).astype(np.int64)
    in_hist_range = detected & (bin_idx >= 0) & (bin_idx < NZ_OII)
    np.add.at(histogram, bin_idx[in_hist_range], 1)
    # gsetc.c:2160-2165 QUIRK, preserved verbatim: `ngtot++` sits *inside*
    # the `if (j>=0 && j<NZ_OII)` histogram-range check, so a detected
    # object whose redshift falls outside [ZMIN_OII, ZMIN_OII+NZ_OII*DZ_OII)
    # is written to the output catalog but counted in neither ngal[] nor
    # ngtot -- i.e. n_targets == histogram.sum(), which may be less than
    # len(catalog). Reachable: with the packaged LR config (blue arm
    # lmin=380nm) an [OII] emitter at z ~ 0.02-0.099 passes the arm gate
    # (373.8*(1+z) > 380) yet bins to j < 0.
    n_targets = int(np.count_nonzero(in_hist_range))

    catalog = Table(
        {
            "obj_id": obj_id[detected],
            "z": z[detected],
            "reff": r_eff[detected] * u.arcsec,
            "flux_oii": foii[detected] * u.erg / u.s / u.cm**2,
            "snr": snrtot[detected],
        }
    )
    catalog.meta["redshift_histogram"] = {
        "zmin": ZMIN_OII,
        "dz": DZ_OII,
        "n_bins": NZ_OII,
        "counts": histogram.tolist(),
        "n_targets": n_targets,
    }
    return catalog, histogram, n_targets


# --- Top-level API -----------------------------------------------------


def run_etc(params: EtcParams) -> EtcResults:
    """Run the full ETC pipeline for `params` and return the results;
    writes no files (see `run_etc_files` for that).

    Port of `main` (gsetc.c:1690-2177), minus argument parsing (folded
    into `EtcParams`) and file I/O (folded into `run_etc_files`/`io.py`).
    """
    params.validate()

    config_path = _resolve_config_path(params)
    degrade = resolve_degrade(params)
    spectro = load_spectrograph_config(config_path, degrade)
    if spectro.N_arms != MAXARM:
        raise NotImplementedError(
            "run_etc: the SNL/[OII]-curve ECSV schemas hardcode 3 arm "
            f"columns (snr_b/snr_r-or-snr_m/snr_n); got N_arms="
            f"{spectro.N_arms} from {config_path}"
        )
    magspec = MagSpec(params.mag, params.mag_file)

    # Encircled energy in fiber at 800nm (gsetc.c:1959-1960).
    aperture_factor_800_target = float(
        psf.geometric_throughput(
            spectro,
            800.0,
            params.reff,
            params.fiber_offset,
            params.field_ang,
            params.seeing,
        )
    )
    aperture_factor_800_point = float(
        psf.geometric_throughput(
            spectro, 800.0, 0.0, params.fiber_offset, params.field_ang, params.seeing
        )
    )

    # Generate and write noise vector, or reload (gsetc.c:1979-2017).
    if params.noise_reused:
        if params.outfile_noise is None:
            raise ValueError(
                "params.noise_reused=True requires params.outfile_noise to "
                "name the (already-written) noise ECSV to reload"
            )
        noise_path = Path(params.outdir) / params.outfile_noise
        noise_table = io.read_table(noise_path)
        noise_arrays, sky_arrays = _noise_arrays_from_table(spectro, noise_table)
    else:
        noise_result = compute_noise(params, spectro)
        noise_arrays = [arm.noise for arm in noise_result.arms]
        sky_arrays = [arm.sky for arm in noise_result.arms]
        noise_table = _build_noise_table(spectro, noise_result)

    oii_curve = (
        _compute_oii_curve(params, spectro, noise_arrays)
        if params.outfile_oii is not None
        else None
    )
    snl = (
        _compute_snl(params, spectro, noise_arrays, magspec)
        if params.outfile_snl is not None
        else None
    )
    snc = (
        _compute_snc(params, spectro, noise_arrays, sky_arrays, magspec)
        if params.outfile_snc is not None
        else None
    )

    oii_catalog = None
    oii_histogram = None
    oii_n_targets = None
    if params.oii_cat_in is not None:
        oii_catalog, oii_histogram, oii_n_targets = _compute_oii_catalog(
            params, spectro, noise_arrays
        )

    return EtcResults(
        noise=noise_table,
        snc=snc,
        snl=snl,
        oii_curve=oii_curve,
        oii_catalog=oii_catalog,
        aperture_factor_800_target=aperture_factor_800_target,
        aperture_factor_800_point=aperture_factor_800_point,
        oii_histogram=oii_histogram,
        oii_n_targets=oii_n_targets,
    )


def run_etc_files(params: EtcParams) -> EtcResults:
    """Run `run_etc(params)` and write each requested output as an ECSV
    file (`io.write_table`), under `params.outdir` for the
    `outfile_noise`/`outfile_snc`/`outfile_snl`/`outfile_oii` fields
    (documented as outdir-relative in `EtcParams`); `oii_cat_out`, like the
    legacy `INFILE_OIICat`/`OUTFILE_OIICat` pair it replaces, is used as
    given (not joined with `outdir`) since it is not documented as
    outdir-relative and the legacy wrapper never treated it as such.

    `params.overwrite=False` existence check, a verbatim port of
    `pfsspecsim/pfsetc.py:236-249`: only `outfile_noise`/`outfile_snc`/
    `outfile_snl` are checked for pre-existence, and the check runs (and
    raises `FileExistsError`) *before* `run_etc` does any computation --
    the old wrapper aborted before ever invoking the (expensive) C
    subprocess. `outfile_oii` and `oii_cat_out` were never covered by that
    guard in the old wrapper either (the check simply never inspected
    them), so they are not checked here -- a preserved quirk, not an
    oversight: once past the guard, every write below always overwrites
    (matching the C engine's unconditional `fopen(path, "w")`) -- except
    `outfile_noise` when `params.noise_reused` is set, which is *never*
    (re)written: the C engine's reload branch (gsetc.c:1980-2001) only
    reads the noise file, leaving its content and mtime untouched. Corollary
    (also verbatim legacy behavior): `overwrite=False` combined with
    `noise_reused=True` always raises, since the noise ECSV must already
    exist to be reloaded, and its existence is exactly what trips the guard.
    """
    outdir = Path(params.outdir)
    noise_path = (
        outdir / params.outfile_noise if params.outfile_noise is not None else None
    )
    snc_path = outdir / params.outfile_snc if params.outfile_snc is not None else None
    snl_path = outdir / params.outfile_snl if params.outfile_snl is not None else None
    oii_curve_path = (
        outdir / params.outfile_oii if params.outfile_oii is not None else None
    )
    oii_cat_out_path = (
        Path(params.oii_cat_out) if params.oii_cat_out is not None else None
    )

    if not params.overwrite:
        existing = [
            p for p in (noise_path, snc_path, snl_path) if p is not None and p.exists()
        ]
        if existing:
            raise FileExistsError(
                "params.overwrite=False and output file(s) already exist: "
                + ", ".join(str(p) for p in existing)
            )

    results = run_etc(params)
    config_path = _resolve_config_path(params)

    # gsetc.c:1980-2001 (flag_reused==1) only *reads* the noise file, it
    # never rewrites it -- so skip the write when reloading. Rewriting
    # here would be a no-op for the data (the reloaded table is written
    # back unchanged) but would stamp the *current* call's params/meta and
    # mtime onto a file whose contents were produced by an earlier run
    # with possibly different params (e.g. the legacy chained pattern
    # `make_noise_model(); set_param(...); make_snc()`). The
    # overwrite=False pre-check above is unaffected (with
    # noise_reused=True it necessarily raised already -- verbatim legacy
    # behavior, see the docstring).
    if noise_path is not None and not params.noise_reused:
        io.write_table(results.noise, noise_path, params, config_path, overwrite=True)
    if snc_path is not None and results.snc is not None:
        io.write_table(results.snc, snc_path, params, config_path, overwrite=True)
    if snl_path is not None and results.snl is not None:
        io.write_table(results.snl, snl_path, params, config_path, overwrite=True)
    if oii_curve_path is not None and results.oii_curve is not None:
        io.write_table(
            results.oii_curve, oii_curve_path, params, config_path, overwrite=True
        )
    if oii_cat_out_path is not None and results.oii_catalog is not None:
        io.write_table(
            results.oii_catalog, oii_cat_out_path, params, config_path, overwrite=True
        )

    return results
