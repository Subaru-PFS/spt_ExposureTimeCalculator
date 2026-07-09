# -*- coding: utf-8 -*-
"""Deprecated compatibility layer for the legacy stdin-driven C ETC engine.

`Etc` used to drive `bin/gsetc[_omp].x` as a subprocess, feeding it 27
values over stdin and reading its plain-text output files back with
`np.genfromtxt`. The C engine has been replaced by the pure-Python
`pfsspecsim.etc` package (`pfsspecsim.etc.engine.run_etc_files`); this
module is now a thin translation shim from the old ALL_CAPS `params` dict
(and its 'Y'/'N'/'-' string conventions) to `pfsspecsim.etc.params.EtcParams`
and back, so that existing scripts/notebooks built around
``Etc(); set_param(...); run(); get_noise()`` keep working unmodified.

New code should use `pfsspecsim.etc` directly (`EtcParams`/`load_params`,
`engine.run_etc_files`) or the `pfs-etc` CLI -- see the project README's
old -> new parameter migration table. This module (and the `Etc` class in
particular) is deprecated and may be removed in a future release; every
`Etc()` construction emits a `DeprecationWarning`.

The output *files* are still written at exactly the paths the caller
configures via `OUTFILE_NOISE`/`OUTFILE_SNC`/`OUTFILE_SNL`/`OUTFILE_OII`
(unchanged naming), but their *content* is now Astropy ECSV (with a
`table.meta['params']` block recording every resolved input), not the old
whitespace-delimited plain text.
"""

from __future__ import annotations

import time
import warnings
from pathlib import Path

import numpy as np

from .etc import engine
from .etc.params import EtcParams, calc_obscuration

__all__ = ["Etc", "calc_obscuration"]

# Historical hardcoded constants from the old subprocess wrapper
# (pfsetc.py:17-18). The new engine turns both into `EtcParams` fields
# (`sky_type`, `fiber_offset`), but the old ALL_CAPS `params` dict never
# exposed them as settable parameters, so this compatibility layer keeps
# using the same fixed values -- which happen to equal `EtcParams`'s own
# defaults for those two fields.
SKYMODELS = "11006"
OFFSET_FIB = "0.10"


def _to_bool(value) -> bool:
    """Old 'Y'/'N' (also accepts 'yes'/'no') string convention -> bool."""
    return str(value).strip().lower() in ("y", "yes")


def _to_path_or_none(value):
    """Old '-' sentinel convention (`INFILE_OIICat`, `OUTFILE_*`, ...) ->
    `None`; anything else -> `Path`."""
    return None if str(value) == "-" else Path(value)


class Etc:
    """Deprecated compatibility wrapper around `pfsspecsim.etc`.

    .. deprecated::
        `Etc` reproduces the surface of the old subprocess-driven ETC
        wrapper (`params` dict with ALL_CAPS keys, `set_param`/`run`/
        `make_snc`/etc.) on top of the pure-Python `pfsspecsim.etc` engine.
        New code should call `pfsspecsim.etc.params.load_params` +
        `pfsspecsim.etc.engine.run_etc_files` directly, or use the
        `pfs-etc` console script. Constructing `Etc` emits a
        `DeprecationWarning`.

    See https://github.com/Subaru-PFS/spt_ExposureTimeCalculator for
    details.

    Parameters
    ----------
    omp_num_threads : `int` (default: 16)
        Deprecated: the pure-Python engine has no OpenMP thread pool to
        size, but `omp_num_threads` is still honored, mapped to the new
        engine's `pfsspecsim.etc.params.EtcParams.n_workers` (arm-level
        `ThreadPoolExecutor` parallelism), capped at 3 -- the number of
        spectrograph arms, beyond which more workers cannot help.

    Examples
    ----------
    """

    def __init__(self, omp_num_threads=16):
        warnings.warn(
            "pfsspecsim.pfsetc.Etc is deprecated; it now wraps the "
            "pure-Python pfsspecsim.etc engine. Use "
            "pfsspecsim.etc.params.load_params + "
            "pfsspecsim.etc.engine.run_etc_files (or the `pfs-etc` CLI) "
            "in new code. `omp_num_threads` is now mapped to the new "
            "engine's `n_workers` (arm-parallel thread count, capped at 3) "
            "instead of being ignored.",
            DeprecationWarning,
            stacklevel=2,
        )
        self.params = {
            "SEEING": "0.80",
            "ZENITH_ANG": "35.0",
            "GALACTIC_EXT": "0.00",
            "MOON_ZENITH_ANG": "30.0",
            "MOON_TARGET_ANG": "60.0",
            "MOON_PHASE": "0.125",
            "EXP_TIME": "900",
            "EXP_NUM": "4",
            "FIELD_ANG": "0.45",
            "MAG_FILE": "22.5",
            "REFF": "0.3",
            "LINE_FLUX": "1.0e-17",
            "LINE_WIDTH": "70",
            "NOISE_REUSED": "N",
            "MR_MODE": "N",
            "OVERWRITE": "Y",
            "INFILE_OIICat": "-",
            "OUTFILE_OIICat": "-",
            "minSNR": "9.0",
            "degrade": "1.0",
            "SKY_SUB_FLOOR": "0.01",
            "DIFFUSE_STRAY": "0.02",
            "throughput_model": "20240714",
            "spectrograph": "ave",
            "OUTDIR": "out",
            "TMPDIR": "tmp",
            "BINDIR": "bin",
            "obscFoVDep": "Y",
        }
        outdir = Path(self.params["OUTDIR"])
        self.params["OUTFILE_NOISE"] = str(outdir / "ref.noise.dat")
        self.params["OUTFILE_SNC"] = str(outdir / "ref.snc.dat")
        self.params["OUTFILE_SNL"] = str(outdir / "ref.snl.dat")
        self.params["OUTFILE_OII"] = str(outdir / "ref.sno2.dat")

        # Deprecated attribute, kept readable for backward compatibility;
        # `_to_new_params` maps it to `EtcParams.n_workers` (arm-parallel
        # thread count, clamped to [1, 3]).
        self.omp_num_threads = omp_num_threads

    def set_param(self, param_name, param_value):
        """Set ETC parameters

        Parameters
        ----------
        param_name : `str` (the parameter name)
        param_value:       (the parameter value)

        Returns
        -------

        Examples
        ----------
        """
        if param_name in self.params:
            try:
                self.params[param_name] = str(param_value)
            except Exception:
                print("Error!")
        else:
            print(f"param_name {param_name} can not be recognized ...")
        return 0

    def load_param_file(self, filename):
        with open(filename, "r") as f:
            for line in f:
                a = line.split()
                if line[0] != "#" and len(a) > 0:
                    self.params[a[0]] = a[1]
        return 0

    def _to_new_params(self, **field_overrides) -> EtcParams:
        """Translate the legacy `self.params` dict into an `EtcParams`.

        Follows the old -> new parameter migration table (see the project
        README / task brief): `MAG_FILE` splits into the mutually
        exclusive `mag` (if it parses as a `float`) / `mag_file` (the path,
        otherwise) -- reproducing pfsetc.py:217-233's old
        float-succeeds-vs-fails branch, but without ever writing a
        temporary resampled file (the new `MagSpec` handles that
        in-memory). `'-'` -> `None`; `'Y'`/`'N'` (and `'yes'`/`'no'`) ->
        `bool`. `TMPDIR`/`BINDIR` have no effect and are not consulted --
        historical no-ops preserved only as dict keys.

        `outdir` is always `'.'` and every `OUTFILE_*`/`OUTFILE_OIICat`
        path is passed through unchanged (including any leading `OUTDIR`,
        e.g. the default `'out/ref.snc.dat'`), so files land at exactly
        the path the caller configured -- matching the old wrapper, which
        always honored the literal `OUTFILE_*` string regardless of
        `OUTDIR`.

        `field_overrides` lets `make_noise_model`/`make_snc`/`make_snl`/
        `make_sno2` reproduce the old wrapper's per-method output
        selection and `NOISE_REUSED` overrides (pfsetc.py:397-490,
        548-578, 627-657, 705-734) on top of the common translation.

        `omp_num_threads` (the `__init__` argument, stored on `self`) maps
        to `EtcParams.n_workers`, capped at 3 (the number of spectrograph
        arms -- see `__init__`'s docstring).
        """
        p = self.params
        mag_raw = p["MAG_FILE"]
        try:
            mag = float(mag_raw)
            mag_file = None
        except (TypeError, ValueError):
            mag = None
            mag_file = Path(mag_raw)

        kwargs = dict(
            seeing=float(p["SEEING"]),
            zenith_ang=float(p["ZENITH_ANG"]),
            galactic_ext=float(p["GALACTIC_EXT"]),
            field_ang=float(p["FIELD_ANG"]),
            fiber_offset=float(OFFSET_FIB),
            moon_zenith_ang=float(p["MOON_ZENITH_ANG"]),
            moon_target_ang=float(p["MOON_TARGET_ANG"]),
            moon_phase=float(p["MOON_PHASE"]),
            exp_time=float(p["EXP_TIME"]),
            exp_num=int(float(p["EXP_NUM"])),
            mag=mag,
            mag_file=mag_file,
            reff=float(p["REFF"]),
            line_flux=float(p["LINE_FLUX"]),
            line_width=float(p["LINE_WIDTH"]),
            mr_mode=_to_bool(p["MR_MODE"]),
            throughput_model=p["throughput_model"],
            spectrograph=p["spectrograph"],
            degrade=float(p["degrade"]),
            obsc_fov_dep=_to_bool(p["obscFoVDep"]),
            sky_type=SKYMODELS,
            sky_sub_floor=float(p["SKY_SUB_FLOOR"]),
            diffuse_stray=float(p["DIFFUSE_STRAY"]),
            oii_cat_in=_to_path_or_none(p["INFILE_OIICat"]),
            oii_cat_out=_to_path_or_none(p["OUTFILE_OIICat"]),
            min_snr=float(p["minSNR"]),
            noise_reused=_to_bool(p["NOISE_REUSED"]),
            overwrite=_to_bool(p["OVERWRITE"]),
            outdir=Path("."),
            outfile_noise=_to_path_or_none(p["OUTFILE_NOISE"]),
            outfile_snc=_to_path_or_none(p["OUTFILE_SNC"]),
            outfile_snl=_to_path_or_none(p["OUTFILE_SNL"]),
            outfile_oii=_to_path_or_none(p["OUTFILE_OII"]),
            # Clamp to [1, 3]: legacy callers could pass any int (0 or
            # negative was a harmless no-op in the subprocess era) and must
            # not trip EtcParams.validate()'s n_workers >= 1 check.
            n_workers=max(1, min(self.omp_num_threads, 3)),
        )
        kwargs.update(field_overrides)
        return EtcParams(**kwargs)

    # -- Old-attribute restoration -------------------------------------
    #
    # Ports pfsetc.py's np.genfromtxt column unpacking (run: 297-384,
    # make_noise_model: 494-510, make_snc: 581-605, make_snl: 661-681,
    # make_sno2: 738-760) to read the same columns back out of the new
    # engine's Astropy Tables instead of re-parsing plain text. Column
    # order/semantics verified 1:1 against gsetc.c's fprintf calls
    # (gsetc.c:2010, 2048-2053, 2084-2087, 2109-2110).

    def _restore_noise(self, results):
        tbl = results.noise
        self.nsm_arms = np.asarray(tbl["arm"])
        self.nsm_pixs = np.asarray(tbl["pixel"])
        self.nsm_lams = np.asarray(tbl["wavelength"])
        self.nsm_nois = np.asarray(tbl["variance"])
        self.nsm_skys = np.asarray(tbl["sky"])

    def _restore_snc(self, results):
        tbl = results.snc
        if tbl is None:
            return
        self.snc_arms = np.asarray(tbl["arm"])
        self.snc_pixs = np.asarray(tbl["pixel"])
        self.snc_lams = np.asarray(tbl["wavelength"])
        self.snc_sncs = np.asarray(tbl["snr"])
        self.snc_sigs = np.asarray(tbl["signal"])
        self.snc_nois_mobj = np.asarray(tbl["noise_variance"])
        self.snc_nois = np.asarray(tbl["noise_variance_tot"])
        self.snc_spin = np.asarray(tbl["input_mag"])
        self.snc_conv = np.asarray(tbl["conversion_factor"])
        self.snc_samp = np.asarray(tbl["sampling_factor"])
        self.snc_skys = np.asarray(tbl["sky"])

    def _restore_snl(self, results):
        tbl = results.snl
        if tbl is None:
            return
        self.snl_lams = np.asarray(tbl["wavelength"])
        self.snl_fcov = np.asarray(tbl["fiber_aperture_factor"])
        self.snl_effa = np.asarray(tbl["effective_area"])
        self.snl_sna0 = np.asarray(tbl["snr_b"])
        mid_col = "snr_m" if "snr_m" in tbl.colnames else "snr_r"
        self.snl_sna1 = np.asarray(tbl[mid_col])
        self.snl_sna2 = np.asarray(tbl["snr_n"])
        self.snl_snls = np.asarray(tbl["snr_tot"])

    def _restore_sno2(self, results):
        tbl = results.oii_curve
        if tbl is None:
            return
        self.sno2_zsps = np.asarray(tbl["z"])
        self.sno2_lam1 = np.asarray(tbl["wavelength0"])
        self.sno2_lam2 = np.asarray(tbl["wavelength1"])
        self.sno2_fcov = np.asarray(tbl["fiber_aperture_factor"])
        self.sno2_effa = np.asarray(tbl["effective_area"])
        self.sno2_sna0 = np.asarray(tbl["snr_b"])
        self.sno2_sna1 = np.asarray(tbl["snr_r"])
        self.sno2_sna2 = np.asarray(tbl["snr_n"])
        self.sno2_sno2 = np.asarray(tbl["snr_tot"])

    def run(self):
        start = time.time()

        print("##### starting to run ETC ... (it takes a few min.) #####")
        params = self._to_new_params()
        results = engine.run_etc_files(params)

        self._restore_noise(results)
        self._restore_snc(results)
        self._restore_snl(results)
        self._restore_sno2(results)

        elapsed_time = time.time() - start
        print(f"##### finished (elapsed_time: {elapsed_time:.1f}[sec]) #####")

        return 0

    def make_noise_model(self):
        start = time.time()

        print("##### starting to make a noise model ... (it takes about 2 min.) #####")
        # QUIRK, preserved verbatim (pfsetc.py:478 hardcoded '0' regardless
        # of `self.params['NOISE_REUSED']`): making a fresh noise model
        # always (re)computes it, never reloads.
        params = self._to_new_params(
            noise_reused=False,
            outfile_snc=None,
            outfile_snl=None,
            outfile_oii=None,
            oii_cat_in=None,
            oii_cat_out=None,
        )
        results = engine.run_etc_files(params)
        self._restore_noise(results)

        elapsed_time = time.time() - start
        print(f"##### finished (elapsed_time: {elapsed_time:.1f}[sec]) #####")
        return 0

    def proc_multi(self, inputs):
        self.params[inputs[0]] = str(inputs[1])
        for outFileName in [
            "OUTFILE_NOISE",
            "OUTFILE_SNC",
            "OUTFILE_SNL",
            "OUTFILE_OII",
        ]:
            if self.params[outFileName] != "-":
                self.params[outFileName] += ".".join(["", inputs[0], inputs[1]])
        self.run()
        return 0

    def run_multi(self, nproc, param_name, param_values):
        from multiprocessing import Pool

        p = Pool(nproc)
        result = p.map(self.proc_multi, [(param_name, v) for v in param_values])
        return 0

    def get_noise(self):
        return self.nsm_lams, self.nsm_nois

    def make_snc(self):
        # QUIRK, preserved verbatim (pfsetc.py:566 hardcoded '1' regardless
        # of `self.params['NOISE_REUSED']`): making an SNC curve always
        # reloads the noise vector previously written to `OUTFILE_NOISE`
        # (by `run()` or `make_noise_model()`) rather than recomputing it.
        start = time.time()
        print("##### starting to make an SNC model ... (it takes about 1 min.) #####")
        params = self._to_new_params(
            noise_reused=True,
            outfile_snl=None,
            outfile_oii=None,
            oii_cat_in=None,
            oii_cat_out=None,
        )
        results = engine.run_etc_files(params)
        self._restore_snc(results)

        elapsed_time = time.time() - start
        print(f"##### finished (elapsed_time: {elapsed_time:.1f}[sec]) #####")

        return 0

    def get_snc(self):
        return self.snc_lams, self.snc_sncs

    def make_snl(self):
        # QUIRK, preserved verbatim (pfsetc.py:645, same as `make_snc`):
        # always reloads the noise vector rather than recomputing it.
        start = time.time()
        print("##### starting to make an SNL model ... (it takes about 1 min.) #####")
        params = self._to_new_params(
            noise_reused=True,
            outfile_snc=None,
            outfile_oii=None,
            oii_cat_in=None,
            oii_cat_out=None,
        )
        results = engine.run_etc_files(params)
        self._restore_snl(results)

        elapsed_time = time.time() - start
        print(f"##### finished (elapsed_time: {elapsed_time:.1f}[sec]) #####")

        return 0

    def get_snl(self):
        return self.snl_lams, self.snl_snls

    def make_sno2(self):
        # QUIRK, preserved verbatim (pfsetc.py:722, same as `make_snc`):
        # always reloads the noise vector rather than recomputing it.
        start = time.time()
        print("##### starting to make an OII model ... (it takes about 2 min.) #####")
        params = self._to_new_params(
            noise_reused=True,
            outfile_snc=None,
            outfile_snl=None,
            oii_cat_in=None,
            oii_cat_out=None,
        )
        results = engine.run_etc_files(params)
        self._restore_sno2(results)

        elapsed_time = time.time() - start
        print(f"##### finished (elapsed_time: {elapsed_time:.1f}[sec]) #####")

        return 0

    def get_sno2(self):
        return self.sno2_zsps, self.sno2_sno2
