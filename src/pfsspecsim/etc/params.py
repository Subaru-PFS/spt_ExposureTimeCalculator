"""ETC input parameter model: ``EtcParams``, TOML loading, and ``MagSpec``.

See the "旧->新パラメータ名移行表" (old -> new parameter name migration
table) in the project plan for the mapping from the legacy stdin-driven C
engine's ALL_CAPS parameters to these snake_case fields.
"""

from __future__ import annotations

import dataclasses
import tomllib
from pathlib import Path
from typing import Any

import numpy as np

# Fields whose values are filesystem paths; used by `load_params` to coerce
# TOML/CLI string values (and by nothing else -- the dataclass itself always
# stores `Path | None`, never a bare string).
_PATH_FIELDS = frozenset(
    {
        "mag_file",
        "instr_config",
        "oii_cat_in",
        "oii_cat_out",
        "outdir",
        "outfile_noise",
        "outfile_snc",
        "outfile_snl",
        "outfile_oii",
    }
)

_VALID_SPECTROGRAPHS = frozenset({"ave", "sm1", "sm2", "sm3", "sm4"})


@dataclasses.dataclass
class EtcParams:
    """Resolved ETC input parameters (one flat, picklable, snake_case set).

    Defaults reproduce the legacy engine's defaults (see
    `pfsspecsim.pfsetc.Etc.__init__` / `scripts/run_etc.defaults`).
    """

    # Observing conditions
    seeing: float = 0.80  # arcsec FWHM @800nm
    zenith_ang: float = 35.0  # deg
    galactic_ext: float = 0.00  # E(B-V) mag
    field_ang: float = 0.45  # deg
    fiber_offset: float = 0.10  # arcsec (formerly hardcoded OFFSET_FIB)
    moon_zenith_ang: float = 30.0  # deg; >=90 disables moonlight (gsetc.c:956)
    moon_target_ang: float = 60.0  # deg
    moon_phase: float = 0.125  # 0=new, 0.5=full

    # Exposure
    exp_time: float = 900.0  # s / exposure
    exp_num: int = 4

    # Target
    mag: float | None = 22.5  # AB mag (flat); XOR with mag_file
    mag_file: Path | None = None  # 2-column wavelength[nm] mag file
    reff: float = 0.3  # arcsec
    line_flux: float = 1.0e-17  # erg/cm2/s
    line_width: float = 70.0  # km/s (sigma)

    # Instrument / model selection
    mr_mode: bool = False
    throughput_model: str = "20240714"
    spectrograph: str = "ave"  # 'ave' | 'sm1'..'sm4'
    instr_config: Path | None = None  # explicit override; takes precedence
    degrade: float = 1.0
    obsc_fov_dep: bool = True  # apply calc_obscuration correction to degrade
    sky_type: str = "11006"  # hex bitmask (formerly hardcoded SKYMODELS)
    hgcdte_sutr: bool = True  # formerly compile flag -DHGCDTE_SUTR

    # Systematics
    sky_sub_floor: float = 0.01
    diffuse_stray: float = 0.02

    # [OII] catalog
    oii_cat_in: Path | None = None
    oii_cat_out: Path | None = None
    min_snr: float = 9.0

    # Execution control / output
    # Thread cap per parallel stage (arm loops, output products); 1=serial;
    # result is unchanged regardless of value.
    n_workers: int = 3
    noise_reused: bool = False
    overwrite: bool = True
    outdir: Path = Path("out")
    outfile_noise: Path | None = Path("ref.noise.ecsv")  # relative to outdir
    outfile_snc: Path | None = Path("ref.snc.ecsv")
    outfile_snl: Path | None = Path("ref.snl.ecsv")
    outfile_oii: Path | None = None

    def validate(self) -> None:
        """Cross-field and range checks; raises ``ValueError`` on failure."""
        if (self.mag is None) == (self.mag_file is None):
            raise ValueError(
                "Exactly one of `mag` (flat AB magnitude) or `mag_file` "
                "(wavelength-dependent magnitude spectrum) must be set, got "
                f"mag={self.mag!r} mag_file={self.mag_file!r}"
            )
        if self.seeing <= 0:
            raise ValueError(f"seeing must be positive, got {self.seeing}")
        if not (0.0 <= self.zenith_ang < 90.0):
            raise ValueError(
                f"zenith_ang must be in [0, 90) deg, got {self.zenith_ang}"
            )
        if self.galactic_ext < 0:
            raise ValueError(f"galactic_ext must be >= 0, got {self.galactic_ext}")
        if self.field_ang < 0:
            raise ValueError(f"field_ang must be >= 0, got {self.field_ang}")
        if self.fiber_offset < 0:
            raise ValueError(f"fiber_offset must be >= 0, got {self.fiber_offset}")
        if not (0.0 <= self.moon_zenith_ang <= 180.0):
            raise ValueError(
                "moon_zenith_ang must be in [0, 180] deg, got "
                f"{self.moon_zenith_ang}"
            )
        if not (0.0 <= self.moon_target_ang <= 180.0):
            raise ValueError(
                "moon_target_ang must be in [0, 180] deg, got "
                f"{self.moon_target_ang}"
            )
        if not (0.0 <= self.moon_phase <= 1.0):
            raise ValueError(f"moon_phase must be in [0, 1], got {self.moon_phase}")
        if self.exp_time <= 0:
            raise ValueError(f"exp_time must be positive, got {self.exp_time}")
        if self.exp_num < 1:
            raise ValueError(f"exp_num must be >= 1, got {self.exp_num}")
        if self.reff < 0:
            raise ValueError(f"reff must be >= 0, got {self.reff}")
        if self.line_width <= 0:
            raise ValueError(f"line_width must be positive, got {self.line_width}")
        if self.degrade <= 0:
            raise ValueError(f"degrade must be positive, got {self.degrade}")
        if self.sky_sub_floor < 0:
            raise ValueError(f"sky_sub_floor must be >= 0, got {self.sky_sub_floor}")
        if self.diffuse_stray < 0:
            raise ValueError(f"diffuse_stray must be >= 0, got {self.diffuse_stray}")
        if self.min_snr <= 0:
            raise ValueError(f"min_snr must be positive, got {self.min_snr}")
        if self.n_workers < 1:
            raise ValueError(f"n_workers must be >= 1, got {self.n_workers}")
        if self.spectrograph.lower() not in _VALID_SPECTROGRAPHS:
            raise ValueError(
                "spectrograph must be one of "
                f"{sorted(_VALID_SPECTROGRAPHS)}, got {self.spectrograph!r}"
            )
        try:
            int(self.sky_type, 16)
        except ValueError as exc:
            raise ValueError(
                f"sky_type must be a hex string, got {self.sky_type!r}"
            ) from exc


def load_params(
    toml_path: str | Path | None = None,
    overrides: dict[str, Any] | None = None,
) -> EtcParams:
    """Build an `EtcParams`, merging defaults, an optional TOML file, and
    explicit overrides, with priority CLI (`overrides`) > TOML > defaults.

    `toml_path` is a flat TOML file with the same snake_case keys as
    `EtcParams` fields. Unknown keys, in either the TOML file or
    `overrides`, raise `ValueError`.
    """
    field_names = {f.name for f in dataclasses.fields(EtcParams)}
    values: dict[str, Any] = {}

    if toml_path is not None:
        with open(toml_path, "rb") as fh:
            toml_data = tomllib.load(fh)
        unknown = set(toml_data) - field_names
        if unknown:
            raise ValueError(
                f"Unknown key(s) in TOML file {toml_path}: {sorted(unknown)}"
            )
        values.update(toml_data)

    if overrides:
        unknown_overrides = set(overrides) - field_names
        if unknown_overrides:
            raise ValueError(f"Unknown override key(s): {sorted(unknown_overrides)}")
        values.update(overrides)

    # If exactly one of `mag`/`mag_file` was actually supplied (by the TOML
    # file and/or overrides) and the other was not mentioned at all, force
    # the unmentioned one to `None` -- otherwise it would silently fall back
    # to `EtcParams`'s own default (`mag=22.5`), spuriously colliding with
    # the one the caller *did* specify and tripping `validate()`'s XOR
    # check.
    if ("mag" in values) != ("mag_file" in values):
        if "mag" in values:
            values["mag_file"] = None
        else:
            values["mag"] = None

    for name in _PATH_FIELDS:
        if values.get(name) is not None:
            values[name] = Path(values[name])

    params = EtcParams(**values)
    params.validate()
    return params


def calc_obscuration(dist: float) -> tuple[float, float]:
    """Obscuration by the detector and surrounding structures.

    Verbatim port of `pfsspecsim.pfsetc.calc_obscuration` (see
    https://pfspipe.ipmu.jp/jira/browse/PIPE2D-980). `dist` is the distance
    from the field center on the telescope focal plane (deg, i.e. same
    quantity as `EtcParams.field_ang`).

    Returns
    -------
    obsc : float
        Fraction of light obscured at `dist`.
    corr : float
        Throughput correction factor relative to the obscuration assumed
        inside the ETC's own throughput model (independent of FoV location).
    """
    obsc_etc = 0.19  # obscuration assumed in ETC (independent of PFI location)
    obsc_inner = 0.27  # obscuration at center (PIPE2D-980)
    obsc_outer = 0.36  # obscuration at edge (PIPE2D-980)
    dist_edge = 0.675  # radius of FoV in deg (FIELD_ANG)
    obsc = (obsc_outer - obsc_inner) / dist_edge * dist + obsc_inner
    corr = (1 - obsc) / (1 - obsc_etc)
    return obsc, corr


def resolve_degrade(params: EtcParams) -> float:
    """Resolve the effective throughput `degrade` factor.

    When `params.obsc_fov_dep` is set, applies the `calc_obscuration`
    correction for `params.field_ang` on top of `params.degrade`; otherwise
    returns `params.degrade` unchanged.
    """
    if not params.obsc_fov_dep:
        return params.degrade
    _, corr = calc_obscuration(params.field_ang)
    return params.degrade * corr


class MagSpec:
    """A continuum AB-magnitude spectrum: a flat scalar or a tabulated file.

    Ports the file-padding logic of the legacy C reader (gsetc.c:1891-1952)
    and replaces both the C engine's sentinel-flagged (`mag == -99.9`)
    manual per-pixel interpolation (gsetc.c:1194-1210, 1401-1416) and the
    old Python wrapper's tmp-file resampling step (`pfsspecsim.pfsetc.Etc`,
    formerly pfsetc.py:217-233 -- a double linear interpolation that was a
    mathematical no-op) with a single `numpy.interp` call.
    """

    def __init__(
        self, mag: float | None = None, mag_file: str | Path | None = None
    ) -> None:
        if (mag is None) == (mag_file is None):
            raise ValueError(
                "Exactly one of `mag` or `mag_file` must be given, got "
                f"mag={mag!r} mag_file={mag_file!r}"
            )
        self._scalar_mag = mag
        if mag_file is not None:
            lam, mag_arr = np.loadtxt(mag_file, usecols=(0, 1), unpack=True)
            if np.any(np.diff(lam) <= 0):
                raise ValueError(
                    f"mag_file {mag_file!r} must have strictly increasing "
                    "wavelengths in column 0"
                )
            self._lam, self._mag = self._pad(lam, mag_arr)
        else:
            self._lam = None
            self._mag = None

    @staticmethod
    def _pad(lam: np.ndarray, mag: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """Pad the (lambda, mag) tables exactly as gsetc.c:1891-1952 does.

        Per-point sentinel: any `mag <= 0` in the *interior* (i.e. real
        input rows) is replaced with 99.9 (gsetc.c:1932), signalling
        "effectively no flux" downstream. The leading/trailing pad points
        (300nm / 1300nm boundary, or one grid step beyond the data if it
        already extends past those bounds) are a verbatim quirk: they reuse
        the *raw*, unthresholded first/last data value when extending the
        data rather than hard-clamping to 300/1300nm (gsetc.c:1918-1926,
        1940-1947) -- i.e. the `mag <= 0` -> 99.9 substitution is *not*
        re-applied to that reused endpoint value.
        """
        n = lam.size
        lam2 = np.empty(n + 2, dtype=float)
        mag2 = np.empty(n + 2, dtype=float)

        lam2[1:-1] = lam
        mag2[1:-1] = np.where(mag <= 0.0, 99.9, mag)

        if lam[0] > 300.0:
            lam2[0] = 300.0
            mag2[0] = 99.9
        else:
            lam2[0] = lam[0] - 1.0
            mag2[0] = mag[0]

        if lam[-1] < 1300.0:
            lam2[-1] = 1300.0
            mag2[-1] = 99.9
        else:
            lam2[-1] = lam[-1] + 1.0
            mag2[-1] = mag[-1]

        return lam2, mag2

    def __call__(self, lam_vec: np.ndarray) -> np.ndarray:
        """Evaluate the AB magnitude at the given wavelengths (nm)."""
        lam_vec = np.asarray(lam_vec, dtype=float)
        if self._lam is None:
            return np.full(lam_vec.shape, self._scalar_mag, dtype=float)
        return np.interp(lam_vec, self._lam, self._mag)
