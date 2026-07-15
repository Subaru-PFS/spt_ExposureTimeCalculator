"""Instrument configuration file parser -> `Spectrograph` dataclass.

Ports the `PFS.*.dat` spectrograph-attribute file reader
(``gsReadSpectrographConfig``, gsetc.c:1477-1686) plus two small pieces of
supporting logic that the C engine repeats verbatim at several call sites:

* :func:`spectro_arm` <- ``spectro_arm`` (gsetc.c:76-83): maps the internal
  0-based arm loop index ``ia`` to the *output* arm id, which differs only
  in medium-resolution (MR) mode, where the red arm (``ia == 1``) reports as
  arm id 3 instead of 1.
* :func:`field_interp` <- the 5-node field-angle interpolation that gsetc.c
  open-codes three times (gsetc.c:517-533, 584-594, 768-779) for
  ``rms_spot``/``EFL``, ``vignette``, and (again) ``EFL``. Implemented once
  here; callers pass whichever of the 5-element attribute arrays they need.

:func:`find_config_file` reproduces the file-naming convention of the old
subprocess wrapper (``pfsspecsim.pfsetc.Etc.run``, formerly pfsetc.py:186-199),
resolving packaged config files under ``pfsspecsim/config/`` via
``importlib.resources`` instead of a `HOME_DIR`-relative path.
"""

from __future__ import annotations

import dataclasses
from importlib import resources
from pathlib import Path

import numpy as np

from .constants import ADJUST_THROUGHPUT_LR, ADJUST_THROUGHPUT_MR

#: Number of field-angle sample nodes in the 5-element attribute arrays
#: (rms_spot, EFL, vignette): gsetc.c's SPECTRO_ATTRIB struct fixes this at 5.
_N_FIELD_NODES = 5

#: Maximum number of spectrograph arms (gsetc.c:26, `#define MAXARM 3`).
MAXARM = 3


@dataclasses.dataclass
class Spectrograph:
    """Parsed spectrograph attributes; mirrors C's `SPECTRO_ATTRIB` struct.

    Field-angle-dependent quantities (`rms_spot`, `EFL`, `vignette`) keep the
    C struct's fixed-5-node layout (field angle = 0, rfov/4, rfov/2,
    3*rfov/4, rfov); see :func:`field_interp`. Per-arm quantities are
    length-`N_arms` arrays indexed by the internal arm index `ia`
    (0-based; see :func:`spectro_arm` for the output arm id mapping).
    Throughput-grid quantities (`l`, `T`) are concatenated across all arms,
    with `istart[ia]:istart[ia + 1]` giving the slice for arm `ia`.

    A plain, picklable dataclass; all array fields are `numpy.ndarray`.
    """

    # --- OPTICS / SPOT / VIGNET / FIBER (field-angle-dependent, 5 nodes) ---
    D_outer: float  # Outer diameter, meters
    centobs: float  # Central obscuration (fraction)
    rfov: float  # Field of view radius, degrees
    EFL: np.ndarray  # Effective focal length, meters; shape (5,)
    rms_spot: np.ndarray  # RMS spot size, microns; shape (5,)
    vignette: np.ndarray  # Vignetting factor; shape (5,)

    # --- FIBER ---
    fiber_ent_rad: float  # Fiber entrance radius, microns

    # --- ARMS / MEDIUM_RESOLUTION ---
    N_arms: int
    MR: bool

    # --- PARAM (per arm) ---
    lmin: np.ndarray  # Min wavelength, nm; shape (N_arms,)
    lmax: np.ndarray  # Max wavelength, nm; shape (N_arms,)
    npix: np.ndarray  # Number of pixels; shape (N_arms,), int
    dl: np.ndarray  # Wavelength spacing, nm/pix, derived; shape (N_arms,)
    width: np.ndarray  # Trace width used in analysis, pixels; shape (N_arms,), int

    # --- CAMERA (per arm) ---
    fratio: np.ndarray  # Camera f/ratio; shape (N_arms,)
    thick: np.ndarray  # Detector thickness, microns; shape (N_arms,)
    pix: np.ndarray  # Pixel scale, microns; shape (N_arms,)
    temperature: np.ndarray  # Camera temperature, K; shape (N_arms,)
    rms_cam: np.ndarray  # Camera RMS spot per axis, microns; shape (N_arms,)
    diam: np.ndarray  # Geometric fiber image diameter, microns; shape (N_arms,)
    dark: np.ndarray  # Dark current, e/pix/s; shape (N_arms,)
    read: np.ndarray  # Read noise, e rms; shape (N_arms,)
    sep: np.ndarray  # Trace spacing, microns; shape (N_arms,)

    # --- NLINES (per arm) ---
    nline: np.ndarray  # Effective number of lines; shape (N_arms,)

    # --- HGCDTE (per arm) ---
    Dtype: np.ndarray  # 0 = Si, 1 = HgCdTe; shape (N_arms,), int

    # --- THRPUT (concatenated grid across arms) ---
    N_thr: int
    istart: np.ndarray  # Arm boundaries into l/T; shape (N_arms + 1,), int
    l: np.ndarray  # Wavelength at throughput grid points, nm; shape (N_thr,)
    T: np.ndarray  # Throughput at grid points (excl. atmosphere & fiber geom)

    # --- Set later by the caller (not read from the config file) ---
    sysfrac: float = 0.0  # Sky-subtraction systematic (gsetc.c:1684)
    diffuse_stray: float = 0.0


def spectro_arm(ia: int, mr: bool) -> int:
    """Map the internal 0-based arm loop index to the output arm id.

    Verbatim port of `spectro_arm` (gsetc.c:76-83): in low-resolution (LR)
    mode the output id equals `ia`; in medium-resolution (MR) mode the red
    arm (`ia == 1`) reports as arm id 3 (all other arms unchanged).
    """
    if not mr:
        return ia
    return 3 if ia == 1 else ia


def field_interp(values5: np.ndarray, fieldang: float, rfov: float) -> float:
    """Interpolate a 5-node field-angle-dependent attribute at `fieldang`.

    `values5` holds samples at field angle 0, rfov/4, rfov/2, 3*rfov/4, rfov
    (degrees). Ports the piecewise-linear interpolation that gsetc.c
    open-codes three times: gsetc.c:517-533 (`rms_spot`/`EFL` in
    `gsGeometricThroughput`), gsetc.c:584-594 (`vignette` in `gsAeff`), and
    gsetc.c:768-779 (`EFL` again, in `gsGetNoise`). Below node 0 or above
    node 4 the value is clamped to the corresponding endpoint (the C code
    special-cases `i<0` and `i>=4` this way rather than extrapolating).
    """
    x = 4.0 * fieldang / rfov
    i = int(np.floor(x))
    if i < 0:
        return float(values5[0])
    if i >= _N_FIELD_NODES - 1:
        return float(values5[_N_FIELD_NODES - 1])
    frac = x - i
    return float(values5[i] + (values5[i + 1] - values5[i]) * frac)


def load_spectrograph_config(path: str | Path, degrade: float) -> Spectrograph:
    """Parse a `PFS.*.dat` spectrograph configuration file.

    Port of `gsReadSpectrographConfig` (gsetc.c:1477-1686). `degrade` is a
    throughput scale factor folded into `T` at read time (gsetc.c:1608-1613),
    together with the empirical `ADJUST_THROUGHPUT_LR`/`ADJUST_THROUGHPUT_MR`
    per-arm factors (indexed by the internal arm index `ia`, *not* the output
    arm id -- see the "index 規約" note in the task brief).

    Section keywords recognized (in any order, one or more lines each,
    except THRPUT, which owns all subsequent lines until as many `D` marker
    lines as there are arms have been seen):
    OPTICS, SPOT, VIGNET, FIBER, ARMS, MEDIUM_RESOLUTION, PARAM (one line per
    arm), CAMERA (one line per arm), THRPUT, HGCDTE (one line per HgCdTe
    arm), NLINES (one line per arm).

    Raises `ValueError` on missing/malformed sections (the C code prints to
    stderr and calls `exit(1)`; ints from `MAXARM`/`MAXPIX`/`MAXNTHR` static
    bounds are not enforced here since numpy arrays are sized dynamically).
    """
    path = Path(path)
    lines = path.read_text().splitlines()

    D_outer = -1.0
    centobs = 0.0
    rfov = 0.0
    EFL = np.full(_N_FIELD_NODES, np.nan)
    rms_spot = np.full(_N_FIELD_NODES, -1.0)
    vignette = np.ones(_N_FIELD_NODES)
    fiber_ent_rad = -1.0
    N_arms = 0
    MR = False

    lmin = lmax = npix = dl = width = None
    fratio = thick = pix = temperature = rms_cam = diam = dark = read = sep = None
    nline = None
    Dtype = None
    N_thr = 0
    istart = None
    thr_l: list[float] = []
    thr_T: list[float] = []

    idx = 0
    n_lines = len(lines)
    while idx < n_lines:
        raw = lines[idx]
        idx += 1
        if raw.startswith("OPTICS"):
            vals = [float(v) for v in raw[7:].split()]
            if len(vals) != 8:
                raise ValueError(
                    f"{path}: OPTICS line has {len(vals)}/8 arguments: {raw!r}"
                )
            D_outer, centobs, rfov = vals[0], vals[1], vals[2]
            EFL = np.array(vals[3:8], dtype=float)
        elif raw.startswith("SPOT"):
            vals = [float(v) for v in raw[5:].split()]
            if len(vals) != 5:
                raise ValueError(
                    f"{path}: SPOT line has {len(vals)}/5 arguments: {raw!r}"
                )
            rms_spot = np.array(vals, dtype=float)
        elif raw.startswith("VIGNET"):
            vals = [float(v) for v in raw[7:].split()]
            if len(vals) != 5:
                raise ValueError(
                    f"{path}: VIGNET line has {len(vals)}/5 arguments: {raw!r}"
                )
            vignette = np.array(vals, dtype=float)
        elif raw.startswith("FIBER"):
            vals = [float(v) for v in raw[6:].split()]
            if len(vals) != 1:
                raise ValueError(
                    f"{path}: FIBER line has {len(vals)}/1 arguments: {raw!r}"
                )
            fiber_ent_rad = vals[0]
        elif raw.startswith("MEDIUM_RESOLUTION"):
            tok = raw[17:].split()
            if len(tok) != 1 or tok[0] not in ("0", "1"):
                raise ValueError(f"{path}: illegal MEDIUM_RESOLUTION line: {raw!r}")
            MR = bool(int(tok[0]))
        elif raw.startswith("ARMS"):
            tok = raw[5:].split()
            if len(tok) != 1:
                raise ValueError(f"{path}: illegal ARMS line: {raw!r}")
            N_arms = int(tok[0])
            lmin = np.full(N_arms, -1.0)
            lmax = np.full(N_arms, -1.0)
            npix = np.zeros(N_arms, dtype=np.int64)
            dl = np.zeros(N_arms)
            width = np.zeros(N_arms, dtype=np.int64)
            fratio = np.full(N_arms, -1.0)
            thick = np.zeros(N_arms)
            pix = np.zeros(N_arms)
            temperature = np.zeros(N_arms)
            rms_cam = np.zeros(N_arms)
            diam = np.zeros(N_arms)
            dark = np.zeros(N_arms)
            read = np.zeros(N_arms)
            sep = np.zeros(N_arms)
            nline = np.full(N_arms, 1e12)
            Dtype = np.zeros(N_arms, dtype=np.int64)
        elif raw.startswith("PARAM"):
            if N_arms <= 0:
                raise ValueError(f"{path}: PARAM line before ARMS: {raw!r}")
            tok = raw[6:].split()
            if len(tok) != 5:
                raise ValueError(
                    f"{path}: PARAM line has {len(tok)}/5 arguments: {raw!r}"
                )
            i = int(tok[0])
            if not (0 <= i < N_arms):
                raise ValueError(f"{path}: illegal PARAM arm #{i}: {raw!r}")
            a_lmin, a_lmax, a_npix, a_width = (
                float(tok[1]),
                float(tok[2]),
                int(tok[3]),
                int(tok[4]),
            )
            lmin[i] = a_lmin
            lmax[i] = a_lmax
            npix[i] = a_npix
            dl[i] = (a_lmax - a_lmin) / a_npix
            width[i] = a_width
        elif raw.startswith("CAMERA"):
            if N_arms <= 0:
                raise ValueError(f"{path}: CAMERA line before ARMS: {raw!r}")
            tok = raw[7:].split()
            if len(tok) != 10:
                raise ValueError(
                    f"{path}: CAMERA line has {len(tok)}/10 arguments: {raw!r}"
                )
            i = int(tok[0])
            if not (0 <= i < N_arms):
                raise ValueError(f"{path}: illegal CAMERA arm #{i}: {raw!r}")
            (
                fratio[i],
                thick[i],
                pix[i],
                temperature[i],
                rms_cam[i],
                diam[i],
                dark[i],
                read[i],
                sep[i],
            ) = (float(v) for v in tok[1:10])
        elif raw.startswith("HGCDTE"):
            if N_arms <= 0:
                raise ValueError(f"{path}: HGCDTE line before ARMS: {raw!r}")
            i = int(raw[7:].split()[0])
            if not (0 <= i < N_arms):
                raise ValueError(f"{path}: HGCDTE {i}: illegal arm index")
            Dtype[i] = 1
        elif raw.startswith("NLINES"):
            if N_arms <= 0:
                raise ValueError(f"{path}: NLINES line before ARMS: {raw!r}")
            tok = raw[7:].split()
            i, n = int(tok[0]), float(tok[1])
            if not (0 <= i < N_arms):
                raise ValueError(f"{path}: NLINES {i}: illegal arm index")
            if n < 10:
                raise ValueError(f"{path}: NLINES {i}: nlines={n} is illegal.")
            nline[i] = n
        elif raw.startswith("THRPUT"):
            if N_arms <= 0:
                raise ValueError(f"{path}: THRPUT section before ARMS: {raw!r}")
            i = 0
            i_arm = 0
            istart = np.zeros(N_arms + 1, dtype=np.int64)
            adjust = ADJUST_THROUGHPUT_MR if MR else ADJUST_THROUGHPUT_LR
            while True:
                if idx >= n_lines:
                    raise ValueError(f"{path}: unexpected EOF at THRPUT grid point {i}")
                thr_line = lines[idx]
                idx += 1
                if thr_line.startswith("D"):
                    i_arm += 1
                    istart[i_arm] = i
                    if i_arm == N_arms:
                        N_thr = i
                        break
                    continue
                tok = thr_line.split()
                if len(tok) < 6:
                    raise ValueError(
                        f"{path}: illegal throughput table line: "
                        f"{len(tok)}/6 arguments: {thr_line!r}"
                    )
                lam = float(tok[0])
                t1, t2, t3, t4, t5 = (float(v) for v in tok[1:6])
                thr_l.append(lam)
                thr_T.append(t1 * t2 * t3 * t4 * t5 * degrade * adjust[i_arm])
                i += 1

    if D_outer <= 0:
        raise ValueError(f"{path}: illegal outer diameter or no OPTICS line.")
    if np.any(rms_spot < 0):
        raise ValueError(f"{path}: illegal rms spot or no SPOT line.")
    if fiber_ent_rad <= 0:
        raise ValueError(f"{path}: illegal fiber radius or no FIBER line.")
    if N_arms <= 0:
        raise ValueError(f"{path}: {N_arms} arms illegal or no ARMS line.")
    if lmin is None or np.any(lmin < 0):
        raise ValueError(f"{path}: illegal lambda min or missing PARAM line(s).")
    if istart is None:
        raise ValueError(f"{path}: missing THRPUT section.")

    return Spectrograph(
        D_outer=D_outer,
        centobs=centobs,
        rfov=rfov,
        EFL=EFL,
        rms_spot=rms_spot,
        vignette=vignette,
        fiber_ent_rad=fiber_ent_rad,
        N_arms=N_arms,
        MR=MR,
        lmin=lmin,
        lmax=lmax,
        npix=npix,
        dl=dl,
        width=width,
        fratio=fratio,
        thick=thick,
        pix=pix,
        temperature=temperature,
        rms_cam=rms_cam,
        diam=diam,
        dark=dark,
        read=read,
        sep=sep,
        nline=nline,
        Dtype=Dtype,
        N_thr=N_thr,
        istart=istart,
        l=np.array(thr_l, dtype=float),
        T=np.array(thr_T, dtype=float),
        sysfrac=0.0,
        diffuse_stray=0.0,
    )


def find_config_file(
    throughput_model: str, spectrograph: str = "ave", mr_mode: bool = False
) -> Path:
    """Resolve the packaged `PFS.*.dat` config file path for the given model.

    Port of the filename convention in the old subprocess wrapper
    (`pfsspecsim.pfsetc.Etc.run`, formerly pfsetc.py:186-199): for
    `spectrograph in {sm1, sm2, sm3, sm4}` the file is
    `PFS.<throughput_model>.<spectrograph>.dat` (`PFS.redMR.<...>` in MR
    mode); for the average/nominal spectrograph (`spectrograph == "ave"`,
    formerly any other string) it is `PFS.<throughput_model>.dat`
    (`PFS.redMR.<throughput_model>.dat` in MR mode). Files are packaged as
    `pfsspecsim/config/*.dat` and resolved via `importlib.resources` rather
    than a `HOME_DIR`-relative path.

    Raises `FileNotFoundError` if the resolved file does not exist in the
    installed package.
    """
    spectrograph = spectrograph.lower()
    prefix = "PFS.redMR" if mr_mode else "PFS"
    if spectrograph in ("sm1", "sm2", "sm3", "sm4"):
        filename = f"{prefix}.{throughput_model}.{spectrograph}.dat"
    else:
        filename = f"{prefix}.{throughput_model}.dat"

    # `resources.as_file` would, for a zipped/namespace package, extract to a
    # temp file that is removed once its context manager exits -- unsafe to
    # return past that point. The project pins `zip-safe = false`
    # (pyproject.toml), so the package is always installed as a real
    # directory tree and `as_file` yields the on-disk path with no
    # extraction/cleanup step, making it safe to return here.
    config_traversable = resources.files("pfsspecsim") / "config" / filename
    with resources.as_file(config_traversable) as config_path:
        if not config_path.is_file():
            raise FileNotFoundError(
                f"No packaged spectrograph config file found for "
                f"throughput_model={throughput_model!r}, "
                f"spectrograph={spectrograph!r}, mr_mode={mr_mode!r} "
                f"(looked for {filename!r} under pfsspecsim/config/)"
            )
        return config_path
