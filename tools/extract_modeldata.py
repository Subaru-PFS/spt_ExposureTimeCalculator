#!/usr/bin/env python
"""Extract numeric lookup tables baked into the legacy C ETC sources.

This script regex-parses the C source files ``legacy/c_src/gsetc.c`` and
``legacy/c_src/modeldata.h`` for the hard-coded numeric arrays that the C
engine uses as lookup tables (sky spectral atlases, atmospheric transmission
grids, silicon optical constants, and the Milky Way dust extinction law), and
serializes them to a single compressed ``.npz`` archive plus a provenance
``README.md``.

It was originally run once, before the C sources were removed from the
``src/`` tree (see task T15 of the pure-Python port plan); the C sources are
now vendored permanently, frozen, at ``legacy/c_src/`` (see
``legacy/c_src/README.md``). The script itself is committed to document how
the archive was produced and to allow re-extraction for validation.

Usage::

    uv run python tools/extract_modeldata.py

Outputs (relative to the repository root)::

    src/pfsspecsim/etc/data/modeldata.npz
    src/pfsspecsim/etc/data/README.md
"""

from __future__ import annotations

import argparse
import re
import subprocess
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parent.parent
GSETC_C = REPO_ROOT / "legacy" / "c_src" / "gsetc.c"
MODELDATA_H = REPO_ROOT / "legacy" / "c_src" / "modeldata.h"
OUTPUT_DIR = REPO_ROOT / "src" / "pfsspecsim" / "etc" / "data"

# Number of entries in the OH night-sky line table (modeldata.h:5653,
# `#define N_IR_OH_LINE 698`).
N_IR_OH_LINE = 698


def _extract_c_array(
    text: str, name: str, ctype: str = r"(?:double|long)"
) -> np.ndarray:
    """Pull the numeric literals out of a C array initializer.

    Matches declarations of the form ``<ctype> <name>[] = { ... };`` (the
    body may span many lines) and returns the comma-separated numbers as a
    1-D float64 array. Works for both file-scope tables (modeldata.h) and
    tables declared as local variables inside a function body (gsetc.c),
    since the regex only anchors on the declaration/initializer syntax.
    """
    pattern = re.compile(
        rf"\b{ctype}\s+{re.escape(name)}\s*\[\s*\]\s*=\s*\{{(.*?)\}}\s*;",
        re.DOTALL,
    )
    match = pattern.search(text)
    if match is None:
        raise ValueError(f"Could not find array '{name}' in source text")
    body = match.group(1)
    tokens = [tok.strip() for tok in body.split(",")]
    values = [float(tok) for tok in tokens if tok]
    return np.array(values, dtype=np.float64)


def extract_uves_lambda(modeldata_text: str) -> np.ndarray:
    """UVES sky atlas wavelengths (air, nm). modeldata.h:12."""
    arr = _extract_c_array(modeldata_text, "gsSKY_UVES_LAMBDA")
    if arr.shape != (2816,):
        raise ValueError(f"gsSKY_UVES_LAMBDA: expected shape (2816,), got {arr.shape}")
    return arr


def extract_uves_int(modeldata_text: str) -> np.ndarray:
    """UVES sky atlas line intensities (1e-16 erg/cm^2/s/arcsec^2 == 1e-12
    erg/m^2/s/arcsec^2). modeldata.h:2830."""
    arr = _extract_c_array(modeldata_text, "gsSKY_UVES_INT")
    if arr.shape != (2816,):
        raise ValueError(f"gsSKY_UVES_INT: expected shape (2816,), got {arr.shape}")
    return arr


def extract_oh_data(modeldata_text: str) -> np.ndarray:
    """OH night-sky line table, reshaped to (N_IR_OH_LINE, 2) columns of
    [vacuum nm, intensity]. modeldata.h:5655."""
    flat = _extract_c_array(modeldata_text, "OHDATA")
    if flat.size != N_IR_OH_LINE * 2:
        raise ValueError(f"OHDATA: expected {N_IR_OH_LINE * 2} values, got {flat.size}")
    arr = flat.reshape(N_IR_OH_LINE, 2)
    if arr.shape != (698, 2):
        raise ValueError(f"oh_data: expected shape (698, 2), got {arr.shape}")
    return arr


def extract_atm_trans_kp(modeldata_text: str) -> np.ndarray:
    """Kitt Peak atmospheric transmission spectrum.

    Grid: 500 nm + 0.025 nm * i (used at gsetc.c:444, index clamp at
    gsetc.c:447 requires xint+1 <= 39999, i.e. len >= 40000). modeldata.h:6356.
    """
    arr = _extract_c_array(modeldata_text, "AtmTransKP")
    if arr.shape[0] < 40000:
        raise ValueError(f"AtmTransKP: expected len >= 40000, got {arr.shape[0]}")
    return arr


def extract_mk_trans_3mm(modeldata_text: str) -> np.ndarray:
    """Gemini/Mauna Kea simulated transmission spectrum at 3mm precipitable
    water vapor.

    Grid: 900 nm + 0.02 nm * i (used at gsetc.c:463, index clamp at
    gsetc.c:467 requires xint+1 <= 30000, i.e. len >= 30001). modeldata.h:8360.
    """
    arr = _extract_c_array(modeldata_text, "MKTrans_3mm")
    if arr.shape[0] < 30001:
        raise ValueError(f"MKTrans_3mm: expected len >= 30001, got {arr.shape[0]}")
    return arr


def extract_dust_norm(gsetc_text: str) -> np.ndarray:
    """Milky Way dust extinction law A_lambda/E(B-V), inline table `norm[]`
    in gsGalactic_Alambda__EBV(). gsetc.c:347-368.

    Log-spaced grid: index = 100/ln(10) * ln(10/lambda_um), lambda in
    0.1-10 um, clamped to [0, 199] before linear interpolation against
    norm[il] and norm[il+1] (gsetc.c:378-384), hence 201 entries.
    """
    arr = _extract_c_array(gsetc_text, "norm")
    if arr.shape != (201,):
        raise ValueError(f"dust_norm: expected shape (201,), got {arr.shape}")
    return arr


def extract_si_index_table(gsetc_text: str) -> np.ndarray:
    """Silicon complex refractive index (real part) table `si_table[]` in
    gsOP_Si_indexreal(). gsetc.c:295-305.

    Linear grid: 200 nm + 5 nm * i, over lambda in [200, 1100] nm
    (gsetc.c:311-321), 181 entries.
    """
    arr = _extract_c_array(gsetc_text, "si_table")
    if arr.shape != (181,):
        raise ValueError(f"si_index_table: expected shape (181,), got {arr.shape}")
    return arr


def get_git_sha(repo_root: Path) -> str:
    try:
        result = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=repo_root,
            capture_output=True,
            text=True,
            check=True,
        )
        return result.stdout.strip()
    except (subprocess.CalledProcessError, FileNotFoundError):
        return "unknown"


README_TEMPLATE = """\
# modeldata.npz -- extracted lookup tables

This archive was generated by `tools/extract_modeldata.py`, which
regex-parses the numeric lookup tables hard-coded into the legacy C ETC
engine (`legacy/c_src/gsetc.c` and `legacy/c_src/modeldata.h`) and
serializes them with `np.savez_compressed(..., allow_pickle=False)`. It is
committed as a snapshot so the pure-Python engine can load the tables
without depending on the C sources at runtime.

Source files:

- `legacy/c_src/gsetc.c`
- `legacy/c_src/modeldata.h`

Extraction git SHA (repository HEAD at the time this archive was produced):
`{git_sha}`

## Arrays

| npz key | shape | dtype | source | units / notes |
|---|---|---|---|---|
| `uves_lambda` | {uves_lambda_shape} | float64 | `gsSKY_UVES_LAMBDA`, modeldata.h:12 | air nm |
| `uves_int` | {uves_int_shape} | float64 | `gsSKY_UVES_INT`, modeldata.h:2830 | 1e-12 erg/m^2/s/arcsec^2 (== 1e-16 erg/cm^2/s/arcsec^2 as documented in modeldata.h) |
| `oh_data` | {oh_data_shape} | float64 | `OHDATA`, modeldata.h:5655 | columns: [vacuum nm, 1e-16 erg/cm^2/s/arcsec^2], normalized to 14.8 mag Vega/arcsec^2 in J band |
| `atm_trans_kp` | {atm_trans_kp_shape} | float64 | `AtmTransKP`, modeldata.h:6356 | dimensionless transmission; grid = 500 nm + 0.025 nm * i (gsetc.c:444) |
| `mk_trans_3mm` | {mk_trans_3mm_shape} | float64 | `MKTrans_3mm`, modeldata.h:8360 | dimensionless transmission; grid = 900 nm + 0.02 nm * i (gsetc.c:463) |
| `dust_norm` | {dust_norm_shape} | float64 | inline `norm[201]`, gsetc.c:347-368 | A_lambda / E(B-V); log-spaced grid, index = 100/ln(10) * ln(10/lambda_um), lambda in [0.1, 10] um |
| `si_index_table` | {si_index_table_shape} | float64 | inline `si_table[181]`, gsetc.c:295-305 | real part of Si refractive index; linear grid = 200 nm + 5 nm * i, lambda in [200, 1100] nm |

## Upstream references

- UVES sky spectral line atlas: ESO UVES sky atlas (modeldata.h:3-8).
- Kitt Peak atmospheric transmission spectrum (`AtmTransKP`): Kitt Peak
  Observatory transmission spectrum (gsetc.c:443).
- ATRAN-modeled 3mm precipitable water transmission (`MKTrans_3mm`): Gemini
  Observatory simulated atmospheric transmission spectrum for Mauna Kea
  (gsetc.c:459-461).
- Milky Way dust extinction law (`dust_norm`):
  ftp://ftp.astro.princeton.edu/draine/dust/mix/kext_albedo_WD_MW_3.1_60_D03.all
  (retrieved 2011-03-18; gsetc.c:340).
- Silicon complex refractive index (`si_index_table`):
  http://snap.lbl.gov/ccdweb/ccd_data/complex_index.dat (gsetc.c:282); see
  also Don Groom, SPIE 1999, and D. F. Edwards in the Handbook of Optical
  Constants of Solids (1985) for lambda < 750 nm (gsetc.c:284-288).

## Regenerating

```
uv run python tools/extract_modeldata.py
```
"""


def build_readme(git_sha: str, arrays: dict[str, np.ndarray]) -> str:
    return README_TEMPLATE.format(
        git_sha=git_sha,
        uves_lambda_shape=arrays["uves_lambda"].shape,
        uves_int_shape=arrays["uves_int"].shape,
        oh_data_shape=arrays["oh_data"].shape,
        atm_trans_kp_shape=arrays["atm_trans_kp"].shape,
        mk_trans_3mm_shape=arrays["mk_trans_3mm"].shape,
        dust_norm_shape=arrays["dust_norm"].shape,
        si_index_table_shape=arrays["si_index_table"].shape,
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=OUTPUT_DIR,
        help="Directory to write modeldata.npz and README.md into "
        f"(default: {OUTPUT_DIR.relative_to(REPO_ROOT)})",
    )
    args = parser.parse_args()

    gsetc_text = GSETC_C.read_text()
    modeldata_text = MODELDATA_H.read_text()

    arrays = {
        "uves_lambda": extract_uves_lambda(modeldata_text),
        "uves_int": extract_uves_int(modeldata_text),
        "oh_data": extract_oh_data(modeldata_text),
        "atm_trans_kp": extract_atm_trans_kp(modeldata_text),
        "mk_trans_3mm": extract_mk_trans_3mm(modeldata_text),
        "dust_norm": extract_dust_norm(gsetc_text),
        "si_index_table": extract_si_index_table(gsetc_text),
    }

    # Spot-check known values (guards against a regex/parsing regression).
    assert arrays["dust_norm"][0] == 0.24174, arrays["dust_norm"][0]
    assert arrays["dust_norm"][200] == 13.92363, arrays["dust_norm"][200]
    assert arrays["si_index_table"][60] == 4.2975, arrays["si_index_table"][60]

    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    npz_path = output_dir / "modeldata.npz"
    np.savez_compressed(npz_path, allow_pickle=False, **arrays)
    print(f"Wrote {npz_path}")

    git_sha = get_git_sha(REPO_ROOT)
    readme_path = output_dir / "README.md"
    readme_path.write_text(build_readme(git_sha, arrays))
    print(f"Wrote {readme_path}")


if __name__ == "__main__":
    main()
