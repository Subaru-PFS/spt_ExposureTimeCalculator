#!/usr/bin/env python3
"""Custom build backend for pfsspecsim - compiles C executables via Makefile."""

import subprocess
import sys
from pathlib import Path
from setuptools import build_meta as _orig

__all__ = [
    "get_requires_for_build_wheel",
    "get_requires_for_build_sdist",
    "get_requires_for_build_editable",
    "prepare_metadata_for_build_wheel",
    "prepare_metadata_for_build_editable",
    "build_wheel",
    "build_sdist",
    "build_editable",
]

get_requires_for_build_wheel = _orig.get_requires_for_build_wheel
get_requires_for_build_sdist = _orig.get_requires_for_build_sdist
get_requires_for_build_editable = _orig.get_requires_for_build_editable
prepare_metadata_for_build_wheel = _orig.prepare_metadata_for_build_wheel
prepare_metadata_for_build_editable = _orig.prepare_metadata_for_build_editable


def _compile_executables():
    """Compile gsetc and gsetc_omp executables using the Makefile."""
    root_dir = Path(__file__).parent.absolute()

    print("=" * 60)
    print("Building C executables (gsetc.x and gsetc_omp.x)")
    print("=" * 60)

    cmd = ["make", "-f", "Makefile", "all"]

    try:
        result = subprocess.run(
            cmd, cwd=root_dir, check=True, capture_output=True, text=True
        )
        print(result.stdout)
        if result.stderr:
            print(result.stderr, file=sys.stderr)
    except subprocess.CalledProcessError as e:
        print(f"Build failed: {e}", file=sys.stderr)
        print(e.stdout, file=sys.stderr)
        print(e.stderr, file=sys.stderr)
        raise RuntimeError(f"Failed to compile C executables: {e}")

    # Verify executables were created
    bin_dir = root_dir / "python" / "pfsspecsim" / "bin"
    gsetc_serial = bin_dir / "gsetc.x"
    gsetc_omp = bin_dir / "gsetc_omp.x"

    if not gsetc_serial.exists():
        raise RuntimeError(f"Failed to create {gsetc_serial}")
    if not gsetc_omp.exists():
        raise RuntimeError(f"Failed to create {gsetc_omp}")

    print(f"Successfully built: {gsetc_serial}")
    print(f"Successfully built: {gsetc_omp}")
    print("=" * 60)


def build_wheel(wheel_directory, config_settings=None, metadata_directory=None):
    """Build wheel with custom compilation step."""
    _compile_executables()
    return _orig.build_wheel(
        wheel_directory,
        config_settings=config_settings,
        metadata_directory=metadata_directory
    )


def build_sdist(sdist_directory, config_settings=None):
    """Build source distribution (includes source, not compiled binaries)."""
    return _orig.build_sdist(sdist_directory, config_settings=config_settings)


def build_editable(wheel_directory, config_settings=None, metadata_directory=None):
    """Build editable wheel with custom compilation step."""
    _compile_executables()
    return _orig.build_editable(
        wheel_directory,
        config_settings=config_settings,
        metadata_directory=metadata_directory
    )
