"""Tests for pfsspecsim.etc.config (task T4).

Port of `gsReadSpectrographConfig` (gsetc.c:1477-1686), `spectro_arm`
(gsetc.c:76-83), and the 5-node field-angle interpolation that gsetc.c
open-codes three times (gsetc.c:517-533, 584-594, 768-779).
"""

from pathlib import Path

import numpy as np
import pytest

from pfsspecsim.etc.config import (
    find_config_file,
    field_interp,
    load_spectrograph_config,
    spectro_arm,
)

# Protected fixture; never modify (see task brief).
PROTECTED_FIXTURE = (
    Path(__file__).resolve().parents[2] / "tests" / "PFS.20211220.dat"
)


class TestLoadSpectrographConfig:
    def test_lr_config_arm_attributes(self):
        path = find_config_file("20240714")
        sp = load_spectrograph_config(path, degrade=1.0)

        assert sp.N_arms == 3
        assert sp.MR is False
        np.testing.assert_array_equal(sp.Dtype, [0, 0, 1])
        np.testing.assert_allclose(sp.lmin, [380.0, 630.0, 940.0])
        np.testing.assert_array_equal(sp.npix, [4096, 4096, 4096])
        np.testing.assert_allclose(sp.nline, [6.73e4, 5.57e4, 5.91e4])

    def test_lr_config_dl_derived(self):
        # dl = (lmax - lmin) / npix (gsetc.c:1562).
        path = find_config_file("20240714")
        sp = load_spectrograph_config(path, degrade=1.0)
        expected_dl = (sp.lmax - sp.lmin) / sp.npix
        np.testing.assert_allclose(sp.dl, expected_dl)

    @pytest.mark.parametrize("degrade", [1.0, 0.5, 1.37])
    def test_lr_config_first_throughput_value(self, degrade):
        # T = prod(cols 1..5) * degrade * ADJUST_THROUGHPUT_LR[ia]
        # (gsetc.c:1608-1613); blue arm (ia=0) first grid point is
        # lambda=380.0, cols = 0.0302 1.0 1.0 1.0 1.0 (a 6th, unused, column
        # follows). ADJUST_THROUGHPUT_LR[0] = 1.08 (gsetc.c:1458).
        path = find_config_file("20240714")
        sp = load_spectrograph_config(path, degrade=degrade)
        assert sp.l[0] == pytest.approx(380.0)
        assert sp.T[0] == pytest.approx(0.0302 * 1.08 * degrade, rel=1e-12)

    def test_mr_config_flag_and_arm_ids(self):
        path = find_config_file("20240714", mr_mode=True)
        sp = load_spectrograph_config(path, degrade=1.0)

        assert sp.MR is True
        assert [spectro_arm(ia, sp.MR) for ia in range(sp.N_arms)] == [0, 3, 2]

    def test_lr_config_arm_ids_unchanged(self):
        path = find_config_file("20240714")
        sp = load_spectrograph_config(path, degrade=1.0)
        assert [spectro_arm(ia, sp.MR) for ia in range(sp.N_arms)] == [0, 1, 2]

    def test_throughput_grid_concatenated_with_istart_boundaries(self):
        path = find_config_file("20240714")
        sp = load_spectrograph_config(path, degrade=1.0)
        assert sp.istart[0] == 0
        assert sp.istart[-1] == sp.N_thr
        assert sp.istart.shape == (sp.N_arms + 1,)
        assert sp.l.shape == (sp.N_thr,)
        assert sp.T.shape == (sp.N_thr,)
        # Each arm's throughput-grid wavelengths are strictly increasing.
        for ia in range(sp.N_arms):
            lo, hi = sp.istart[ia], sp.istart[ia + 1]
            assert np.all(np.diff(sp.l[lo:hi]) > 0)

    def test_protected_fixture_loads(self):
        # tests/PFS.20211220.dat is a read-only acceptance fixture; must
        # never be modified, only read.
        sp = load_spectrograph_config(PROTECTED_FIXTURE, degrade=1.0)
        assert sp.N_arms == 3
        assert sp.MR is False
        np.testing.assert_allclose(sp.lmin, [380.0, 630.0, 940.0])

    def test_missing_optics_line_raises(self, tmp_path):
        bad = tmp_path / "bad.dat"
        bad.write_text("ARMS 1\nPARAM 0 380 650 4096 5\n")
        with pytest.raises(ValueError, match="OPTICS"):
            load_spectrograph_config(bad, degrade=1.0)


class TestFindConfigFile:
    def test_ave_lr(self):
        path = find_config_file("20240714", spectrograph="ave", mr_mode=False)
        assert path.name == "PFS.20240714.dat"
        assert path.is_file()

    def test_ave_mr(self):
        path = find_config_file("20240714", spectrograph="ave", mr_mode=True)
        assert path.name == "PFS.redMR.20240714.dat"
        assert path.is_file()

    def test_named_spectrograph_lr(self):
        path = find_config_file("20231001", spectrograph="sm1", mr_mode=False)
        assert path.name == "PFS.20231001.sm1.dat"
        assert path.is_file()

    def test_named_spectrograph_mr(self):
        path = find_config_file("20231001", spectrograph="sm1", mr_mode=True)
        assert path.name == "PFS.redMR.20231001.sm1.dat"
        assert path.is_file()

    def test_unknown_model_raises(self):
        with pytest.raises(FileNotFoundError):
            find_config_file("no-such-model")


class TestSpectroArm:
    @pytest.mark.parametrize("ia", [0, 1, 2])
    def test_lr_identity(self, ia):
        assert spectro_arm(ia, mr=False) == ia

    def test_mr_red_arm_remapped(self):
        assert spectro_arm(0, mr=True) == 0
        assert spectro_arm(1, mr=True) == 3
        assert spectro_arm(2, mr=True) == 2


class TestFieldInterp:
    VALUES = np.array([10.0, 20.0, 30.0, 40.0, 50.0])
    RFOV = 0.65

    def test_at_nodes(self):
        for i, expected in enumerate(self.VALUES):
            fieldang = i * self.RFOV / 4.0
            assert field_interp(self.VALUES, fieldang, self.RFOV) == pytest.approx(
                expected
            )

    def test_at_midpoints(self):
        # Halfway between node i and i+1: linear average.
        for i in range(4):
            fieldang = (i + 0.5) * self.RFOV / 4.0
            expected = 0.5 * (self.VALUES[i] + self.VALUES[i + 1])
            assert field_interp(self.VALUES, fieldang, self.RFOV) == pytest.approx(
                expected
            )

    def test_clamped_below_zero(self):
        assert field_interp(self.VALUES, -0.1, self.RFOV) == pytest.approx(
            self.VALUES[0]
        )

    def test_clamped_above_rfov(self):
        assert field_interp(self.VALUES, self.RFOV * 1.5, self.RFOV) == pytest.approx(
            self.VALUES[-1]
        )

    def test_at_rfov_exactly(self):
        # x = 4*fieldang/rfov = 4 -> i = 4 -> clamped to last node (gsetc.c
        # tests `i>=4`, so the top node itself is reached via the clamp
        # branch, not the interpolation branch).
        assert field_interp(self.VALUES, self.RFOV, self.RFOV) == pytest.approx(
            self.VALUES[-1]
        )
