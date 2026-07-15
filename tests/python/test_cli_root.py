"""Tests for the `pfs-spec` umbrella `typer` app (`pfsspecsim.cli.app`).

These are smoke tests for the app-level wiring only -- that `etc` and `sim`
are both registered as subcommands, that the shared `--version` flag works
at the root, and that invoking the bare command shows help rather than
erroring. The subcommands' own option surfaces and behavior are covered by
`test_cli.py` (`etc`) and `test_sim_spec.py` (`sim`).
"""

from __future__ import annotations

from typer.testing import CliRunner

from pfsspecsim.cli import app

runner = CliRunner()


class TestUmbrellaHelpAndVersion:
    def test_help_lists_both_subcommands(self):
        result = runner.invoke(app, ["--help"])
        assert result.exit_code == 0
        assert "etc" in result.output
        assert "sim" in result.output

    def test_version_exits_zero_and_prints_pfs_spec(self):
        result = runner.invoke(app, ["--version"])
        assert result.exit_code == 0
        assert "pfs-spec" in result.output

    def test_bare_invocation_shows_help(self):
        # Click's `no_args_is_help` treats a missing subcommand as a usage
        # error (exit code 2) rather than a clean exit, but it still prints
        # the same help text as `--help` -- this asserts the help content
        # is shown, not that the exit code is 0.
        result = runner.invoke(app, [])
        assert result.exit_code == 2
        assert "etc" in result.output
        assert "sim" in result.output
