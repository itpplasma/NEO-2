import importlib.util
import os
from pathlib import Path
import subprocess
import sys

import pytest


SCRIPT = Path(__file__).resolve().parents[2] / "PythonScripts" / "efit.py"


@pytest.fixture
def entry():
    spec = importlib.util.spec_from_file_location("neo2_efit_entry", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_conversion_entry_uses_checked_helper(entry, monkeypatch, tmp_path):
    from libneo import eqdsk_to_boozer

    source, output = tmp_path / "input.geqdsk", tmp_path / "output.bc"
    calls = []

    def convert(gfile, destination):
        calls.append((gfile, destination))
        destination.write_text("checked conversion\n")
        return {"metadata": "not part of the historical entry return"}

    monkeypatch.setattr(eqdsk_to_boozer, "convert_eqdsk_to_bc", convert)
    assert entry.efit_to_boozer(source, output) is None
    assert calls == [(source, output)]
    assert output.read_text() == "checked conversion\n"


def test_conversion_entry_propagates_rejection(entry, monkeypatch, tmp_path):
    from libneo import eqdsk_to_boozer

    output = tmp_path / "output.bc"
    output.write_text("existing output\n")

    def reject(gfile, destination):
        raise ValueError("equilibrium rejected")

    monkeypatch.setattr(eqdsk_to_boozer, "convert_eqdsk_to_bc", reject)
    with pytest.raises(ValueError, match="equilibrium rejected"):
        entry.efit_to_boozer(tmp_path / "invalid.geqdsk", output)
    assert output.read_text() == "existing output\n"


def test_conversion_cli_is_headless(tmp_path):
    package = tmp_path / "libneo"
    package.mkdir()
    (package / "__init__.py").write_text("")
    (package / "eqdsk_to_boozer.py").write_text(
        "from pathlib import Path\n"
        "def convert_eqdsk_to_bc(gfile, output):\n"
        "    Path(output).write_text(str(gfile) + '\\n')\n"
    )
    source, output = tmp_path / "input.geqdsk", tmp_path / "output.bc"
    env = dict(os.environ, PYTHONPATH=str(tmp_path), MPLBACKEND="Agg")
    env.pop("DISPLAY", None)
    result = subprocess.run(
        [sys.executable, str(SCRIPT), str(source), str(output)],
        env=env, capture_output=True, text=True, timeout=30,
    )
    assert result.returncode == 0, result.stderr
    assert output.read_text() == str(source) + "\n"
