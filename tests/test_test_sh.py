import os
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


def test_test_sh_detects_and_runs_pytest_with_same_python(tmp_path):
    fake_bin = tmp_path / "bin"
    fake_bin.mkdir()
    invocation_log = tmp_path / "python-invocations.txt"
    fake_python = fake_bin / "python"
    fake_python.write_text(
        "#!/bin/sh\n"
        'printf \'%s\\n\' "$*" >> "$INVOCATION_LOG"\n'
        'if [ "$1" = "-c" ]; then exit 0; fi\n'
        'if [ "$1" = "-m" ] && [ "$2" = "pytest" ]; then exit 0; fi\n'
        "exit 64\n",
        encoding="utf-8",
    )
    fake_python.chmod(0o755)
    env = {
        **os.environ,
        "INVOCATION_LOG": str(invocation_log),
        "MHCSEQS_PYTHON": str(fake_python),
        "TEST_SH_MAX": "1",
    }

    result = subprocess.run(
        ["bash", "test.sh", "--sentinel"],
        cwd=ROOT,
        env=env,
        capture_output=True,
        text=True,
        check=False,
    )

    assert result.returncode == 0, result.stderr
    assert invocation_log.read_text(encoding="utf-8").splitlines() == [
        "-c import xdist",
        "-m pytest -n 1 tests/ -v --sentinel",
    ]
