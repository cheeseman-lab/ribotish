import subprocess
import sys
from pathlib import Path


def test_version():
    exe = Path(sys.executable).parent / "ribotish"
    out = subprocess.run([str(exe), "--version"], capture_output=True, text=True, check=True)
    assert "0.2.8+cheeseman.1" in out.stdout + out.stderr
