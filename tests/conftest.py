import os
import shutil
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
WORK = ROOT / "work"


@pytest.hookimpl(tryfirst=True)
def pytest_configure(config):
    # Keep tmp_path off /tmp (shared cluster); pass --basetemp to run sessions in parallel.
    if not config.option.basetemp:
        config.option.basetemp = str(WORK / "pytest")


def _missing(reason):
    if os.environ.get("RIBOTISH_ALLOW_SKIP") == "1":
        pytest.skip(reason)
    pytest.fail(reason + " (set RIBOTISH_ALLOW_SKIP=1 to skip instead)")


def _fork_bin():
    cand = Path(os.sys.executable).parent / "ribotish"
    return str(cand) if cand.exists() else shutil.which("ribotish")


@pytest.fixture(scope="session")
def fork_bin():
    b = _fork_bin()
    if not b:
        _missing("fork ribotish binary not found")
    return b


@pytest.fixture(scope="session")
def upstream_bin():
    b = WORK / "upstream" / ".venv" / "bin" / "ribotish"
    if not b.exists():
        _missing("upstream reference not installed (run scripts/install_upstream_reference.sh)")
    return str(b)


@pytest.fixture(scope="session")
def mini_dir():
    d = WORK / "mini"
    if not (d / "chr22.gtf").exists():
        _missing("mini fixture absent (run scripts/make_mini_fixture.sh)")
    return d


@pytest.fixture(scope="session")
def chr22_dir():
    d = WORK / "chr22"
    if not (d / "chr22.gtf").exists():
        _missing("work/chr22 fixture absent")
    return d


@pytest.fixture(scope="session")
def ref_root():
    return WORK / "ref"
