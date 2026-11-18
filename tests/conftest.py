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


DATA = ROOT / "tests" / "data" / "chr22_19M"
FIXTURES = ("u2os", "chen2020_ipsc")


@pytest.fixture(scope="session", params=FIXTURES)
def fixture_dir(request):
    """Each committed fixture: chr22:19,000,000-21,000,000 (tests/data/README.md)."""
    return DATA / request.param


@pytest.fixture(scope="session")
def mini_dir():
    """The fixture used by the edge-case tests."""
    return DATA / "u2os"


@pytest.fixture(scope="session")
def chr22_dir():
    d = WORK / "chr22"
    if not (d / "chr22.gtf").exists():
        pytest.skip("work/chr22 (whole chr22, local only) absent")
    return d


@pytest.fixture(scope="session")
def ref_root():
    return WORK / "ref"
