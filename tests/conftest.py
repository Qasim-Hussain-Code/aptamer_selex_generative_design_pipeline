import sys
from pathlib import Path
import shutil
import tempfile

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))
_TEMP_RUN=None


def pytest_configure(config):
    global _TEMP_RUN
    if config.option.basetemp is None:
        # Global Windows temp directories can belong to another execution identity.
        # Keep each test run inside the measured project and avoid sharing private caches.
        cache=(ROOT/".cache").resolve()
        cache.mkdir(exist_ok=True)
        _TEMP_RUN=Path(tempfile.mkdtemp(prefix="pytest_run_",dir=cache)).resolve()
        config.option.basetemp=str(_TEMP_RUN)


def pytest_sessionfinish(session,exitstatus):
    if _TEMP_RUN is not None and _TEMP_RUN.resolve().is_relative_to((ROOT/".cache").resolve()):
        shutil.rmtree(_TEMP_RUN,ignore_errors=False)
