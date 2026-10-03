import argparse
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
from common import ROOT, atomic_text


def main():
    p=argparse.ArgumentParser(description="Verify committed code without data caches in an isolated local clone")
    p.parse_args()
    cache=(ROOT/".cache").resolve()
    cache.mkdir(exist_ok=True)
    path=Path(tempfile.mkdtemp(prefix="clean_clone_",dir=cache)).resolve()
    if not path.is_relative_to(cache):raise ValueError("Unsafe clone path")
    try:
        subprocess.run(["git","-c",f"safe.directory={ROOT.as_posix()}","clone","--no-hardlinks","--local",str(ROOT),str(path)],check=True)
        env=os.environ.copy()
        env["PYTHON"]=sys.executable
        env.pop("PYTHONPATH",None)
        # Pass the actual running Bash path; Windows' system bash can invoke a different WSL context.
        bash=env.get("PIPELINE_BASH_EXECUTABLE","bash")
        smoke=subprocess.run([bash,"run_all.sh","--mode","smoke"],cwd=path,env=env,text=True,capture_output=True)
        tests=subprocess.run([sys.executable,"-m","pytest","-q"],cwd=path,env=env,text=True,capture_output=True)
        atomic_text(ROOT/"logs/clean_clone_verification.txt",smoke.stdout+smoke.stderr+tests.stdout+tests.stderr)
        if smoke.returncode or tests.returncode:raise RuntimeError("Clean-clone smoke/tests failed; see logs/clean_clone_verification.txt")
        print("Clean clone: synthetic workflow and full Python suite passed without research data")
    finally:
        # Verify the resolved deletion target remains in the named project cache.
        if path.resolve().is_relative_to(cache):shutil.rmtree(path)


if __name__=="__main__":
    main()
