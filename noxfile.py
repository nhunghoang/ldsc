import os
import pathlib
import nox

os.environ.update({"PDM_IGNORE_SAVED_PYTHON": "1"})

VENV_DIR = pathlib.Path('./.venv').resolve()

@nox.session(tags=["tests"], python=["3.8", "3.9", "3.10", "3.11", "3.12"])
def run_test(session):
    session.run_always("pdm", "install", "--without", "dev", external=True)

    python = os.fsdecode(VENV_DIR.joinpath("bin/python"))
    
    session.run(python, "-m", "unittest", "discover", "-s", "./tests/", "-p", "*.py")