import subprocess, pathlib, importlib

def _bin(name: str) -> pathlib.Path:
    """
    Return the path to a compiled executable in autoseed/bin.
    """
    pkg = importlib.import_module("autoseed")
    basedir = pathlib.Path(pkg.__file__).resolve().parent
    exe = basedir / "bin" / name
    return exe

def localmax(args):
    return subprocess.run(
        [str(_bin("localmax")), *map(str, args)],
        check=True,
        capture_output=True,
        text=True,
    )

def seedext(args):
    return subprocess.run(
        [str(_bin("seedext")), *map(str, args)],
        check=True,
        capture_output=True,
        text=True,
    )

def sim(args):
    return subprocess.run(
        [str(_bin("sim")), *map(str, args)],
        check=True,
        capture_output=True,
        text=True,
    )

def genint(args):
    return subprocess.run(
        [str(_bin("genint")), *map(str, args)],
        check=True,
        capture_output=True,
        text=True,
    )
