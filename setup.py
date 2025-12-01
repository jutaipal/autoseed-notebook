import os, pathlib, subprocess, tarfile, tempfile, urllib.request
from setuptools import setup
from setuptools.command.build_py import build_py

# GitHub repos pinned to HEAD or specify commit hashes if desired
REPOS = [
    ("localmax", "https://github.com/jutaipal/kmercount_with_localmax/archive/refs/heads/main.tar.gz"),
    ("seedext",  "https://github.com/jutaipal/seedextender/archive/refs/heads/main.tar.gz"),
    ("sim",      "https://github.com/jutaipal/motifsimilarity/archive/refs/heads/main.tar.gz"),
    ("genint",   "https://github.com/jutaipal/multinomial_motif_generator/archive/refs/heads/main.tar.gz")
]

def find_c_file(path: pathlib.Path):
    """Return the first .c file inside a repo."""
    for root, dirs, files in os.walk(path):
        for f in files:
            if f.endswith(".c"):
                return pathlib.Path(root) / f
    raise RuntimeError(f"No C file found inside {path}")

class BuildC(build_py):
    def run(self):
        here = pathlib.Path(__file__).resolve().parent
        bin_dir = here / "src" / "autoseed" / "bin"
        bin_dir.mkdir(parents=True, exist_ok=True)

        with tempfile.TemporaryDirectory() as tmp:
            tmp = pathlib.Path(tmp)

            for name, url in REPOS:
                print(f"Downloading: {url}")
                tar_path = tmp / f"{name}.tar.gz"
                urllib.request.urlretrieve(url, tar_path)

                print(f"Extracting: {tar_path}")
                with tarfile.open(tar_path, "r:gz") as tf:
                    tf.extractall(tmp)

                # find extracted dir (name may vary)
                extracted_dirs = [p for p in tmp.iterdir() if p.is_dir() and name in p.name]
                if not extracted_dirs:
                    extracted_dirs = [p for p in tmp.iterdir() if p.is_dir()]
                repo_root = extracted_dirs[0]

                cfile = find_c_file(repo_root)
                outfile = bin_dir / name
                print(f"Compiling {cfile} -> {outfile}")

                subprocess.check_call(["gcc", "-O3", str(cfile), "-o", str(outfile)])

        super().run()


setup(
    cmdclass={"build_py": BuildC},
)
