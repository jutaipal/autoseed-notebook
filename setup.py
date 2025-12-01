import os, pathlib, subprocess, tarfile, tempfile, urllib.request
from setuptools import setup
from setuptools.command.build_py import build_py

# (binary_name, repo_url, c_filename)
REPOS = [
    {
        "name": "localmax-motif",
        "url": "https://github.com/jutaipal/localmax-motif/archive/refs/heads/main.tar.gz",
        "c_file": "localmax-motif.c",
        "output": "localmax-motif",
    },
    {
        "name": "seedextender",
        "url": "https://github.com/jutaipal/seedextender/archive/refs/heads/main.tar.gz",
        "c_file": "seedextender.c",
        "output": "seedextender",
    },
    {
        "name": "motifsimilarity",
        "url": "https://github.com/jutaipal/motifsimilarity/archive/refs/heads/main.tar.gz",
        "c_file": "motifsimilarity.c",
        "output": "motifsimilarity",
    },
    {
        # genint-PWM lives in the genint repo, branch "newmain"
        "name": "genint-PWM",
        "url": "https://github.com/jutaipal/genint/archive/refs/heads/newmain.tar.gz",
        "c_file": "genint-PWM.c",
        "output": "genint-PWM",
    },
]



class BuildC(build_py):
    def run(self):
        here = pathlib.Path(__file__).resolve().parent
        bin_dir = here / "src" / "autoseed" / "bin"
        bin_dir.mkdir(parents=True, exist_ok=True)

        with tempfile.TemporaryDirectory() as tmp:
            tmp = pathlib.Path(tmp)

            for bin_name, url, c_name in REPOS:
                print(f"Downloading: {url}")
                tar_path = tmp / f"{bin_name}.tar.gz"
                urllib.request.urlretrieve(url, tar_path)

                print(f"Extracting: {tar_path}")
                with tarfile.open(tar_path, "r:gz") as tf:
                    tf.extractall(tmp)

                # find extracted dir (first subdirectory)
                extracted_dirs = [p for p in tmp.iterdir() if p.is_dir()]
                repo_root = extracted_dirs[0]

                cfile = None
                for root, dirs, files in os.walk(repo_root):
                    if c_name in files:
                        cfile = pathlib.Path(root) / c_name
                        break
                if cfile is None:
                    raise RuntimeError(f"Could not find {c_name} in {repo_root}")

                outfile = bin_dir / bin_name
                print(f"Compiling {cfile} -> {outfile}")

                subprocess.check_call(["cc", "-O3", str(cfile), "-o", str(outfile)])

        super().run()


setup(
    cmdclass={"build_py": BuildC},
)
