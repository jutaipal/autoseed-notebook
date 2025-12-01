import os
import pathlib
import subprocess
import tempfile
import tarfile
import urllib.request
import shutil

from setuptools import setup
from setuptools.command.build_py import build_py as _build_py
from distutils.core import Command


class BuildC(Command):
    """Download and build C helper binaries from GitHub."""

    description = "Download and build C helper binaries from GitHub"
    user_options = []

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        # Where to put the compiled binaries inside the package
        pkg_root = pathlib.Path(__file__).parent
        bin_dir = pkg_root / "src" / "autoseed" / "bin"
        bin_dir.mkdir(parents=True, exist_ok=True)

        # (executable name, tarball URL, C source file name)
        REPOS = [
            (
                "localmax-motif",
                "https://github.com/jutaipal/localmax-motif/archive/refs/heads/main.tar.gz",
                "localmax-motif.c",
            ),
            (
                "seedextender",
                "https://github.com/jutaipal/seedextender/archive/refs/heads/main.tar.gz",
                "seedextender.c",
            ),
            (
                "motifsimilarity",
                "https://github.com/jutaipal/motifsimilarity/archive/refs/heads/main.tar.gz",
                "motifsimilarity.c",
            ),
            (
                "genint-PWM",
                "https://github.com/jutaipal/genint/archive/refs/heads/newmain.tar.gz",
                "genint-PWM.c",
            ),
        ]

        cc = os.environ.get("CC", "cc")

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = pathlib.Path(tmpdir)

            for exe_name, tar_url, c_name in REPOS:
                print("Downloading:", tar_url, "for", exe_name)
                tar_path = tmpdir / (exe_name + ".tar.gz")
                urllib.request.urlretrieve(tar_url, tar_path)

                # Extract only the wanted .c file
                with tarfile.open(tar_path, "r:gz") as tf:
                    members = [
                        m for m in tf.getmembers()
                        if m.name.endswith("/" + c_name)
                    ]
                    if not members:
                        raise RuntimeError(f"Could not find {c_name} in {tar_url}")
                    member = members[0]
                    src_path = tmpdir / c_name
                    with tf.extractfile(member) as src_f, open(src_path, "wb") as out_f:
                        shutil.copyfileobj(src_f, out_f)

                out_path = bin_dir / exe_name
                print(f"Compiling {c_name} -> {out_path}")
                cmd = [cc, "-O3", "-std=c11", str(src_path), "-o", str(out_path)]
                subprocess.check_call(cmd)

        print("C helper binaries built into", bin_dir)


class build_py(_build_py):
    """Run BuildC before the normal build_py."""

    def run(self):
        self.run_command("build_c")
        super().run()


setup(
    # All metadata / options come from pyproject.toml; we only hook commands here.
    cmdclass={
        "build_py": build_py,
        "build_c": BuildC,
    },
)
