import os
import pathlib
import subprocess
import tempfile
import tarfile
import urllib.request
import shutil

from setuptools import setup, Command
from setuptools.dist import Distribution
from setuptools.command.build_py import build_py as _build_py


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
                "https://github.com/jutaipal/localmax-motif/archive/6a0360774efdcb584a34d37dc939dbcd23e2178c.tar.gz",
                "localmax-motif.c",
            ),
            (
                "seedextender",
                "https://github.com/jutaipal/seedextender/archive/07eb219f800c5ab872281e097b15295c62eef6be.tar.gz",
                "seedextender.c",
            ),
            (
                "motifsimilarity",
                "https://github.com/jutaipal/motifsimilarity/archive/33fd7384f77ea785f595a1276904e4fc12cf5e61.tar.gz",
                "motifsimilarity.c",
            ),
            (
                "genint-PWM",
                "https://github.com/jutaipal/genint/archive/226eb7fb7f4719efd78e6207fcc386a0305bbfad.tar.gz",
                "genint.c",  # v0.9: has -keepfeatures, -threshold, -nofilter, -protect_flanks
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
                cmd = [cc, "-O3", str(src_path), "-o", str(out_path), "-lm", "-pthread"]
                subprocess.check_call(cmd)

        print("C helper binaries built into", bin_dir)


class BinaryDistribution(Distribution):
    """Mark wheel as platform-specific, since it contains compiled binaries."""

    def has_ext_modules(self):
        return True


class build_py(_build_py):
    """Run BuildC before the normal build_py."""

    def run(self):
        self.run_command("build_c")
        super().run()


setup(
    # All metadata / options come from pyproject.toml; we only hook commands here.
    distclass=BinaryDistribution,
    cmdclass={
        "build_py": build_py,
        "build_c": BuildC,
    },
)
