from setuptools import setup, Command
import subprocess
import shutil
from setuptools.command.install import install
from pathlib import Path


class BuildFortran(Command):
    description = "Compile Fortran library with fpm and f2py"
    user_options = []

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        build_dir = Path("..") / "build" / "python"
        link_dir = build_dir.absolute() / "lib"
        incl_dir = build_dir.absolute() / "include"
        this_dir = Path(".")

        subprocess.check_call(
            ["fpm", "install", "--profile", "release", "--prefix", build_dir]
        )

        subprocess.check_call(
            [
                "f2py",
                "-m",
                "yaeos_compiled",
                f"-L{link_dir}",
                f"-I{incl_dir}",
                "-c",
                "yaeos/fortran_wrap/yaeos_c.f90",
                "-lyaeos",
                "--backend",
                "meson",
            ]
        )

        for file in this_dir.glob("yaeos_compiled.*"):
            shutil.move(file, this_dir / "yaeos/fortran_wrap")


class CustomInstall(install):
    def run(self):
        self.run_command("build_fortran")
        install.run(self)


setup(
    name="yaeos",
    version="0.3.5",
    packages=["yaeos"],
    cmdclass={
        "build_fortran": BuildFortran,
        "install": CustomInstall,
    },
    install_requires=["numpy", "fpm"],
)
