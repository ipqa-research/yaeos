from setuptools import setup, Command
import subprocess
import shutil
from setuptools.command.install import install
from pathlib import Path
import sysconfig


# =============================================================================
# Building
# =============================================================================
THIS_DIR = Path(".")
BUILD_DIR = (Path("..") / "build" / "python").absolute()
LINK_DIR = BUILD_DIR / "lib"
INCL_DIR = BUILD_DIR / "include"


def pre_build():

    subprocess.check_call(
        [
            "fpm",
            "install",
            "--profile",
            "release",
            "--flag",
            "-g -fPIC -funroll-loops -fstack-arrays -Ofast -frepack-arrays  -faggressive-function-elimination -fopenmp",
            "--c-flag",
            "-fPIC",
            "--prefix",
            BUILD_DIR,
        ]
    )

    subprocess.check_call(
        [
            "f2py",
            "-m",
            "yaeos_compiled",
            f"-L{LINK_DIR}",
            f"-I{INCL_DIR}",
            "-c",
            "yaeos/fortran_wrap/yaeos_c.f90",
            "-lyaeos",
            "--backend",
            "meson",
        ]
    )


class BuildFortran(Command):
    description = "Compile Fortran library with fpm and f2py"
    user_options = []

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        pre_build()
        site_packages_dir = Path(sysconfig.get_path("purelib"))

        for file in THIS_DIR.glob("yaeos_compiled.*"):
            target_dir = site_packages_dir / "yaeos" / "compiled_module"
            target_dir.mkdir(parents=True, exist_ok=True)
            shutil.move(file, target_dir)


# =============================================================================
# Building for developers (editable installation)
# install the package on the enviroment (pip install -e .)
# Build fortran to the editable command with:
#
# python3 setup.py build_fortran_editable
# =============================================================================
class BuildFortranEditable(BuildFortran):

    def run(self):
        pre_build()
        for file in THIS_DIR.glob("yaeos_compiled.*"):
            target_dir = Path(".") / "yaeos" / "compiled_module"
            target_dir.mkdir(parents=True, exist_ok=True)

            shutil.move(file.absolute(), target_dir.absolute()/file)


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
        "build_fortran_editable": BuildFortranEditable,
        "install": CustomInstall,
    },
    install_requires=["numpy", "fpm"],
    package_data={
        "yaeos": [
            "compiled_module/",
        ],
    },
    include_package_data=True,
    exclude_package_data={
        "yaeos": ["__pycache__", "*.f90", "*.egg-info"],
    },
)
