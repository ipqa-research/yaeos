from setuptools import setup, Command, find_packages
import subprocess
import shutil
from setuptools.command.install import install
from setuptools.command.develop import develop
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


# =============================================================================
# Building for developers (editable installation)
# install the package on the enviroment (pip install -e .)
# Build fortran to the editable command with:
#
# python3 setup.py build_fortran_editable
# =============================================================================
class CustomDevelop(develop):
    def run(self):
        print("this command is editable")
        
        self.run_command("build_fortran")

        for file in THIS_DIR.glob("yaeos_compiled.*"):
            target_dir = Path(".") / "yaeos" / "compiled_module"
            target_dir.mkdir(parents=True, exist_ok=True)

            shutil.move(file.absolute(), target_dir.absolute()/file)

        super().run()


class CustomInstall(install):
    def run(self):
        print("this command is no editable")
        self.run_command("build_fortran")

        site_packages_dir = Path(sysconfig.get_path("purelib"))

        for file in THIS_DIR.glob("yaeos_compiled.*"):
            target_dir = site_packages_dir / "yaeos" / "compiled_module"
            target_dir.mkdir(parents=True, exist_ok=True)
            shutil.move(file, target_dir)

        super().run()


setup(
    name="yaeos",
    version="0.3.5",
    packages= find_packages(),
    cmdclass={
        "build_fortran": BuildFortran,
        "develop": CustomDevelop,
        "install": CustomInstall,
    },
)
