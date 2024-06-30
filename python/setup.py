from setuptools import setup, Command
import subprocess
import shutil
from setuptools.command.install import install
from pathlib import Path
import sysconfig


# =============================================================================
# Building
# =============================================================================
def build_library():
    build_dir = Path("..") / "build" / "python"
    link_dir = build_dir.absolute() / "lib"
    incl_dir = build_dir.absolute() / "include"
    this_dir = Path(".")

    subprocess.check_call(
        [
            "fpm",
            "install",
            "--profile",
            "release",
            "--flags",
            "-fPIC",
            "--c-flags",
            "-fPIC",
            "--prefix",
            build_dir,
        ]
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

    site_packages_dir = Path(sysconfig.get_path("purelib"))

    for file in this_dir.glob("yaeos_compiled.*"):
        target_dir = site_packages_dir / "yaeos" / "compiled_module"
        target_dir.mkdir(parents=True, exist_ok=True)
        shutil.move(file, target_dir)


class BuildFortran(Command):
    description = "Compile Fortran library with fpm and f2py"
    user_options = []

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        build_library()


# =============================================================================
# Building for developers (editable installation)
# install the package on the enviroment (pip install -e .)
# Build fortran to the editable command with:
#
# python3 setup.py build_fortran_editable
# =============================================================================
class BuildFortranEditable(Command):
    description = "Compile Fortran library with fpm and f2py"
    user_options = []

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        build_library()


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
