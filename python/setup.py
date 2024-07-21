import shutil
import subprocess
from pathlib import Path

from setuptools import Command, setup
from setuptools.command.editable_wheel import editable_wheel
from setuptools.command.egg_info import egg_info


# =============================================================================
# Directories and constants
# =============================================================================
THIS_DIR = Path(__file__).parent
BUILD_DIR = (THIS_DIR.parent / "build" / "python").absolute()
LINK_DIR = BUILD_DIR / "lib"
INCL_DIR = BUILD_DIR / "include"
COMPILED_FLAG = THIS_DIR / "compiled_flag"

FFLAGS = "-g -fPIC -funroll-loops -fstack-arrays -Ofast -frepack-arrays -faggressive-function-elimination -fopenmp"  # noqa
CFLAGS = "-fPIC"


# =============================================================================
# Usefull functions
# =============================================================================
def pre_build():
    """Execute fpm and f2py compilations commands."""
    print(THIS_DIR)

    if COMPILED_FLAG.exists():
        return

    subprocess.check_call(
        [
            "fpm",
            "install",
            "--profile",
            "release",
            "--flag",
            f"{FFLAGS}",
            "--c-flag",
            f"{CFLAGS}",
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
            "-lyaeos", "-llapack",
            "--backend",
            "meson",
        ]
    )

    COMPILED_FLAG.touch()


def initial_compiled_clean():
    """Erase all compiled files from development directory"""
    # Clear fpm build
    if BUILD_DIR.exists():
        shutil.rmtree(BUILD_DIR)

    # Clear compiled files on compiled_files
    compiled_module_dir = THIS_DIR / "yaeos" / "compiled_module"

    if compiled_module_dir.exists():
        for so_file in compiled_module_dir.glob("*.so"):
            so_file.unlink()

    # Additionally, clear any .so files in the root directory if present
    for so_file in THIS_DIR.glob("yaeos_compiled*.so"):
        so_file.unlink()


def final_build_clean():
    """Clean the build of setuptools."""

    if (THIS_DIR / "build").exists():
        print((THIS_DIR / "build").absolute())
        shutil.rmtree(THIS_DIR / "build")

    # Clear compiled files on compiled_files
    compiled_module_dir = THIS_DIR / "yaeos" / "compiled_module"

    if compiled_module_dir.exists():
        for so_file in compiled_module_dir.glob("*.so"):
            so_file.unlink()

    # Clean COMPILED_FLAG
    if COMPILED_FLAG.exists():
        COMPILED_FLAG.unlink()


def move_compiled_to_editable_loc():
    """Move compiled files to 'compiled_module' directory"""

    for file in THIS_DIR.glob("yaeos_compiled.*"):
        target_dir = THIS_DIR / "yaeos" / "compiled_module"
        target_dir.mkdir(parents=True, exist_ok=True)

        shutil.move(file.absolute(), (target_dir / file.name).absolute())


def save_editable_compiled():
    """Temporaly save the editable compiled from install and sdist commands"""
    tmp_dir = THIS_DIR / "tmp_editable"

    if not tmp_dir.exists():
        tmp_dir.mkdir()

    compiled_module_dir = THIS_DIR / "yaeos" / "compiled_module"

    if compiled_module_dir.exists():
        for so_file in compiled_module_dir.glob("*.so"):
            if ((tmp_dir / so_file.name).absolute()).exists():
                ((tmp_dir / so_file.name).absolute()).unlink()

            shutil.move(
                so_file.absolute(), (tmp_dir / so_file.name).absolute()
            )


def restore_save_editable_compiled():
    tmp_dir = THIS_DIR / "tmp_editable"

    compiled_module_dir = THIS_DIR / "yaeos" / "compiled_module"

    if tmp_dir.exists():
        for so_file in tmp_dir.glob("*.so"):
            shutil.move(
                so_file.absolute(),
                (compiled_module_dir / so_file.name).absolute(),
            )

        tmp_dir.rmdir()


# =============================================================================
# Build command
# =============================================================================
class BuildFortran(Command):
    description = "Compile Fortran library with fpm and f2py"
    user_options = []

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        dir = str(THIS_DIR.absolute())

        if ("build" in dir) or ("check-manifest" in dir):
            # Do not compile, we are building, the compilation has been already
            # done at this point.
            ...
        else:
            pre_build()


# =============================================================================
# - Building for developers (editable installation)
#      pip install -e .
# =============================================================================
class CustomEditable(editable_wheel):
    def run(self):
        self.run_command("build_fortran")
        move_compiled_to_editable_loc()
        save_editable_compiled()

        # Run base editable_wheel run method
        super().run()


# =============================================================================
# - Custom egg_info command
# =============================================================================
class CustomEgg(egg_info):
    def run(self):
        self.run_command("build_fortran")
        move_compiled_to_editable_loc()
        super().run()


# =============================================================================
# Call setup
# =============================================================================
name = "yaeos"
version = "0.3.0"
author = "Federico E. Benelli"
author_email = "federico.benelli@mi.unc.edu.ar"
maintainer = "Federico E. Benelli"
maintainer_email = "federico.benelli@mi.unc.edu.ar"
lic = "MPL"


save_editable_compiled()

initial_compiled_clean()

setup(
    name=name,
    version=version,
    author=author,
    author_email=author_email,
    maintainer=maintainer,
    maintainer_email=maintainer_email,
    description="",
    license=lic,
    keywords="thermodynamics equation-of-state",
    url="https://github.com/ipqa-research/yaeos",
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research/Engineering",
        "Topic :: Thermodynamics",
        "License :: OSI Approved :: Mozilla Public License 2.0 (MPL 2.0)",
    ],
    install_requires=["numpy"],
    cmdclass={
        "build_fortran": BuildFortran,
        "editable_wheel": CustomEditable,
        "egg_info": CustomEgg
    },
    packages=["yaeos"],
    package_data={"yaeos": ["compiled_module/*.so"]},
    include_package_data=True,
)

final_build_clean()

restore_save_editable_compiled()
