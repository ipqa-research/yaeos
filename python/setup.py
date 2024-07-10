import shutil
import subprocess
import sysconfig
from pathlib import Path

from setuptools import Command, setup
from setuptools.command.editable_wheel import editable_wheel
from setuptools.command.install import install
from setuptools.command.sdist import sdist


# =============================================================================
# Directories
# =============================================================================
THIS_DIR = Path(__file__).parent
BUILD_DIR = (THIS_DIR.parent / "build" / "python").absolute()
LINK_DIR = BUILD_DIR / "lib"
INCL_DIR = BUILD_DIR / "include"

# Signal to skip compilation when building
compilation_skip_signal = not (THIS_DIR / "tox.ini").exists()

print()


# =============================================================================
# Usefull functions
# =============================================================================
def pre_build():
    """Execute fpm and f2py compilations commands."""

    subprocess.check_call(
        [
            "fpm",
            "install",
            "--profile",
            "release",
            "--flag",
            "-g -fPIC -funroll-loops -fstack-arrays -Ofast -frepack-arrays -faggressive-function-elimination -fopenmp", # noqa
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


def clean_editable_compiled():
    """Erase all compiled files from development directory"""
    # Clear fpm build
    if BUILD_DIR.exists():
        shutil.rmtree(BUILD_DIR)

    # Clear compiled files on compiled_files
    compiled_module_dir = THIS_DIR / "yaeos" / "compiled_module"

    if compiled_module_dir.exists():
        for so_file in compiled_module_dir.glob("*.so"):
            so_file.unlink()


def move_compiled_to_editable_loc():
    """Move compiled files to 'compiled_module' directory"""

    for file in THIS_DIR.glob("yaeos_compiled.*"):
        target_dir = THIS_DIR / "yaeos" / "compiled_module"
        target_dir.mkdir(parents=True, exist_ok=True)

        shutil.move(file.absolute(), (target_dir / file.name).absolute())


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
        if compilation_skip_signal:
            # Do not compile, we are building, the compilation has been already
            # done at this point.
            ...
        else:
            pre_build()


# =============================================================================
# - Normal build and installation:
#      pip install .
# =============================================================================
class CustomInstall(install):
    def run(self):
        clean_editable_compiled()

        self.run_command("build_fortran")

        site_packages_dir = Path(sysconfig.get_path("purelib"))

        for file in THIS_DIR.glob("yaeos_compiled.*"):
            target_dir = site_packages_dir / "yaeos" / "compiled_module"
            target_dir.mkdir(parents=True, exist_ok=True)

            if (target_dir / file.name).exists():
                (target_dir / file.name).unlink()

            shutil.move(file, target_dir)

        super().run()


# =============================================================================
# - Building for developers (editable installation)
#      pip install -e .
# =============================================================================
class CustomEditable(editable_wheel):
    def run(self):
        clean_editable_compiled()

        self.run_command("build_fortran")

        move_compiled_to_editable_loc()

        # Run base editable_wheel run method
        super().run()


# =============================================================================
# - Python Build for distribution
#      pip install build
#      python3 -m build
# =============================================================================
class CustomBuild(sdist):
    def run(self):
        # Clean compiled files and recompile as an editable installation
        clean_editable_compiled()

        self.run_command("build_fortran")

        move_compiled_to_editable_loc()

        # Run base sdist run method
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
        "install": CustomInstall,
        "sdist": CustomBuild,
    },
    package_data={"yaeos": ["compiled_module/*.so"]},
    include_package_data=True,
)
