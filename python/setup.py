from setuptools import setup, Command
import subprocess
import shutil
from setuptools.command.install import install
from setuptools.command.sdist import sdist
from setuptools.command.editable_wheel import editable_wheel
from pathlib import Path
import sysconfig


# =============================================================================
# Directories
# =============================================================================
THIS_DIR = Path(__file__).parent
BUILD_DIR = (THIS_DIR.parent / "build" / "python").absolute()
LINK_DIR = BUILD_DIR / "lib"
INCL_DIR = BUILD_DIR / "include"


# =============================================================================
# Usefull functions
# =============================================================================
def pre_build():
    subprocess.check_call(
        [
            "fpm",
            "install",
            "--profile",
            "release",
            "--flag",
            "-g -fPIC -funroll-loops -fstack-arrays -Ofast -frepack-arrays -faggressive-function-elimination -fopenmp",
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
    # Erase all compiled files just in case before building again
    compiled_module_dir = THIS_DIR / "yaeos" / "compiled_module"

    if compiled_module_dir.exists():
        for so_file in compiled_module_dir.glob("*.so"):
            so_file.unlink()


def move_compiled_to_editable_loc():
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
        pre_build()


# =============================================================================
# - Normal build and installation:
#      pip install .
# =============================================================================
class CustomInstall(install):
    def run(self):
        clean_editable_compiled()
        
        try:
            self.run_command("build_fortran")
        except:
            ...

        site_packages_dir = Path(sysconfig.get_path("purelib"))

        for file in THIS_DIR.glob("yaeos_compiled.*"):
            target_dir = site_packages_dir / "yaeos" / "compiled_module"
            target_dir.mkdir(parents=True, exist_ok=True)
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
# =============================================================================
class CustomBuild(sdist):
    def run(self):       
        clean_editable_compiled()

        self.run_command("build_fortran")

        move_compiled_to_editable_loc()

        # Run base sdist run method
        super().run()


# =============================================================================
# Call setup
# =============================================================================
setup(
    name="yaeos",
    version="0.3.5",
    package_data={
        'yaeos': ['compiled_module/*.so']
    },
    include_package_data=True,
    cmdclass={
        "build_fortran": BuildFortran,
        "editable_wheel": CustomEditable,
        "install": CustomInstall,
        "sdist": CustomBuild,
    },
    install_requires=["numpy"],
)
