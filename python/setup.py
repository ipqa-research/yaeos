from setuptools import setup, Command
import subprocess
import os
import glob
import shutil
from setuptools.command.install import install

class BuildFortran(Command):
    description = "Compile Fortran library with fpm and f2py"
    user_options = []

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):      
        # Compile with fpm
        os.chdir("./../")

        subprocess.check_call([
           "fpm", "install", "--profile", "release", "--prefix", "build/python"
        ])
        
        # Compile with f2py
        os.chdir("./python/")

        subprocess.check_call([
           "f2py", "-m", "yaeos_compiled",
           "-L../build/python/lib/",
           "-I../build/python/include/",
           "-c", "yaeos/fortran_wrap/yaeos_c.f90", "-lyaeos"
        ])

        for file in glob.glob("yaeos_compiled.*"):
            shutil.move(file, "yaeos/fortran_wrap")

class CustomInstall(install):
    def run(self):
        self.run_command("build_fortran")
        install.run(self)

setup(
    name="yaeos",
    version="0.1",
    packages=["yaeos"],
    cmdclass={
        "build_fortran": BuildFortran,
        "install": CustomInstall,
    },
    install_requires=[
       "numpy",
       "fpm"
    ],
)
