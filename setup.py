#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from setuptools import setup
import os

# from sys import platform
# from setuptools.command.install import install
from distutils.command.build import build
from distutils.command.clean import clean
from subprocess import call

# from multiprocessing import cpu_count


class GsetcBuild(build):
    def run(self):
        # build gsetc
        build_path = os.path.abspath(self.build_temp)
        print(build_path)

        target_files = ["src/gsetc.x", "src/gsetc_omp.x"]

        cmd_build = ["make", "-f", "Makefile", "all"]
        # cmd_install = ["make", "-f", "Makefile", "install"]

        def compile():
            call(cmd_build)

        # def install():
        # call(cmd_install)

        self.execute(compile, [], "Compiling gsetc and gsetc_omp")
        # self.execute(install, [], "Installing gsetc and gsetc_omp")

        # if not self.dry_run:
        #     for target in target_files:
        #         self.copy_file(target, os.path.join(build_path, "bin"))

        # run original build code
        build.run(self)


class GsetcClean(clean):
    def run(self):
        # run original build code
        clean.run(self)

        # build gsetc
        # build_path = os.path.abspath(self.build_temp)

        cmd = ["make", "-f", "Makefile", "clean"]

        def compile():
            call(cmd)

        self.execute(compile, [], "Cleaning gsetc")


def main():
    setup(
        name="pfsspecsim",
        version="0.1",
        description="PFS ETC and spectral simulator",
        author="Kiyoto Yabe",
        author_email="kiyoto.yabe@ipmu.jp",
        url="",
        install_requires=["scipy", "matplotlib"],
        zip_safe=False,
        include_package_data=True,
        license="",  # temporarily set to an empty string as installation failed.
        package_dir={"": "python"},
        packages=["pfsspecsim"],
        cmdclass={"build": GsetcBuild, "clean": GsetcClean},
    )


if __name__ == "__main__":
    main()
