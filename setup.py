#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup


def main():
    setup(name='pfsspecsim',
          version='0.1',
          description='PFS ETC and spectral simulator',
          author='Kiyoto Yabe',
          author_email='kiyoto.yabe@ipmu.jp',
          url='',
          install_requires=['scipy', 'matplotlib'],
          zip_safe=False,
          include_package_data=True,
          license="",  # temporarily set to an empty string as installation failed.
          package_dir={'': 'python'},
          packages=['pfsspecsim'],
          )

if __name__ == '__main__':
    main()
