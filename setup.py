#!/usr/bin/env python
"""
# Author: RyanYip
# Created Time : Wen 21 Fri 2021 10:42:37 AM CST
# File Name: setup.py
# Description:
"""
from setuptools import setup, find_packages

with open('requirements.txt') as f:
    requirements = f.read().splitlines()

setup(name='scpcat',
      version='1.0.0',
      packages=find_packages(),
      description='Single-Cell RNA-seq Simple Pipeline Analysis',
      long_description='',


      author='RyanYip',
      author_email='ryanyip_@hotmail.com',
      url='https://github.com/RyanYip-Kat/scpCAT',
      scripts=['scpCAT.py'],
      install_requires=requirements,
      python_requires='>3.7.7',
      license='MIT',

      classifiers=[
          'Development Status :: 4 - Beta',
          'Intended Audience :: Science/Research',
          'License :: OSI Approved :: MIT License',
          'Programming Language :: Python :: 3.7',
          'Operating System :: MacOS :: MacOS X',
          'Operating System :: Microsoft :: Windows',
          'Operating System :: POSIX :: Linux',
          'Topic :: Scientific/Engineering :: Bio-Informatics',
     ],
     )
