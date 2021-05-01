#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 29 07:24:23 2021

@author: maryamzaheri
"""

from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    use_scm_version=True,
    setup_requires=['setuptools_scm'],
    install_requires=['setuptools_scm'],
    tests_require=['pytest', 'flake8'],
    name='covid_spike_lineage',
    author='Maryam Zaheri',
    author_email='lastname@gmail.com',
    description='Identify the lineage of SARS-CoV-2 based on a part of Spike gene.',
    packages=find_packages('src'),
    package_dir={'': 'src'},
    entry_points={
        'console_scripts': ['covid_spike_lineage = covid_spike_lineage.cli:main']
    },
    license='MIT',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/medvir/covid_spike_lineage",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)