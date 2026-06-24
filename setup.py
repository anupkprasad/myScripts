#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Setup configuration for myScripts package.

This package contains custom scripts for bioinformatics, machine learning,
and molecular dynamics simulations.
"""

from setuptools import setup, find_packages
from pathlib import Path

# Read README
readme_file = Path(__file__).parent / "README.md"
long_description = readme_file.read_text(encoding="utf-8") if readme_file.exists() else ""

setup(
    name="myScripts",
    version="1.0.0",
    author="Anup K. Prasad",
    author_email="anupkprasad121@gmail.com",
    description="Custom scripts for bioinformatics, ML, and molecular dynamics simulations",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/yourusername/myScripts",
    packages=find_packages(exclude=["tests", "docs", "*.tests", "*.tests.*"]),
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Chemistry",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
    ],
    python_requires=">=3.8",
    install_requires=[
        "biopython>=1.80",
        "MDAnalysis>=2.0.0",
        "numpy>=1.21.0",
        "scipy>=1.7.0",
        "matplotlib>=3.5.0",
        "Pillow>=9.0.0",
        "PyYAML>=6.0",
        "dash>=2.0.0",
        "dash-bio>=1.0.0",
    ],
    extras_require={
        "ml": [
            "pyrosetta>=3.13",
            "torch>=1.10.0",
        ],
        "simulations": [
            "MDAnalysis>=2.0.0",
            "netCDF4>=1.5.0",
        ],
        "dev": [
            "pytest>=7.0.0",
            "pytest-cov>=4.0.0",
            "black>=22.0.0",
            "flake8>=4.0.0",
            "isort>=5.10.0",
            "mypy>=0.950",
            "sphinx>=4.5.0",
        ],
    },
    entry_points={
        "console_scripts": [
            # Add command-line scripts here if needed
            # "myscripts-cli=myScripts.cli:main",
        ],
    },
    include_package_data=True,
    package_data={
        "myScripts": [
            "config.example.yaml",
            "data/*.txt",
            "data/*.mapping",
        ],
    },
    zip_safe=False,
)
