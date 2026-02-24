"""
Setuptools fallback for legacy tooling.

Primary package metadata lives in pyproject.toml. This file provides
equivalent core metadata for older installers that still execute setup.py.
"""

from pathlib import Path

from setuptools import find_packages, setup


ROOT = Path(__file__).resolve().parent
README = ROOT / "README.md"


def _read_readme() -> str:
    try:
        return README.read_text(encoding="utf-8")
    except FileNotFoundError:
        return "PyExoCross: A Python package for generating molecular and atomic spectra and cross sections."


if __name__ == "__main__":
    setup(
        name="pyexocross",
        version="1.0.0.dev1",
        description="A Python package for generating molecular and atomic spectra and cross sections",
        long_description=_read_readme(),
        long_description_content_type="text/markdown",
        author="Jingxin Zhang",
        author_email="jingxin.zhang.19@ucl.ac.uk",
        url="https://github.com/Beryl-Jingxin/PyExoCross",
        license="GPL-3.0-only",
        python_requires=">=3.8,<3.13",
        packages=find_packages(where="src", include=["pyexocross*"]),
        package_dir={"": "src"},
        include_package_data=True,
        package_data={"pyexocross": ["*.txt", "*.md", "database/meta/*.txt"]},
        install_requires=[
            "numpy>=1.20.0,<2.0.0",
            "pandas>=1.4.0,<3.0.0",
            "scipy>=1.7.0,<2.0.0",
            "matplotlib>=3.5.0,<4.0.0",
            "tqdm>=4.64.0,<5.0.0",
            "tabulate>=0.8.9,<1.0.0",
            "numexpr>=2.7.0,<3.0.0",
            "astropy>=5.0.0,<7.0.0",
            "dask>=2022.5.0,<2025.0.0",
            "pandarallel>=1.6.5,<2.0.0",
            "requests>=2.25.1,<3.0.0",
        ],
        extras_require={
            "dev": [
                "pytest",
                "pytest-cov",
                "black",
                "flake8",
            ]
        },
        entry_points={"console_scripts": ["pyexocross=pyexocross.cli:main"]},
        classifiers=[
            "Development Status :: 4 - Beta",
            "Intended Audience :: Science/Research",
            "Topic :: Scientific/Engineering :: Physics",
            "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
            "Programming Language :: Python :: 3",
            "Programming Language :: Python :: 3.8",
            "Programming Language :: Python :: 3.9",
            "Programming Language :: Python :: 3.10",
            "Programming Language :: Python :: 3.11",
            "Programming Language :: Python :: 3.12",
        ],
        keywords=["spectroscopy", "molecule", "atom", "cross section", "exomol", "hitran"],
    )
