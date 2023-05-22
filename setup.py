#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import setuptools
from glob import glob
setuptools.setup(name="miBrain",
                 version="0.0.1",
                 author="Benjamin James",
                 author_email="benjames@mit.edu",
                 url="https://github.com/KellisLab/miBrain",
                 license="GPL",
                 install_requires=[
                     "numpy",
                     "scipy",
                     "scikit-learn",
                     "pandas",
                     "cellphonedb",
                     "umap-learn",
                     "matplotlib",
                     "seaborn",
                     "tqdm"
                 ],
                 packages=setuptools.find_packages("."),
                 test_suite="test",
                 scripts=glob("scripts/*.py") + glob("scripts/*.sh")
                 )
