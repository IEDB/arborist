#!/usr/bin/env python3

from pathlib import Path
from setuptools import setup

directory = Path(__file__).resolve().parent
with open(directory / 'README.md', encoding='utf-8') as f:
  long_description = f.read()

setup(
  name='protein_tree',
  version='0.1.0',
  description='Assigning IEDB source antigens and epitopes to their genes and proteins.',
  author='Daniel Marrama',
  author_email='dmarrama@lji.org',
  license='NPOSL-3.0',
  long_description=long_description,
  long_description_content_type='text/markdown',
  url='https://github.com/danielmarrama/protein_tree',
  packages = ['protein_tree'],
  classifiers=[
    "Programming Language :: Python :: 3",
    "License :: OSI Approved :: Non-Profit Open Software License 3.0 (Non-Profit OSL 3.0)",
  ],
  install_requires=[
    'biopython>=1.78',
    'lxml>=4.9.2',
    'mysql-connector-python>=8.0.32',
    'pandas>=1.3.0',
    'pepmatch>=0.9.4',
    'sqlalchemy>=2.0.4',
    'pytest>=7.3.2',
    'requests>=2.31.0',
    'bio-arc @ git+https://github.com/IEDB/ARC.git@master#egg=bio-arc'
  ],
  python_requires='>=3.8')