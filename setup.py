"""Pip installation script.

TODO:
    - On installation, ask where to put option file directory (with default)
"""
from setuptools import find_packages, setup

setup(
    name='matsim',
    version="0.1",
    description=(
        "Materials simulation management code."
    ),
    author='Adam J. Plowman, Maria S. Yankova',
    packages=find_packages(),
    install_requires=[
        'plotly',
        'spglib',
        'numpy',
        'matplotlib',
        'mendeleev',
        'dropbox',
        'PyYAML',
        'beautifultable',
        'passlib',
        'pymysql',
    ]
)
