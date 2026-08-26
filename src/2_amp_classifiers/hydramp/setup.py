from setuptools import find_packages
from setuptools import setup

setup(
    name='HydrAMP',
    version='1.2.0',
    description='Python package for peptide generation',
    author='Paulina Szymczak',
    author_email='szymczak.pau@gmail.com',
    url='https://hydramp.mimuw.edu/',
    packages=find_packages(),
    python_requires='>=3.9,<3.12',
    install_requires=[
        'tensorflow>=2.11,<2.12',
        'tensorflow-probability>=0.19,<0.20',
        'cloudpickle>=2.2',
        'numpy>=1.23,<2.0',
        'pandas>=1.5',
        'scikit-learn>=1.2',
        'modlamp>=4.3',
        'matplotlib>=3.7',
        'protobuf>=3.20,<5',
        'seaborn>=0.12',
        'setuptools>=65',
        'joblib>=1.2',
        'tqdm>=4.66',
        'biopython==1.83',
        'gdown>=5.2.0',
    ],
    setup_requires=['wheel'],
)
