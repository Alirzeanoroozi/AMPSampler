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
        # Conda-friendly modern TF stack (Python 3.9-3.11).
        'tensorflow>=2.13,<2.14',
        'tensorflow-probability>=0.21,<0.22',
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
        'argparse',
        'tqdm>=4.66',
        'biopython==1.83',
        'gdown>=5.2.0'
    ],
    setup_requires=['wheel']
)
