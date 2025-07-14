from setuptools import setup, find_packages

setup(
    name="halogenase_miner",
    version="0.1.20",
    description="A Python package to identify categorize and mine for halogenases",
    url="https://github.com/JudSze/hal_miner",
    author="Judit Szenei",
    author_email="szenei@dtu.dk",
    license="Affero",
    packages=find_packages(),
    install_requires=["pyhmmer", "pandas"],
    package_data={
    'halogenase_miner': [
        'motif_db/specific_enzymes.json',
        'phmm_db/*.hmm',
    ],
},
)