from setuptools import setup

setup(
    name="halogenase_detection",
    version="0.1.0",
    description="A Python package to identify categorize and mine for halogenases",
    url="https://github.com/JudSze/hal_miner",
    author="Judit Szenei",
    author_email="szenei@dtu.dk",
    license="Affero",
    packages=["miner"],
    install_requires=["pyhmmer", "pandas"]
)