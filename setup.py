import setuptools

with open("README.md") as fh:
    long_description=fh.read()

setuptools.setup(
    name="ngs-primerplex",
    version="1.1",
    author="Andrey Kechin",
    author_email="aa_kechin@niboch.nsc.ru",
    description="High-throughput tool for mupltiplex primer design",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/aakechin/NGS-PrimerPlex",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: Unix",
        ],
    python_requires=">=3.0",
)
