from setuptools import setup, find_packages

setup(
    name="resequencer",
    version="0.1.1",
    packages=find_packages(),
    install_requires=["biopandas", "pandas"],
    entry_points={
        "console_scripts": [
            "resequencer = resequencer.main:main",
        ]
    },
    author="Micah Smith",
    author_email="mykahsmith21@gmail.com",
    description="Replace Residues in PDB files with a simple command-line call!",
)
