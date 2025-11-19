from setuptools import setup, Extension, find_packages
from Cython.Build import cythonize
from pathlib import Path
import sys

BASE_DIR = Path(__file__).resolve().parent

# Path to the fiber module
x3dna_dir = BASE_DIR / "resequencer" / "external" / "x3dna"

fiber_sources = [
    str(x3dna_dir / "ana_fncs.c"),
    str(x3dna_dir / "app_fncs.c"),
    str(x3dna_dir / "cmn_fncs.c"),
    str(x3dna_dir / "fncs_slre.c"),
    str(x3dna_dir / "nrutil.c"),
    str(x3dna_dir / "reb_fncs.c"),
    # str(x3dna_dir / "fiber.c"),
]
# [str(p) for p in x3dna_dir.glob("*.c")]
fiber_sources.append(str(x3dna_dir / "fiber_bindings.pyx"))

fiber_ext = Extension(
    "resequencer.external.x3dna.fiber_bindings",  # dotted module name
    sources=fiber_sources,
    include_dirs=[str(x3dna_dir)],
    define_macros=[
        ("_CRT_SECURE_NO_WARNINGS", None),
        ("X3DNA_DIR", f'"{str(x3dna_dir.as_posix())}"'),
    ]
    if sys.platform.startswith("win")
    else [],
    # extra_compile_args=[""],
)

setup(
    name="resequencer",
    version="0.2.0",
    packages=find_packages(),
    install_requires=["biopandas", "pandas", "Cython"],
    ext_modules=cythonize([fiber_ext]),
    entry_points={
        "console_scripts": [
            "resequencer = resequencer.main:main",
        ]
    },
    author="Micah Smith",
    author_email="mykahsmith21@gmail.com",
    description="Replace Residues in PDB files with a simple command-line call!",
)
