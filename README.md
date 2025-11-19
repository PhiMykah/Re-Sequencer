# Re-Sequencer

The **Re-Sequencer** project is designed to streamline and automate the modification of PDBs for use in molecular dynamics (MD) simulations. It provides tools for residue substitution, and resdiue extension of proteins, making it easier for researchers to test modified PDBs. The project supports various input formats and integrates with popular MD analysis libraries to ensure compatibility and extensibility.

## Key Features

- Automated re-sequencing of MD trajectory frames
- User-friendly command-line interface
- Modular and extensible codebase
- Python 3.10 Support 

## Building from Source

Set up a Python virtual environment using either `venv` or `conda` to manage dependencies and isolate your development environment.

### Requirements

Re-Sequencer requires the following packages/applications:
- `biopandas` (with `pandas` and `numpy`)
- `x3DNA`
- `pymol`
- `Linux/Mac` (windows not yet supported)

### Installing Dependencies

It is recommended that you use [`conda`](https://www.anaconda.com/docs/getting-started/miniconda/install) to install the dependecies for resequencer. 

#### Manual `pymol` installation with venv
If you choose to use [*venv*](https://docs.python.org/3/library/venv.html), you will need to install `pymol` through the website [here](https://www.pymol.org/). Pymol must be added to your path and able to run through the command line. This *has not* been thoroughly tested, so proceed at your own risk.

### Installation

1. Either clone the repository with the following line:
    ```bash
    git clone https://github.com/PhiMykah/Re-Sequencer
    ```
    or [download the main zip file](https://github.com/PhiMykah/Re-Sequencer/archive/refs/heads/main.zip) on Github and extract the Re-Sequencer Folder.

2. After installing `conda`, execute the following code, making sure to type `y` when prompted:
    ```bash
    conda create -n resequencer python=3.10
    conda activate resequencer
    conda config --append channels conda-forge
    conda install biopandas
    conda install -c conda-forge -c schrodinger pymol-bundle
    ```

3. Run the following command to install Re-Sequencer to your environment:
    ```bash
    cd Re-Sequencer # or cd Re-Sequencer-main
    pip install -e .
    ```
    This installs the current directory as a python package, and the `-e` flag marks it as editable. This means if you update Re-Sequencer in the future using `git pull` or `git pull --rebase`, it will update the package automatically

4. Test that **Re-Sequencer** has been successfully installed:
    ```bash
    resequencer --help
    ```
    It is recommended you do excute this command in another directory to ensure it is not only working in the main directory.

After activating your environment, you can manage the project without affecting your global Python installation.
