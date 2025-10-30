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

### Using `venv`

```bash
git clone https://github.com/PhiMykah/Re-Sequencer
cd Re-Sequencer
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate
pip install -r requirements.txt
pip install -e . 
```

### Using `conda`

```bash
git clone https://github.com/PhiMykah/Re-Sequencer
conda env create -f resequencer_env.yml
conda activate resequencer
pip install -e . 
```

After activating your environment, you can manage the project without affecting your global Python installation.