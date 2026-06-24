# MyScripts - Bioinformatics, ML & Simulations Toolkit

A comprehensive collection of Python scripts for protein bioinformatics, machine learning (AlphaFold3/PyRosetta), and molecular dynamics simulations. This toolkit provides reusable modules for sequence analysis, structure prediction, and computational chemistry workflows.

[![Python 3.8+](https://img.shields.io/badge/Python-3.8+-blue.svg)](https://www.python.org/)
[![License](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

---

## 📋 Table of Contents

- [Overview](#overview)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Project Structure](#project-structure)
- [Modules](#modules)
- [Configuration](#configuration)
- [Usage Examples](#usage-examples)
- [Development](#development)
- [Contributing](#contributing)
- [License](#license)

---

## 🔍 Overview

**MyScripts** is a modular toolkit designed for researchers working with:

- **Protein Sequence Analysis**: MSA processing, alignment refinement, conservation scoring
- **Protein Structure Prediction**: AlphaFold3 input preparation, output processing
- **Protein Manipulation**: Structure refinement, energy minimization (PyRosetta)
- **Molecular Dynamics**: GROMACS topology management, trajectory analysis
- **Data Visualization**: Interactive viewers, plotting utilities, image processing

### Key Features

✅ Production-ready bioinformatics algorithms  
✅ Integration with industry-standard tools (GROMACS, PyRosetta, AlphaFold3)  
✅ Modular design for easy integration into workflows  
✅ Comprehensive documentation and examples  
✅ Configuration-driven (no hardcoded paths)  
✅ Type hints and docstrings throughout  

---

## 💾 Installation

### Prerequisites

- Python 3.8 or higher
- pip or conda

### Basic Installation

```bash
# Clone the repository
git clone https://github.com/yourusername/myScripts.git
cd myScripts

# Install in development mode
pip install -e .

# Or install with all dependencies
pip install -e ".[dev]"
```

### Install with ML Support

For machine learning features (PyRosetta):

```bash
pip install -e ".[ml]"
```

### Install with Simulation Support

For MD simulation features:

```bash
pip install -e ".[simulations]"
```

### Development Installation

For development with testing and documentation tools:

```bash
pip install -r requirements-dev.txt
```

### Using Conda (Optional)

```bash
conda create -n myscripts python=3.10
conda activate myscripts
pip install -e .
```

---

## ⚡ Quick Start

### 1. Configuration Setup

Create a `config.yaml` file from the example:

```bash
cp config.example.yaml config.yaml
```

Edit `config.yaml` with your paths:

```yaml
paths:
  data_root: /path/to/your/data
  alphafold_db: /path/to/alphafold/databases
  gromacs_templates: /path/to/gromacs/templates

tools:
  jackhmmer_path: /usr/bin/jackhmmer
  gmx_path: /usr/bin/gmx
```

### 2. Basic Usage

```python
from myScripts.bioinformatics.sequence import seq_manager
from myScripts.config import Config

# Load configuration
config = Config('config.yaml')

# Convert FASTA format
seq_manager.multi2SingleLineFasta('input.fasta', 'output.fasta')

# Convert A3M to FASTA
seq_manager.a3mtoFasta('alignment.a3m', 'alignment.fasta')
```

### 3. Process Alignment

```python
from myScripts.bioinformatics.sequence import alignment

# View and analyze alignment
alignment.main(
    aln_file='alignment.fasta',
    query_id='NP_001317000'
)
```

### 4. Prepare AlphaFold3 Inputs

```python
from myScripts.machineLearning.alphaFold3 import inputPrepAF3

# Generate AF3 input JSON
inputPrepAF3.generateAF3InputJson(
    name='my_protein',
    sequence='MKVLWAALLVTFLAGCAKAKAV...',
    model_seeds=[1, 2, 3],
    output_file='af3_input.json',
    ligand='ATP'
)
```

---

## 📁 Project Structure

```
myScripts/
├── README.md                          # This file
├── PROJECT_REVIEW.md                  # Comprehensive code review
├── CONTRIBUTING.md                    # Contributing guidelines
├── setup.py                           # Package setup
├── requirements.txt                   # Core dependencies
├── requirements-dev.txt               # Development dependencies
├── config.example.yaml                # Example configuration
├── .gitignore                         # Git ignore rules
│
├── myScripts/                         # Main package
│   ├── __init__.py
│   ├── config.py                      # Configuration management
│   │
│   ├── bioinformatics/                # Bioinformatics module
│   │   ├── __init__.py
│   │   ├── clustalOutputReformat.py   # CLUSTAL output processing
│   │   ├── esbgTools.tcl              # TCL utilities
│   │   ├── sequence/                  # Sequence analysis
│   │   │   ├── __init__.py
│   │   │   ├── alignment.py           # MSA analysis
│   │   │   ├── seq_manager.py         # Format conversion
│   │   │   ├── seq_clustering.py      # Sequence clustering
│   │   │   ├── msaCleaning.py         # MSA refinement
│   │   │   └── uniProtApi.py          # UniProt integration
│   │   ├── pymol/                     # PyMOL scripts
│   │   │   ├── imageCartoon.py
│   │   │   ├── imagePutty.py
│   │   │   ├── mutateProtein.py
│   │   │   └── orientKinase.py
│   │   └── hexDock/                   # HexDock integration
│   │       ├── hexdock.sh
│   │       └── hexdock_template.mac
│   │
│   ├── machineLearning/               # ML module
│   │   ├── __init__.py
│   │   ├── alphaFold3/                # AlphaFold3 tools
│   │   │   ├── __init__.py
│   │   │   ├── AF3MSA.py              # MSA generation
│   │   │   ├── inputPrepAF3.py        # Input preparation
│   │   │   ├── mma2a3m.py             # MSA format conversion
│   │   │   ├── outputAF3.py           # Result processing
│   │   │   └── UniProtAPI.py          # UniProt queries
│   │   └── pyrosetta/                 # PyRosetta tools
│   │       ├── freeEnePyrosetta.py    # Free energy calculations
│   │       └── pyrosettminimization.py # Structure minimization
│   │
│   ├── simulations/                   # MD simulations module
│   │   ├── __init__.py
│   │   ├── groTopManager.py           # GROMACS topology
│   │   ├── createMappedPDB.py         # Atom mapping
│   │   ├── createTopology.py          # Topology generation
│   │   ├── createAtomtypeMapping.py   # Atom type mapping
│   │   ├── simSetupCheck.py           # Setup validation
│   │   ├── gmx/                       # GROMACS utilities
│   │   │   └── simulationSetup/
│   │   ├── mdAnalysis/                # Trajectory analysis
│   │   │   ├── __init__.py
│   │   │   ├── aggregateCalculation.py
│   │   │   ├── centroidQresRMSD.py
│   │   │   └── RMSF.py
│   │   ├── trajectory/                # Trajectory tools
│   │   │   └── trajClustering.py
│   │   └── vmdTcl/                    # VMD visualization
│   │       └── motifReader.tcl
│   │
│   └── utils/                         # Utility module
│       ├── __init__.py
│       ├── basicFunction.py           # Common functions
│       ├── plotStyle.py               # Matplotlib utilities
│       ├── renameFiles.py             # File operations
│       ├── imageMovie/                # Image processing
│       │   ├── imageFunctions.py
│       │   ├── imageToGiforVid.py
│       │   └── animation.py
│       └── arialFonts/                # Font resources
│
├── docs/                              # Documentation
│   ├── installation.md
│   ├── quickstart.md
│   └── modules/
│
├── tests/                             # Unit tests
│   ├── __init__.py
│   ├── test_seq_manager.py
│   └── test_alignment.py
│
└── scripts/                           # Example scripts
    ├── process_alignment.py
    └── prepare_af3.py
```

---

## 📦 Modules

### Bioinformatics (`myScripts.bioinformatics`)

**Sequence Analysis and Processing**

- **alignment.py**: Query-centric MSA viewing, conservation scoring, difference highlighting
- **seq_manager.py**: Format conversion (FASTA/A3M/CMA), structure refinement
- **seq_clustering.py**: Sequence clustering and classification
- **msaCleaning.py**: MSA quality assessment and refinement
- **uniProtApi.py**: UniProt database queries

**Structure Visualization**

- **pymol/**: PyMOL scripts for structure visualization and manipulation
- **hexDock/**: Protein-protein docking integration

### Machine Learning (`myScripts.machineLearning`)

**AlphaFold3 Integration**

- **inputPrepAF3.py**: Generate AlphaFold3 input JSON configurations
- **AF3MSA.py**: MSA generation with JackHMMER
- **mma2a3m.py**: Multiple MSA format conversions
- **outputAF3.py**: Parse and analyze AlphaFold3 predictions
- **UniProtAPI.py**: Query UniProt for sequence databases

**PyRosetta (Energy Minimization)**

- **pyrosettminimization.py**: Structure refinement and optimization
- **freeEnePyrosetta.py**: Free energy calculations

### Simulations (`myScripts.simulations`)

**GROMACS Workflow**

- **groTopManager.py**: Topology file editing and generation
- **createMappedPDB.py**: Atom name and order mapping
- **createTopology.py**: Topology file creation
- **createAtomtypeMapping.py**: Atom type assignment

**Trajectory Analysis**

- **mdAnalysis/**: Analysis tools for MD trajectories
  - RMSF calculation
  - RMSD centroid analysis
  - Secondary structure assignment
  - Maxwell-Boltzmann distribution

**Visualization**

- **vmdTcl/**: VMD/TCL scripts for visualization

### Utilities (`myScripts.utils`)

- **basicFunction.py**: Common utility functions
- **plotStyle.py**: Matplotlib functionality (fonts, styles)
- **renameFiles.py**: Batch file operations
- **imageMovie/**: Image processing and animation

---

## ⚙️ Configuration

Create a `config.yaml` file (see `config.example.yaml`):

```yaml
paths:
  # Root directory for data
  data_root: /path/to/data
  
  # AlphaFold3 database location
  alphafold_db: /path/to/alphafold3/database
  
  # GROMACS template files
  gromacs_templates: /path/to/gromacs/templates
  
  # Output directories
  output_dir: ./outputs
  temp_dir: /tmp/myscripts

tools:
  # Tool paths
  jackhmmer_path: /usr/bin/jackhmmer
  gmx_path: /usr/bin/gmx
  pymol_path: /opt/pymol
  
logging:
  level: INFO
  file: logs/myscripts.log
```

Access configuration in code:

```python
from myScripts.config import Config

config = Config('config.yaml')
data_root = config.get('paths.data_root')
```

---

## 💡 Usage Examples

### Example 1: Process Multiple Sequence Alignment

```python
from myScripts.bioinformatics.sequence import alignment
from Bio import AlignIO

# View alignment with query highlighted
alignment.main(
    aln_file='sequences.fasta',
    query_id='NP_001317000'
)

# Get conservation scores
aln = AlignIO.read('sequences.fasta', 'fasta')
conservation = alignment.compute_conservation(aln)
print("Conservation:", conservation)
```

### Example 2: Prepare AlphaFold3 Input

```python
from myScripts.machineLearning.alphaFold3 import inputPrepAF3

# Generate input with custom MSAs
inputPrepAF3.generateAF3InputJson(
    name='my_kinase',
    sequence='MTKLTAVLLSVALVAFFVFCGSSDDS...',
    model_seeds=[1, 2, 3, 4, 5],
    unpaired_msa='# A3M format MSA',
    paired_msa=None,
    templates=None,
    ligand='ATP',
    output_file='af3_input.json'
)
```

### Example 3: GROMACS Topology Management

```python
from myScripts.simulations import groTopManager

# Edit topology to include ligand
groTopManager.edit_topology(
    top_path='topol.top',
    lig_itp='ATP.itp',
    metal_ion='MG',
    metal_ion_count=2
)

# Convert PDB to GRO with proper handling
groTopManager.pdbToGr0(
    pdb_file='protein.pdb',
    path='./',
    ligand='ATP'
)
```

### Example 4: Trajectory Analysis

```python
from myScripts.simulations.mdAnalysis import RMSF
import numpy as np

# Calculate RMSF
rmsf_values = RMSF.calculate_rmsf(
    trajectory='MD_trajectory.xtc',
    topology='topol.pdb'
)

# Get per-residue RMSF
residue_rmsf = np.mean(rmsf_values, axis=1)
```

---

## 🧪 Development

### Running Tests

```bash
# Run all tests
pytest

# With coverage
pytest --cov=myScripts

# Verbose output
pytest -v

# Specific test file
pytest tests/test_seq_manager.py
```

### Code Style & Quality

```bash
# Format code with Black
black myScripts/

# Check with Flake8
flake8 myScripts/

# Sort imports
isort myScripts/

# Type checking with mypy
mypy myScripts/

# Linting
pylint myScripts/
```

### Building Documentation

```bash
cd docs/
sphinx-build -b html . _build/html
```

---

## 🤝 Contributing

Contributions are welcome! Please see [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

### Development Workflow

1. Fork the repository
2. Create a feature branch: `git checkout -b feature/my-feature`
3. Make changes and test: `pytest`
4. Format code: `black . && isort .`
5. Commit with clear messages
6. Push to branch: `git push origin feature/my-feature`
7. Submit a pull request

---

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

## 📞 Contact & Support

- **Author**: Anup K. Prasad
- **Email**: anupkprasad121@gmail.com
- **Issues**: Open an issue on GitHub for bugs or feature requests

---

## 🙏 Acknowledgments

- BioPython community
- MDAnalysis developers
- AlphaFold team
- PyRosetta developers
- GROMACS community

---

## 📖 Additional Resources

- [PROJECT_REVIEW.md](PROJECT_REVIEW.md) - Comprehensive code review and improvements
- [Official BioPython Documentation](https://biopython.org/)
- [MDAnalysis Tutorials](https://userguide.mdanalysis.org/)
- [AlphaFold3 Documentation](https://alphafold.isab.ac.uk/)
- [GROMACS Manual](https://manual.gromacs.org/)

---

**Last Updated**: February 18, 2026  
**Current Version**: 1.0.0
