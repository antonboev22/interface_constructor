## Contents
- [About](#about)
- [Installation](#installation)
- [Examples](#examples)


### About

Interface Constructor is a Python library for automated construction and analysis of coherent interfaces between two crystalline solids.

The package provides:

lattice matching based on Miller indices

strain evaluation using the von Mises criterion

automated slab and interface generation

physical filtering (surface density, charge balance, size limits)

high-throughput–ready workflow

It is designed for computational materials science studies, including thin films, coatings, and solid–solid interfaces.

### Installation

```
git clone https://github.com/cmg-storion/interface_constructor
cd interface_constructor
pip install .
```

### Examples

Example workflows and notebooks are provided in the examples/ directory.
<!-- - [UMLIP-NEB percolation analysis](examples/umlip_neb_percolation_analysis.ipynb) -->
<!-- - [Li|Anode coating|Electrolyte|Cathode sandwich stability](examples/sandwich_stability.ipynb) (using only the Materials Project's entries) -->
<!-- - [High-throughput pipeline for screening the MP database](examples/workflow_MP.ipynb) -->

Typical usage:

from pymatgen.core import Structure
from interface_constructor import InterfaceConstructor, InterfaceConfig

substrate = Structure.from_file("substrate.cif")
film = Structure.from_file("film.cif")

config = InterfaceConfig()

ic = InterfaceConstructor(substrate, film, config)

matches = ic.run_matching()
metadata, interfaces = ic.build_all_interfaces()
ic.save_summary_csv()