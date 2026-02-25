import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[1]))

from pymatgen.core import Structure
from interface_constructor import InterfaceConstructor, InterfaceConfig

# -------------------------------
# Read substrate structure
# -------------------------------
substrate_file = "LiCoO2.cif"
substrate = Structure.from_file(substrate_file)
substrate.add_oxidation_state_by_element({
    "Li": +1,
    "O": -2,
    "Co": +3,
})

# -------------------------------
# Read film structure
# -------------------------------
film_file = "Li.cif"
film = Structure.from_file(film_file)
film.add_oxidation_state_by_element({
    "Li": +1,
})



# -------------------------------
# Configure interface search
# -------------------------------
config = InterfaceConfig(
    # Miller indices control
    sub_miller=(1, 0, 4),     # fixed substrate Miller index (optional)
    # film_miller=(1, 1, 1),   # fixed film Miller index (optional)

    # max_sub_miller=1,         # used if sub_miller=None
    max_film_miller=1,        # used if film_miller=None

    # Matching filters
    von_mises_limit=0.05,     # maximum allowed strain
    num_sites_limit=200,      # maximum number of atoms

    # Optional filters
    # density_limit=0.01,
    # charge_limit=0.2,

    # Geometry parameters
    film_thickness=4,
    substrate_thickness=4,
    surface_thickness=0.8,
    gap=2.0,
    vacuum_over_film=15.0,

    # Output
    output_folder="interfaces",
)

# -------------------------------
# Initialize constructor
# -------------------------------
ic = InterfaceConstructor(substrate, film, config)

# -------------------------------
# 1) Find lattice matches
# -------------------------------

'''
    What happens here:
    - scans Miller indices
    - evaluates lattice compatibility
    - filters by von Mises strain
    - removes duplicates 
'''
matches = ic.run_matching()
print(f"Found {len(matches['filtered'])} suitable matches.")

# -------------------------------
# 2) Build interfaces
# -------------------------------
'''
    This step:
    - constructs slab models
    - aligns film and substrate
    - applies gap and vacuum
    - saves structures to disk
'''
metadata, structures = ic.build_all_interfaces(
    output_dir=config.output_folder
)
print(f"Built {len(metadata)} interfaces.")

# -------------------------------
# 3) Save CSV summary
# -------------------------------
ic.save_summary_csv(output_dir=config.output_folder)
print("CSV summary saved.")
