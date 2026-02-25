import numpy as np


def _select_surface_atoms(structure, select="top", layer_thickness=1.0):
    z_coords = np.array([s.coords[2] for s in structure])
    z_max, z_min = z_coords.max(), z_coords.min()

    if select == "top":
        return [s for s in structure if (z_max - s.coords[2] <= layer_thickness)]
    elif select == "bottom":
        return [s for s in structure if (s.coords[2] - z_min <= layer_thickness)]
    else:
        raise ValueError("select must be 'top' or 'bottom'")
    
def compute_surface_density(structure, select="top", layer_thickness=1.0):
    a_vector, b_vector = structure.lattice.matrix[:2]
    surface_area = np.linalg.norm(np.cross(a_vector, b_vector))
    surface_atoms = _select_surface_atoms(structure, select, layer_thickness)
    return len(surface_atoms) / surface_area

def compute_surface_charge_density(structure, select="top", layer_thickness=1.0):
    a_vector, b_vector = structure.lattice.matrix[:2]
    surface_area = np.linalg.norm(np.cross(a_vector, b_vector))
    surface_atoms = _select_surface_atoms(structure, select, layer_thickness)
    try:
        charge = sum(site.specie.oxi_state for site in surface_atoms)
    except AttributeError:
        raise ValueError("Structure must have oxidation states assigned.")
    return charge / surface_area