import numpy as np


def compute_surface_density(structure, select="top", layer_thickness=1.0):
    """Surface atomic density (atoms / Å²)."""
    a_vector = structure.lattice.matrix[0]
    b_vector = structure.lattice.matrix[1]
    surface_area = np.linalg.norm(np.cross(a_vector, b_vector))

    z_coords = np.array([site.coords[2] for site in structure])
    z_max, z_min = z_coords.max(), z_coords.min()

    if select == "top":
        surface_atoms = [s for s in structure if (z_max - s.coords[2] <= layer_thickness)]
    elif select == "bottom":
        surface_atoms = [s for s in structure if (s.coords[2] - z_min <= layer_thickness)]
    else:
        raise ValueError("select must be 'top' or 'bottom'")

    return len(surface_atoms) / surface_area


def compute_surface_charge_density(structure, select="top", layer_thickness=1.0):
    """Surface charge density (e / Å²)."""

    a_vector = structure.lattice.matrix[0]
    b_vector = structure.lattice.matrix[1]
    surface_area = np.linalg.norm(np.cross(a_vector, b_vector))

    z_coords = np.array([site.coords[2] for site in structure])
    z_max, z_min = z_coords.max(), z_coords.min()

    if select == "top":
        surface_atoms = [s for s in structure if (z_max - s.coords[2] <= layer_thickness)]
    elif select == "bottom":
        surface_atoms = [s for s in structure if (s.coords[2] - z_min <= layer_thickness)]
    else:
        raise ValueError("select must be 'top' or 'bottom'")

    try:
        charge = sum(site.specie.oxi_state for site in surface_atoms)
    except AttributeError as e:
        raise ValueError(
            "Structure must have oxidation states assigned."
        ) from e

    return charge / surface_area