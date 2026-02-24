from dataclasses import dataclass
from typing import Optional, Tuple


@dataclass
class InterfaceConfig:
    # matching
    max_sub_miller: int = 2
    max_film_miller: int = 2
    # misfit_limit: float = 5.0  # %
    von_mises_limit: float = 0.05

    sub_miller: Optional[Tuple[int, int, int]] = None
    film_miller: Optional[Tuple[int, int, int]] = None

    max_supercell: int = 10  # <-- добавлено

    # interface construction
    num_sites_limit: int = 300
    density_limit: float = 0.01
    charge_limit: float = 0.2

    gap: float = 2.0
    vacuum_over_film: float = 15.0
    film_thickness: int = 10
    substrate_thickness: int = 10
    surface_thickness: float = 0.8

    output_folder: str = "interfaces"
