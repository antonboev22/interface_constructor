import logging
from pathlib import Path
from typing import Dict, Any, List
from pymatgen.io.vasp import Poscar
from pymatgen.analysis.interfaces.zsl import ZSLGenerator
from pymatgen.analysis.interfaces.coherent_interfaces import CoherentInterfaceBuilder
from .analysis import compute_surface_density, compute_surface_charge_density
import re

logger = logging.getLogger(__name__)

# -------------------- helper functions --------------------

def safe_termination_str(termination) -> str:
    """Convert termination into a filesystem-safe string."""
    term_str = "_".join(map(str, termination)) if isinstance(termination, tuple) else str(termination)
    return re.sub(r"[^0-9a-zA-Z_]+", "_", term_str)


def _write_interface_files(interface, substrate, film, out_dir: Path, filename: str):
    """Write interface and slab structures to POSCAR files."""
    out_dir.mkdir(parents=True, exist_ok=True)
    Poscar(interface).write_file(out_dir / f"{filename}.POSCAR")
    Poscar(substrate).write_file(out_dir / f"{filename}_substrate.POSCAR")
    Poscar(film).write_file(out_dir / f"{filename}_film.POSCAR")


def save_csv(metadata, substrate, film, output_dir: str):
    """Save CSV summary of generated interfaces."""
    import csv

    if not metadata:
        logger.warning("No metadata to save — CSV skipped.")
        return

    csv_path = (
        Path(output_dir)
        / f"{substrate.composition.reduced_formula}_{film.composition.reduced_formula}_summary.csv"
    )
    csv_path.parent.mkdir(parents=True, exist_ok=True)

    with open(csv_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(metadata[0].keys()))
        writer.writeheader()
        writer.writerows(metadata)

    logger.info(f"CSV summary saved: {csv_path}")


# -------------------- core processing --------------------

def _process_match(substrate_bulk, film_bulk, match_record, config, output_dir: str, logger):
    results_all = []
    accepted_all = []

    sub_hkl = match_record["hkl_sub"]
    film_hkl = match_record["hkl_film"]

    zsl = ZSLGenerator(
        max_area=400,
        max_area_ratio_tol=0.05,
        max_length_tol=0.05,
        max_angle_tol=1,
        bidirectional=False,
    )

    cib = CoherentInterfaceBuilder(
        film_structure=film_bulk,
        substrate_structure=substrate_bulk,
        film_miller=film_hkl,
        substrate_miller=sub_hkl,
        zslgen=zsl,
        filter_out_sym_slabs=True,
    )

    sub_formula = substrate_bulk.composition.reduced_formula
    film_formula = film_bulk.composition.reduced_formula

    for termination in cib.terminations:
        logger.info(f"  Termination: {termination}")

        term_interfaces = [
            iface
            for iface in cib.get_interfaces(
                termination=termination,
                gap=config.gap,
                vacuum_over_film=config.vacuum_over_film,
                film_thickness=config.film_thickness,
                substrate_thickness=config.substrate_thickness,
                in_layers=False,
            )
            if iface.num_sites <= config.num_sites_limit
        ]

        if not term_interfaces:
            continue

        # Physical filtering
        filtered_interfaces = []
        for iface in term_interfaces:
            sub_density = compute_surface_density(
                iface.substrate,
                select="top",
                layer_thickness=config.surface_thickness,
            )
            film_density = compute_surface_density(
                iface.film,
                select="bottom",
                layer_thickness=config.surface_thickness,
            )

            sub_charge = compute_surface_charge_density(
                iface.substrate,
                select="top",
                layer_thickness=config.surface_thickness,
            )
            film_charge = compute_surface_charge_density(
                iface.film,
                select="bottom",
                layer_thickness=config.surface_thickness,
            )

            total_charge = abs(sub_charge + film_charge)

            passed = (
                sub_density > config.density_limit
                and film_density > config.density_limit
                and total_charge < config.charge_limit
            )

            filtered_interfaces.append(
                (
                    iface,
                    sub_density,
                    film_density,
                    sub_charge,
                    film_charge,
                    total_charge,
                    passed,
                )
            )

        if not filtered_interfaces:
            continue

        # Select interface with minimal number of atoms
        iface, sub_density, film_density, sub_charge, film_charge, total_charge, passed = min(
            filtered_interfaces,
            key=lambda x: x[0].num_sites,
        )

        term_str = safe_termination_str(termination)

        filename = (
            f"{sub_formula}_{film_formula}_"
            f"{''.join(map(str, sub_hkl))}_{''.join(map(str, film_hkl))}_"
            f"{term_str}_{iface.num_sites}at"
        )

        out_subdir = (
            Path(output_dir)
            / f"{film_formula}_on_{sub_formula}"
            / f"{''.join(map(str, sub_hkl))}_{''.join(map(str, film_hkl))}"
        )

        _write_interface_files(iface, iface.substrate, iface.film, out_subdir, filename)

        results_all.append(
            {
                "substrate": sub_formula,
                "film": film_formula,
                "hkl_sub": sub_hkl,
                "hkl_film": film_hkl,
                "termination": str(termination),
                "num_sites": iface.num_sites,
                "substrate_density": round(sub_density, 4),
                "film_density": round(film_density, 4),
                "substrate_charge_density": round(sub_charge, 4),
                "film_charge_density": round(film_charge, 4),
                "abs_charge_density": round(total_charge, 4),
                "von_mises": match_record["von_mises"],
                "slab": filename,
                "passed_filters": passed,
            }
        )

        if passed:
            accepted_all.append(iface)

    return results_all, accepted_all


# -------------------- main builder --------------------

def build_interfaces(
    substrate_bulk,
    film_bulk,
    match_record,
    config,
    output_dir="interfaces",
    logger=None,
):
    """Build interfaces for a single match record."""
    if logger is None:
        logger = logging.getLogger(__name__)

    results_all, accepted_all = _process_match(
        substrate_bulk,
        film_bulk,
        match_record,
        config,
        output_dir,
        logger,
    )

    # Sorting priority:
    # 1) lowest total charge
    # 2) largest num_sites
    # 3) lowest von Mises strain
    results_all.sort(
        key=lambda x: (x["abs_charge_density"], -x["num_sites"], x["von_mises"])
    )

    # Save CSV
    save_csv(results_all, substrate_bulk, film_bulk, output_dir)

    return results_all, accepted_all