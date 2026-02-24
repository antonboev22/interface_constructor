import logging
from pathlib import Path
from typing import Dict, Any, List
from multiprocessing import Pool, cpu_count
import csv

from pymatgen.io.vasp import Poscar
from pymatgen.analysis.interfaces.zsl import ZSLGenerator
from pymatgen.analysis.interfaces.coherent_interfaces import CoherentInterfaceBuilder

from .analysis import compute_surface_density, compute_surface_charge_density
import re

logger = logging.getLogger(__name__)


# ============================================================
# I/O
# ============================================================

def _write_interface_files(interface, substrate, film, out_dir: Path, filename: str):
    out_dir.mkdir(parents=True, exist_ok=True)
    Poscar(interface).write_file(out_dir / f"{filename}.POSCAR")
    Poscar(substrate).write_file(out_dir / f"{filename}_substrate.POSCAR")
    Poscar(film).write_file(out_dir / f"{filename}_film.POSCAR")

def _estimate_interface_sites(match, substrate_bulk, film_bulk, config):
    """
    Rough estimate of number of atoms in the interface supercell.
    Works BEFORE building the interface.
    """
    try:
        sc_sub = match.get("supercell_substrate", (1, 1))
        sc_film = match.get("supercell_film", (1, 1))

        mult_sub = abs(int(sc_sub[0]) * int(sc_sub[1]))
        mult_film = abs(int(sc_film[0]) * int(sc_film[1]))

        n_sub = len(substrate_bulk)
        n_film = len(film_bulk)

        est = (
            mult_sub * n_sub * config.substrate_thickness
            + mult_film * n_film * config.film_thickness
        )

        return est

    except Exception:
        return 0




def safe_termination_str(termination) -> str:
    """
    Преобразует объект терминации в строку для имени файла:
    убирает все небезопасные символы, оставляет только буквы, цифры и _
    """
    if isinstance(termination, tuple):
        term_str = "_".join(map(str, termination))
    else:
        term_str = str(termination)
    # заменяем все, что не буквы/цифры/_ на _
    term_str = re.sub(r"[^0-9a-zA-Z_]+", "_", term_str)
    return term_str

# ============================================================
# WORKER: один match
# ============================================================

def _process_match(args):
    substrate_bulk, film_bulk, match_record, config, output_dir = args

    results_all = []
    accepted_all = []

    sub_hkl = match_record["hkl_sub"]
    film_hkl = match_record["hkl_film"]

    try:
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

    except Exception as e:
        logger.warning(f"Failed to build CIB for {sub_hkl}/{film_hkl}: {e}")
        return [], []

    sub_formula = substrate_bulk.composition.reduced_formula
    film_formula = film_bulk.composition.reduced_formula

    # ========================================================
    # LOOP по терминациям
    # ========================================================

    for termination in cib.terminations:

        # --- EARLY SIZE REJECT (BIG SPEEDUP) ---
        est_sites = _estimate_interface_sites(
            match_record, substrate_bulk, film_bulk, config
        )

        if est_sites > config.num_sites_limit * 1.3:
            logger.debug(
                f"Skip large match early: est {est_sites} atoms"
            )
            continue

        # собираем все интерфейсы для этой терминации
        term_interfaces = [
            interface
            for interface in cib.get_interfaces(
                termination=termination,
                gap=config.gap,
                vacuum_over_film=config.vacuum_over_film,
                film_thickness=config.film_thickness,
                substrate_thickness=config.substrate_thickness,
                in_layers=False,
            )
            if interface.num_sites <= config.num_sites_limit
        ]

        if not term_interfaces:
            continue

        # --- физика ---
        filtered_interfaces = []
        for interface in term_interfaces:
            sub_density = compute_surface_density(
                interface.substrate,
                select="top",
                layer_thickness=config.surface_thickness,
            )

            film_density = compute_surface_density(
                interface.film,
                select="bottom",
                layer_thickness=config.surface_thickness,
            )

            sub_charge = compute_surface_charge_density(
                interface.substrate,
                select="top",
                layer_thickness=config.surface_thickness,
            )

            film_charge = compute_surface_charge_density(
                interface.film,
                select="bottom",
                layer_thickness=config.surface_thickness,
            )

            total_charge = abs(sub_charge + film_charge)

            passed = (
                sub_density > config.density_limit
                and film_density > config.density_limit
                and total_charge < config.charge_limit
            )

            if passed:
                filtered_interfaces.append(
                    (interface, sub_density, film_density, sub_charge, film_charge, total_charge)
                )
            else:
                logger.warning(f"Density or charge condition has not passed")


        if not filtered_interfaces:
            continue

        # ========================================================
        # берем интерфейс с минимальным num_sites для этой терминации
        # ========================================================

        interface, sub_density, film_density, sub_charge, film_charge, total_charge = min(
            filtered_interfaces,
            key=lambda x: x[0].num_sites
        )

        term_str = safe_termination_str(termination)

        filename = (
            f"{sub_formula}_{film_formula}_"
            f"{''.join(map(str, sub_hkl))}_"
            f"{''.join(map(str, film_hkl))}_"
            f"{term_str}_"
            f"{interface.num_sites}at"
        )

        out_subdir = (
            Path(output_dir)
            / f"{film_formula}_on_{sub_formula}"
            / f"{''.join(map(str, sub_hkl))}_{''.join(map(str, film_hkl))}"
        )

        _write_interface_files(
            interface,
            interface.substrate,
            interface.film,
            out_subdir,
            filename,
        )

        results_all.append(
            dict(
                substrate=sub_formula,
                film=film_formula,
                hkl_sub=sub_hkl,
                hkl_film=film_hkl,
                termination=str(termination),
                num_sites=interface.num_sites,
                substrate_density=round(sub_density, 4),
                film_density=round(film_density, 4),
                substrate_charge_density=round(sub_charge, 4),
                film_charge_density=round(film_charge, 4),
                abs_charge_density=round(total_charge, 4),
                von_mises=match_record["von_mises"],
                slab=filename,
            )
        )

        accepted_all.append(interface)

    return results_all, accepted_all


# ============================================================
# MAIN
# ============================================================

def build_interfaces(
    substrate_bulk,
    film_bulk,
    match_record,
    # matches: List[Dict[str, Any]],
    config,
    output_dir: str = "interfaces",
    ):
    """
    ULTRA-FAST version:
    multiprocessing по MATCHES (самый эффективный вариант)
    """
    if isinstance(match_record, dict):
        class _Match:
            pass
        match = _Match()
        for k, v in match_record.items():
            setattr(match, k, v)
    else:
        match = match_record

    sub_hkl = match_record["hkl_sub"]
    film_hkl = match_record["hkl_film"]
    von_mises = match_record["von_mises"]

    logger.info(
        f"Processing match: {sub_hkl}/{film_hkl}, "
        f"von Mises = {von_mises:.2f}%"
    )



   


    args_list = [(substrate_bulk, film_bulk, match_record, config, output_dir)]

    # ========================================================
    # multiprocessing по матчам
    # ========================================================

    nproc = min(cpu_count(), len(args_list))

    if nproc == 1:
        pool_results = [_process_match(args_list[0])]
    else:
        from multiprocessing import Pool
        with Pool(processes=nproc) as pool:
            pool_results = pool.map(_process_match, args_list)

    # ========================================================
    # collect
    # ========================================================

    results_all = []
    accepted_all = []

    for res, interfaces in pool_results:
        results_all.extend(res)
        accepted_all.extend(interfaces)

    # ========================================================
    # deduplicate: keep minimal atoms per orientation+termination
    # ========================================================

    def _dedup_interfaces(results, interfaces):
        """
        Keep only the smallest interface for each
        (hkl_sub, hkl_film, termination_sub, termination_film)
        """
        best = {}

        for r, iface in zip(results, interfaces):
            key = (
                tuple(r["hkl_sub"]),
                tuple(r["hkl_film"]),
                r.get("termination"),
            )

            n_atoms = r["num_sites"]

            if key not in best or n_atoms < best[key][0]["num_sites"]:
                best[key] = (r, iface)

        new_results = []
        new_interfaces = []

        for r, iface in best.values():
            new_results.append(r)
            new_interfaces.append(iface)

        return new_results, new_interfaces


    before = len(results_all)
    results_all, accepted_all = _dedup_interfaces(results_all, accepted_all)
    after = len(results_all)

    logger.info(f"Deduplicated interfaces: {before} → {after}")

    # ========================================================
    # ranking
    # ========================================================

    results_all.sort(
        key=lambda x: (
            x["abs_charge_density"],
            -x["num_sites"],
            x["von_mises"],
        )
    )

    # ========================================================
    # CSV
    # ========================================================

    if results_all:
        csv_path = (
            Path(output_dir)
            / f"{substrate_bulk.composition.reduced_formula}_"
              f"{film_bulk.composition.reduced_formula}_summary.csv"
        )

        csv_path.parent.mkdir(parents=True, exist_ok=True)

        with open(csv_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=list(results_all[0].keys()))
            writer.writeheader()
            writer.writerows(results_all)

        logger.info(f"CSV summary saved: {csv_path}")
        logger.info(f"Total accepted interfaces: {len(results_all)}")
    else:
        logger.warning("No interfaces passed filters — CSV not written.")

    return results_all, accepted_all

