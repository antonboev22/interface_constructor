import itertools
import numpy as np
from typing import List, Tuple, Dict, Any, Optional

from pymatgen.analysis.interfaces.substrate_analyzer import SubstrateAnalyzer

from collections import Counter

def log_matches_compact(matches, logger=None):
    """
    Log repeated hkl(sub/film) pairs with occurrence counting.
    matches: list of dicts with keys 'hkl_sub', 'hkl_film', 'von_mises'
    """
    # Build tuples (hkl_sub, hkl_film, von_mises) for counting
    counter = Counter(
        (tuple(m['hkl_sub']), tuple(m['hkl_film']), round(m['von_mises'], 4))
        for m in matches
    )

    # Print in compact form
    for (hkl_sub, hkl_film, von), count in counter.items():
        msg = f"hkl(sub/film) = {hkl_sub}/{hkl_film} | von Mises = {von*100:.1f}%"
        if count > 1:
            msg += f" × {count}"
        if logger:
            logger.info(msg)
        else:
            print(msg)


def generate_millers(max_index: int):
    """Generate Miller indices up to a given maximum index."""
    return [
        (h, k, l)
        for h, k, l in itertools.product(range(max_index + 1), repeat=3)
        if (h, k, l) != (0, 0, 0)
    ]


def _deduplicate_matches(matches, strain_tol=1e-4):
    """Remove duplicate orientation matches (robust version)."""
    unique = []
    seen = set()

    for m in matches:

        if isinstance(m, dict):
            hkl_sub = tuple(m.get("hkl_sub") or m.get("hkl_substrate"))
            hkl_film = tuple(m.get("hkl_film"))

            sc_sub = tuple(m.get("supercell_substrate", ()))
            sc_film = tuple(m.get("supercell_film", ()))

            strain = round(
                float(
                    m.get("von_mises")
                    or m.get("von_mises_strain", 0)
                ),
                4,
            )
        else:
            hkl_sub = tuple(getattr(m, "hkl_substrate"))
            hkl_film = tuple(getattr(m, "hkl_film"))

            sc_sub = tuple(getattr(m, "supercell_substrate", ()))
            sc_film = tuple(getattr(m, "supercell_film", ()))

            strain = round(
                float(getattr(m, "von_mises_strain", 0)),
                4,
            )

        key = (hkl_sub, hkl_film, sc_sub, sc_film, strain)

        if key not in seen:
            seen.add(key)
            unique.append(m)

    return unique


def find_matches(
    substrate,
    film,
    config,
    logger: Optional[Any] = None
) -> Dict[str, List[Dict[str, Any]]]:
    """Find all possible matches between substrate and film."""

    if logger:
        logger.info("Start matching...")

    analyzer = SubstrateAnalyzer(
        film_max_miller=config.max_film_miller,
        bidirectional=False,
        max_area_ratio_tol=0.09,
        max_area=400,
        max_length_tol=0.03,
        max_angle_tol=0.01,
    )

    results_all = []
    results_filtered = []

    # -------- define Miller sets depending on the mode -------- #
    sub_millers = (
        generate_millers(config.max_sub_miller)
        if config.sub_miller is None
        else [config.sub_miller]
    )

    film_millers_allowed = (
        generate_millers(config.max_film_miller)
        if config.film_miller is None
        else [config.film_miller]
    )

    mode_sub = "auto" if config.sub_miller is None else "fixed"
    mode_film = "auto" if config.film_miller is None else "fixed"

    if logger:
        logger.info(f"Matching mode: substrate={mode_sub}, film={mode_film}")

    for sub_hkl in sub_millers:

        matches = analyzer.calculate(
            film=film,
            substrate=substrate,
            substrate_millers=[sub_hkl],
        )

        for match in matches:
            film_hkl = match.film_miller
            if film_hkl not in film_millers_allowed:
                continue

            # ---------------- corrected misfit (optional) ---------------- #
            # misfit_fx, misfit_fy, (nx, ny), (mx, my) = compute_misfit_optimal(
            #     match, substrate, film, max_supercell=10
            # )

            record = {
                "hkl_sub": sub_hkl,
                "hkl_film": film_hkl,
                # "misfit_f": (misfit_fx, misfit_fy),
                # "supercell_sub": (nx, ny),
                # "supercell_film": (mx, my),
                "von_mises": round(match.von_mises_strain, 4),
                "match_obj": match,
            }

            von_mises = match.von_mises_strain
            results_all.append(record)

            if von_mises <= config.von_mises_limit:
                results_filtered.append(record)
            else:
                continue

            # Verbose logging (optional)
            # if logger:
            #     logger.info(
            #         f"hkl(sub/film) = {sub_hkl}/{film_hkl} | "
            #         f"von Mises = {von_mises*100:.1f}% "
            #     )

    if logger:
        logger.info("=== Compact match summary ===")
    log_matches_compact(results_filtered, logger=logger)

    # ============================================================
    # DEDUPLICATION (very important!)
    # ============================================================

    before_all = len(results_all)
    before_filtered = len(results_filtered)

    results_all = _deduplicate_matches(results_all)
    results_filtered = _deduplicate_matches(results_filtered)

    if logger:
        logger.info(
            f"Deduplicated matches: "
            f"all {before_all} → {len(results_all)}, "
            f"filtered {before_filtered} → {len(results_filtered)}"
        )

    return {
        "all": results_all,
        "filtered": results_filtered,
    }