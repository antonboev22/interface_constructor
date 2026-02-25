import logging
from pathlib import Path
from typing import List, Dict, Any, Optional

from pymatgen.core import Structure

from .config import InterfaceConfig
from .matching import find_matches
from .builder import build_interfaces, save_csv


def get_logger(name="InterfaceConstructor", level=logging.DEBUG):
    logger = logging.getLogger(name)
    if not logger.hasHandlers():
        handler = logging.StreamHandler()
        handler.setFormatter(
            logging.Formatter("%(asctime)s [%(levelname)s] %(message)s", datefmt="%H:%M:%S")
        )
        logger.addHandler(handler)
        logger.setLevel(level)
    return logger


class InterfaceConstructor:
    """High-level API for automated interface construction."""

    @staticmethod
    def _get_logger(name="InterfaceConstructor", level=logging.INFO):
        logger = logging.getLogger(name)
        if not logger.hasHandlers():
            handler = logging.StreamHandler()
            handler.setFormatter(
                logging.Formatter("%(asctime)s [%(levelname)s] %(message)s", datefmt="%H:%M:%S")
            )
            logger.addHandler(handler)
            logger.setLevel(level)
        return logger

    def __init__(
        self,
        substrate: Structure,
        film: Structure,
        config: Optional[InterfaceConfig] = None,
    ):
        self.substrate = substrate
        self.film = film
        self.config = config or InterfaceConfig()
        self.logger = self._get_logger()
        self.matches: Optional[Dict[str, List[Dict[str, Any]]]] = None
        self.interfaces_metadata: List[Dict[str, Any]] = []
        self.interfaces_structures: List[Any] = []

    def run_matching(self):
        """Run orientation matching and filter by von Mises strain."""
        self.logger.info("🔍 Starting orientation matching...")
        self.matches = find_matches(
            substrate=self.substrate,
            film=self.film,
            config=self.config,
            logger=self.logger,
        )
        self.logger.info(
            f"✅ Matching completed: {len(self.matches['all'])} total matches, "
            f"{len(self.matches['filtered'])} passed the von Mises filter"
        )
        return self.matches

    def build_all_interfaces(self, output_dir: str = "interfaces"):
        """Build interfaces for all filtered matches with detailed logging."""
        if not self.matches or not self.matches["filtered"]:
            raise RuntimeError("No matches found. Run `run_matching()` first.")

        self.logger.info(
            f"Building filtered interfaces using the following criteria:\n"
            f"  - Surface density limit: {self.config.density_limit}\n"
            f"  - Max total charge: {self.config.charge_limit}\n"
            f"  - Max num_sites: {self.config.num_sites_limit}\n"
            f"  - Surface thickness: {self.config.surface_thickness}\n"
            f"  - Film thickness: {self.config.film_thickness}, "
            f"Substrate thickness: {self.config.substrate_thickness}\n"
            f"  - Gap: {self.config.gap}, Vacuum above film: {self.config.vacuum_over_film}\n\n"
        )

        for idx, match in enumerate(self.matches["filtered"], start=1):
            self.logger.info(
                f"=== Match {idx}/{len(self.matches['filtered'])}: "
                f"{match['hkl_sub']}/{match['hkl_film']} | "
                f"von Mises={match['von_mises']*100:.2f}% ==="
            )

            results, interfaces = build_interfaces(
                substrate_bulk=self.substrate,
                film_bulk=self.film,
                match_record=match,
                config=self.config,
                output_dir=output_dir,
                logger=self.logger,
            )

            # Per-interface detailed logging
            for j, r in enumerate(results, start=1):
                status = "✅ accepted" if r["passed_filters"] else "❌ rejected"
                self.logger.info(
                    f"\tInterface {j}: {r['slab']}\n"
                    f"\t\t\t\tnum_sites = {r['num_sites']}, "
                    f"density = {r['substrate_density']:.2f}/{r['film_density']:.2f} (s/f), "
                    f"total_charge = {r['abs_charge_density']:.2f} → {status}\n"
                )

            # Match summary
            accepted_count = sum(r["passed_filters"] for r in results)
            self.logger.info(
                f"\n\n  → Match {idx} summary: "
                f"{accepted_count}/{len(results)} interfaces accepted\n"
                f"--------------\n"
            )

            self.interfaces_metadata.extend(results)
            self.interfaces_structures.extend(interfaces)

        self.logger.info(
            f"🏁 Total accepted interfaces: {len(self.interfaces_structures)}"
        )
        return self.interfaces_metadata, self.interfaces_structures

    def save_summary_csv(self, output_dir: str = "interfaces"):
        """Save CSV summary of all generated interfaces."""
        save_csv(self.interfaces_metadata, self.substrate, self.film, output_dir)