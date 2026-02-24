import logging
from pathlib import Path
from typing import List, Dict, Any, Optional

from pymatgen.core import Structure

from .config import InterfaceConfig
from .matching import find_matches
from .builder import build_interfaces





def get_logger(name="InterfaceConstructor", level=logging.INFO):
    logger = logging.getLogger(name)
    if not logger.handlers:
        ch = logging.StreamHandler()
        formatter = logging.Formatter(
            "%(asctime)s [%(levelname)s] %(message)s", datefmt="%H:%M:%S"
        )
        ch.setFormatter(formatter)
        logger.addHandler(ch)
        logger.setLevel(level)
    return logger


class InterfaceConstructor:
    """
    High-level API for interface construction.
    """

    def __init__(
        self,
        substrate: Structure,
        film: Structure,
        config: Optional[InterfaceConfig] = None,
    ):
        self.substrate = substrate
        self.film = film
        self.config = config or InterfaceConfig()
        self.logger = get_logger()
        self.matches: Optional[Dict[str, List[Dict[str, Any]]]] = None
        self.interfaces_metadata: List[Dict[str, Any]] = []
        self.interfaces_structures: List[Any] = []

    def run_matching(self):
        """Run orientation matching and filter by von Mises strain."""
        self.logger.info("Running orientation matching...")
        self.matches = find_matches(
            substrate=self.substrate,
            film=self.film,
            config=self.config,
            logger=self.logger,
        )
        self.logger.info(
            f"Total matches: {len(self.matches['all'])}, "
            f"Filtered matches: {len(self.matches['filtered'])}"
        )
        

        return self.matches

    def build_all_interfaces(self, output_dir: str = "interfaces"):
        """Build all interfaces for all filtered matches."""
        if not self.matches:
            raise RuntimeError(
                "No matches found. Run `run_matching()` first."
            )

        self.logger.info("Building all filtered interfaces...")
        for match in self.matches["filtered"]:
            self.logger.info(
                f"Processing match: {match['hkl_sub']}/{match['hkl_film']}, "
                f"von Mises = {match['von_mises']*100:.2f}%"
            )

            results, interfaces = build_interfaces(
                substrate_bulk=self.substrate,
                film_bulk=self.film,
                match_record=match,
                config=self.config,
                output_dir=output_dir,
            )

            self.interfaces_metadata.extend(results)
            self.interfaces_structures.extend(interfaces)

        self.logger.info(
            f"Total accepted interfaces: {len(self.interfaces_metadata)}"
        )
        return self.interfaces_metadata, self.interfaces_structures

    def save_summary_csv(self, output_dir: str = "interfaces"):
        """Save a CSV summary if not already saved by build_interfaces."""
        import csv

        if not self.interfaces_metadata:
            self.logger.warning("No interfaces metadata to save.")
            return

        csv_path = Path(output_dir) / f"{self.substrate.composition.reduced_formula}_{self.film.composition.reduced_formula}_summary.csv"
        with open(csv_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=list(self.interfaces_metadata[0].keys()))
            writer.writeheader()
            writer.writerows(self.interfaces_metadata)

        self.logger.info(f"CSV summary saved: {csv_path}")


    