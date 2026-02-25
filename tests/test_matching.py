import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[1]))

from pymatgen.core import Structure
from interface_constructor import InterfaceConstructor, InterfaceConfig

# -------------------------------
# Загружаем структуры
# -------------------------------
# substrate_file = "Li3GaN2.cif"
substrate_file = "Li8ZrO6.cif"
substrate = Structure.from_file(substrate_file)

# substrate.add_oxidation_state_by_element({"Li": +1, "N": -3, "Ga": +3})
substrate.add_oxidation_state_by_element({"Li": +1, "O": -2, "Zr": +4})

film_file = "Li.cif"


film = Structure.from_file(film_file)
film.add_oxidation_state_by_element({"Li": +1})

# -------------------------------
# Настройки интерфейсов
# -------------------------------
config = InterfaceConfig(
    # sub_miller=(1, 1, 0),  # фиксированный Miller для подложки
    # film_miller=(1, 1, 0),  # фиксированный Miller для пленки
    max_sub_miller=1,  # фиксированный Miller для подложки
    max_film_miller=1,  # фиксированный Miller для пленки
    von_mises_limit=0.05,  # максимум von Mises
    num_sites_limit=200,  # отсечка по числу атомов
    density_limit=0.01,  # мин плотность
    charge_limit=0.2,    # макс заряд
    film_thickness=4,
    substrate_thickness=4,
    surface_thickness=0.8,
    gap=2.0,
    vacuum_over_film=15.0,
    output_folder="interfaces_test"
)

# -------------------------------
# Инициализация конструктора
# -------------------------------
ic = InterfaceConstructor(substrate, film, config)

# -------------------------------
# 1) Поиск совпадений
# -------------------------------
matches = ic.run_matching()
print(f"Найдено {len(matches['filtered'])} подходящих совпадений.")

# -------------------------------
# 2) Построение интерфейсов
# -------------------------------
metadata, structures = ic.build_all_interfaces(output_dir=config.output_folder)
print(f"Построено {len(metadata)} интерфейсов.")

# -------------------------------
# 3) Сохранение CSV summary
# -------------------------------
ic.save_summary_csv(output_dir=config.output_folder)
print("CSV summary сохранен.")
