"""Read ABACUS calculation directory."""

from collections import Counter
from decimal import Decimal
import math
from pathlib import Path
import re

from .source import Source
from .spin import cartesian_rotation, flatten_su2, rotation_to_su2
from .structure import fractional_coordinates, multiply_row, reciprocal_vectors

__all__ = ["AbacusSource", "discover_abacus"]

# ABACUS source_base/constants.h: BOHR_TO_A.
# Only used when reading from OUT.ABACUS/running_scf.log
BOHR_TO_ANGSTROM = 0.5291770
SECTIONS = {
    "ATOMIC_SPECIES", "NUMERICAL_ORBITAL", "NUMERICAL_DESCRIPTOR",
    "PAW_FILES", "LATTICE_CONSTANT", "LATTICE_VECTORS", "LATTICE_PARAMETERS",
    "ATOMIC_POSITIONS",
}


class AbacusSource(Source):
    """Convert the automatically discovered files in an ABACUS directory."""

    process_stru = True
    stru_filename = "stru_out.txt"

    process_bz_sampling = True
    bz_sampling_filename = "bz_sample.txt"

    def __init__(self, directory):
        super().__init__(directory)

    def _extract_stru(self, *, include_symmetry=True) -> dict:
        """Read ABACUS geometry and optional spatial symmetry operations."""
        source, reader, outdir = self._select_structure()
        structure = reader(source.read_text())
        symfile = None
        if include_symmetry:
            candidates = list(dict.fromkeys([self.directory / "symrot_R.txt",
                                            outdir / "symrot_R.txt"]))
            found = [p for p in candidates if p.is_file()]
            if len(found) > 1:
                raise ValueError("Multiple spatial symmetry exports found; use a directory containing one calculation")
            symfile = found[0] if found else None
        operations, antiunitary = ([], [])
        if symfile:
            operations, antiunitary = read_symrot(symfile.read_text())
            # symrot_R.txt exports the cell's full point group, while
            # stru_out.txt kept only the operations that map the sampled
            # k-mesh onto itself. Apply that filter so LibRPA's k-star
            # reconstruction sees a usable set.
            mesh = self._sampled_mesh()
            if mesh is not None:
                operations = [operation for operation in operations
                              if preserves_mesh(operation[0], mesh)]
                antiunitary = [operation for operation in antiunitary
                               if preserves_mesh(operation[0], mesh, (0.0, 0.0, 0.0),
                                                 antiunitary=True)]
            structure["operations"] = operations
        settings = read_input_settings(self.directory)
        structure["nspin"] = int(settings.get("nspin", 1))
        if structure["nspin"] == 4 and structure.get("operations"):
            # SOC needs the spinor representation of every operation. ABACUS
            # writes that block into stru_out.txt; rebuild it here from the
            # spatial operations so stru_out.txt is no longer required.
            def spin_rows(operations, antiunitary_flag):
                return [(antiunitary_flag,
                         flatten_su2(rotation_to_su2(
                             cartesian_rotation(structure["lattice"],
                                                [rotation[3*i:3*i+3] for i in range(3)]))))
                        for rotation, _ in operations]

            structure["spin_symmetry"] = {
                "grey_group": 0 if antiunitary else 1,
                "spin_source": 1,
                "operations": (spin_rows(operations, 0)
                               + spin_rows(antiunitary, 1)),
            }
        structure["source_paths"] = (source, symfile) if symfile else (source,)
        if reader is read_log:
            structure["notes"] = ("geometry precision is limited by the running log; prefer STRU_FINAL when available.",)
        return structure

    def _sampled_mesh(self):
        """Mesh divisions from KPT.info, or None when it cannot be read.

        The structure conversion runs first, so a missing or unreadable
        KPT.info must not block it; the symmetry filter is simply skipped.
        """
        try:
            settings = read_input_settings(self.directory)
            outdir = self.directory / ("OUT." + settings.get("suffix", "ABACUS"))
            return read_kpt_info((outdir / "KPT.info").read_text())["mesh"]
        except (OSError, ValueError, KeyError):
            return None

    def _extract_bz_sampling(self) -> dict:
        settings = read_input_settings(self.directory)
        outdir = self.directory / ("OUT." + settings.get("suffix", "ABACUS"))
        kpt_info = outdir / "KPT.info"
        sampling = read_kpt_info(kpt_info.read_text())
        sampling["source_paths"] = (kpt_info,)
        return sampling

    def _select_structure(self):
        """Select geometry and its parser when structure extraction is requested."""
        directory = self.directory
        if not directory.is_dir():
            raise ValueError(f"Expected a calculation directory: {directory}")
        settings = read_input_settings(directory)
        stru = directory / settings.get("stru_file", "STRU")
        outdir = directory / ("OUT." + settings.get("suffix", "ABACUS"))
        final_stru = outdir / "STRU_FINAL"
        if final_stru.is_file():
            return final_stru, read_stru, outdir
        if stru.is_file():
            return stru, read_stru, outdir
        calculation = settings.get("calculation")
        logs = ([outdir / f"running_{calculation}.log"] if calculation
                else sorted(outdir.glob("running_*.log")))
        logs = [path for path in logs if path.is_file()]
        if not logs:
            logs = sorted(directory.glob("running_*.log"))
        if len(logs) == 1:
            return logs[0], read_log, logs[0].parent
        if len(logs) > 1:
            raise ValueError(
                "Multiple running logs found; specify calculation in INPUT or use a directory containing one log")
        raise ValueError("No ABACUS structure or running log found; check stru_file and suffix in INPUT")


def _is_magnetic_nspin4(directory):
    """True when the run carries magnetization, i.e. antiunitary operations exist.

    ABACUS builds the Shubnikov group only for nspin=4 with a non-zero moment.
    In STRU the moment is the number between each species label and its atom
    count inside ATOMIC_POSITIONS; INPUT's nupdown overrides it.
    """
    directory = Path(directory)
    settings = read_input_settings(directory)
    if "nupdown" in settings:
        try:
            if abs(float(settings["nupdown"])) > 0.0:
                return True
        except ValueError:
            pass
    stru = directory / settings.get("stru_file", "STRU")
    if not stru.is_file():
        return False

    sections = {}
    current = None
    for raw in stru.read_text().splitlines():
        line = re.split(r"#|//", raw, maxsplit=1)[0].strip()
        if not line:
            continue
        if line in SECTIONS:
            current = sections[line] = []
        elif current is not None:
            current.append(line)

    species = [row.split()[0] for row in sections.get("ATOMIC_SPECIES", [])]
    rows = iter(sections.get("ATOMIC_POSITIONS", []))
    next(rows, None)  # Coordinate mode.
    magnetizations = []
    for label in species:
        if next(rows, None) != label:
            return False  # Let read_stru report the malformed block.
        try:
            magnetizations.append(float(next(rows).replace("D", "E").replace("d", "e")))
            int(next(rows))
        except (StopIteration, ValueError, AttributeError):
            return False
    return any(abs(x) > 0.0 for x in magnetizations)


def numbers(line, count=3):
    tokens = line.split()
    if len(tokens) < count:
        raise ValueError(f"Expected {count} numbers: {line!r}")
    values = [float(x.replace("D", "E").replace("d", "e"))
              for x in tokens[:count]]
    if not all(math.isfinite(x) for x in values):
        raise ValueError(f"Non-finite number: {line!r}")
    return values


def read_stru(text):
    """Read explicit lattice vectors and preserve ABACUS species/atom order."""
    sections = {}
    current = None
    for raw in text.splitlines():
        line = re.split(r"#|//", raw, maxsplit=1)[0].strip()
        if not line:
            continue
        if line in SECTIONS:
            if line in sections:
                raise ValueError(f"Duplicate STRU section: {line}")
            current = sections[line] = []
        elif current is not None:
            current.append(line)
        else:
            raise ValueError(f"Expected a STRU section: {line!r}")

    for key in ("ATOMIC_SPECIES", "LATTICE_CONSTANT", "LATTICE_VECTORS", "ATOMIC_POSITIONS"):
        if not sections.get(key):
            raise ValueError(f"Missing {key}; use an explicit-lattice STRU or a running log")
    if len(sections["LATTICE_CONSTANT"]) != 1 or len(sections["LATTICE_VECTORS"]) != 3:
        raise ValueError("Expected one lattice constant and three lattice vectors")
    scale = numbers(sections["LATTICE_CONSTANT"][0], 1)[0]
    if scale <= 0:
        raise ValueError("LATTICE_CONSTANT must be positive")
    lattice = [[scale*x for x in numbers(row)] for row in sections["LATTICE_VECTORS"]]
    reciprocal = reciprocal_vectors(lattice)
    species = [row.split()[0] for row in sections["ATOMIC_SPECIES"]]
    if len(set(species)) != len(species):
        raise ValueError("Duplicate ATOMIC_SPECIES label")
    rows = iter(sections["ATOMIC_POSITIONS"])
    mode = next(rows).lower()
    factors = {"cartesian": scale, "cartesian_au": 1.0,
               "cartesian_angstrom": 1 / BOHR_TO_ANGSTROM}
    shift = [0.0]*3
    if mode.startswith("cartesian_angstrom_center_"):
        center = mode[len("cartesian_angstrom_center_"):]
        if center not in ("xy", "xz", "yz", "xyz"):
            raise ValueError(f"Unsupported coordinate mode: {mode}")
        # ABACUS selects Cartesian components of the full cell center.
        shift = [sum(row[i] for row in lattice)/2 if axis in center else 0.0
                 for i, axis in enumerate("xyz")]
        mode = "cartesian_angstrom"
    if mode != "direct" and mode not in factors:
        raise ValueError(f"Unsupported coordinate mode: {mode}")
    positions, types = [], []
    try:
        for it, label in enumerate(species, 1):
            if next(rows) != label:
                raise ValueError(f"Expected ATOMIC_POSITIONS species {label}")
            numbers(next(rows), 1)  # Species magnetization.
            count = int(next(rows))
            if count <= 0:
                raise ValueError(f"Atom count for {label} must be positive")
            for _ in range(count):
                xyz = numbers(next(rows))
                if mode == "direct":
                    frac = xyz
                else:
                    xyz = [x*factors[mode] + s for x, s in zip(xyz, shift)]
                    frac = fractional_coordinates(xyz, reciprocal)
                # Translate fractional coordinates into [0, 1) before Cartesian conversion.
                frac = [x if 0.0 <= x < 1.0 else x % 1.0 for x in frac]
                xyz = multiply_row(frac, lattice)
                positions.append(xyz)
                types.append(it)
    except StopIteration as exc:
        raise ValueError("Incomplete ATOMIC_POSITIONS section") from exc
    if next(rows, None) is not None:
        raise ValueError("Unexpected data after ATOMIC_POSITIONS")
    return {"lattice": lattice, "positions": positions, "types": types}


def read_log(text):
    """Read the setup geometry from a log with one atomic coordinate table.

    ABACUS calls its lattice scaling factor a_0 in the lattice-table header.
    Thus these vectors, and Cartesian atom coordinates, need that factor.
    """
    lines = text.splitlines()

    def unique_match(pattern):
        matches = list(re.finditer(pattern, text, re.MULTILINE | re.IGNORECASE))
        if len(matches) != 1:
            raise ValueError(f"Expected one log entry matching {pattern!r}; found {len(matches)}")
        return matches[0]

    # A relaxation log can contain later positions without another setup table.
    ionic_steps = [int(x) for x in re.findall(r"\bION\s*=\s*(\d+)", text)]
    if any(x > 1 for x in ionic_steps) or re.search(r"\b(?:STEP OF RELAXATION|MD STEP)\b", text):
        raise ValueError("Use a snapshot STRU for relaxation or molecular-dynamics output")
    scale = numbers(unique_match(r"^\s*lattice constant \(Bohr\)\s*=\s*(\S+)").group(1), 1)[0]
    if scale <= 0:
        raise ValueError("Lattice constant must be positive")
    nat = int(unique_match(r"^\s*TOTAL ATOM NUMBER\s*=\s*(\d+)").group(1))
    blocks = re.split(r"(?im)^\s*READING ATOM TYPE\s+(\d+)\s*$", text)
    types, labels = [], []
    for index in range(1, len(blocks), 2):
        it = int(blocks[index])
        body = blocks[index + 1]
        label = re.search(r"(?i)atom label\s*=\s*(\S+)", body)
        count = re.search(r"(?i)number of atoms? for this type\s*=\s*(\d+)", body)
        if not label or not count or it != len(labels) + 1:
            raise ValueError("Incomplete or out-of-order atom-type information in log")
        labels.append(label.group(1))
        types.extend([it] * int(count.group(1)))
    if nat <= 0 or len(types) != nat:
        raise ValueError("Atom-type counts disagree with TOTAL ATOM NUMBER")
    lattice_headers = [i for i, line in enumerate(lines)
                       if line.strip() == "Lattice vectors: (Cartesian coordinate: in unit of a_0)"]
    coord_headers = [i for i, line in enumerate(lines)
                     if re.match(r"^\s*(DIRECT|CARTESIAN) COORDINATES(?:\s|$)", line)]
    if len(lattice_headers) != 1 or len(coord_headers) != 1:
        raise ValueError("Expected one lattice table and one atomic coordinate table; use a snapshot STRU")
    ilat = lattice_headers[0]
    if len(lines[ilat+1:ilat+4]) != 3:
        raise ValueError("Incomplete lattice table")
    lattice = [[scale*x for x in numbers(row)] for row in lines[ilat+1:ilat+4]]
    ic = coord_headers[0]
    direct = lines[ic].strip().startswith("DIRECT")
    rows = lines[ic+2:ic+2+nat]  # Skip column headings.
    if len(rows) != nat:
        raise ValueError("Incomplete atomic coordinate table")
    positions = []
    for row, it in zip(rows, types):
        fields = row.split()
        label = labels[it - 1]
        if not fields or (fields[0] != label and
                          not re.fullmatch(r"tau[cd]_" + re.escape(label) + r"\d+", fields[0])):
            raise ValueError(f"Expected coordinates for {label}: {row!r}")
        xyz = numbers(" ".join(fields[1:]))
        positions.append(multiply_row(xyz, lattice) if direct else [scale*x for x in xyz])
    return {"lattice": lattice, "positions": positions, "types": types}


def reciprocal_rotation(direct_rotation):
    """Return the reciprocal-space matrix of a fractional row-vector rotation.

    ABACUS records ``k' = k * kg`` for row-vector k, so ``kg`` is the inverse
    transpose of the direct-space operation.
    """
    a, b, c, d, e, f, g, h, i = direct_rotation
    det = a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g)
    if abs(det) != 1:
        raise ValueError("Symmetry rotation must be unimodular")
    return [det * (e * i - f * h), det * (f * g - d * i), det * (d * h - e * g),
            det * (c * h - b * i), det * (a * i - c * g), det * (b * g - a * h),
            det * (b * f - c * e), det * (c * d - a * f), det * (a * e - b * d)]


def preserves_mesh(direct_rotation, mesh, offset=(0.0, 0.0, 0.0), antiunitary=False):
    """Whether an operation maps a Monkhorst-Pack mesh onto itself.

    ABACUS keeps only such operations in ``stru_out.txt`` because the others
    cannot be used to reconstruct k-stars. ``antiunitary`` applies the
    time-reversal sign that Theta*g carries.
    """
    reciprocal = reciprocal_rotation(direct_rotation)
    sign = -1 if antiunitary else 1
    for source in range(3):
        for target in range(3):
            if sign * reciprocal[3 * source + target] * mesh[target] % mesh[source] != 0:
                return False
    origin = [offset[index] / mesh[index] for index in range(3)]
    for target in range(3):
        mapped = sign * sum(origin[source] * reciprocal[3 * source + target]
                            for source in range(3))
        shift = (mapped - origin[target]) * mesh[target]
        if abs(shift - round(shift)) > 1e-8:
            return False
    return True


def read_symrot(text):
    """Extract fractional row-convention operations from the ABACUS spatial symmetry export."""
    lines = [line.strip() for line in text.splitlines() if line.strip()]
    lmax = []
    for label in ("AOs", "ABFs"):
        matches = re.findall(r"Lmax of " + label + r":\s*(\d+)", text)
        if len(matches) != 1:
            raise ValueError(f"Expected one Lmax of {label} in the spatial symmetry export")
        lmax.append(int(matches[0]))
    start = next((i for i, line in enumerate(lines) if line == "0"), None)
    if start is None:
        raise ValueError("No symmetry operations in the spatial symmetry export")
    # Each angular momentum contributes 2*l+1 orbital-matrix rows. A
    # spatial-only export (nspin=4, where the collinear orbital rotations do
    # not exist) omits them and says so in the header.
    orbital_rows = 0 if "Spatial-only export" in text else (max(lmax) + 1)**2
    operations = []
    pos = start
    while pos < len(lines):
        # The optional antiunitary trailer terminates the unitary block.
        if lines[pos].startswith("Antiunitary operations:"):
            break
        if lines[pos] != str(len(operations)) or pos + 5 + orbital_rows > len(lines):
            raise ValueError("Incomplete or out-of-order symmetry operation")
        rotation = [x for row in lines[pos+1:pos+4] for x in numbers(row)]
        if any(abs(x-round(x)) > 1e-8 for x in rotation):
            raise ValueError("Noninteger fractional symmetry rotation")
        rotation = [int(round(x)) for x in rotation]
        translation = lines[pos+4]
        if not (translation.startswith("(") and translation.endswith(")")):
            raise ValueError("Expected a parenthesized symmetry translation")
        translation = numbers(translation[1:-1])
        operations.append((rotation, translation))
        pos += 5 + orbital_rows
    identity = [1, 0, 0, 0, 1, 0, 0, 0, 1]
    if not any(r == identity and all(abs(x-round(x)) < 1e-8 for x in t)
               for r, t in operations):
        raise ValueError("Symmetry operations do not contain the identity")
    return operations, read_symrot_antiunitary(lines[pos:])


def read_symrot_antiunitary(trailer):
    """Read the optional antiunitary trailer appended after the unitary block.

    ABACUS writes ``Antiunitary operations: N`` followed by N entries of the
    same shape as the unitary ones, but without orbital rotation matrices.
    Older exports have no trailer and yield an empty list.
    """
    start = next((i for i, line in enumerate(trailer)
                  if line.startswith("Antiunitary operations:")), None)
    if start is None:
        return []
    count = int(trailer[start].split(":", 1)[1])
    operations = []
    pos = start + 1
    while len(operations) < count:
        if pos + 4 > len(trailer) or trailer[pos] != str(len(operations)):
            raise ValueError("Incomplete or out-of-order antiunitary operation")
        rotation = [x for row in trailer[pos+1:pos+4] for x in numbers(row)]
        if any(abs(x-round(x)) > 1e-8 for x in rotation):
            raise ValueError("Noninteger fractional antiunitary rotation")
        translation = trailer[pos+4]
        if not (translation.startswith("(") and translation.endswith(")")):
            raise ValueError("Expected a parenthesized antiunitary translation")
        operations.append(([int(round(x)) for x in rotation],
                           numbers(translation[1:-1])))
        pos += 5
    return operations


def discover_abacus(directory) -> bool:
    """Identify an ABACUS directory from file markers without selecting inputs."""
    directory = Path(directory)
    if not directory.is_dir():
        return False
    if any((directory / name).is_file() for name in ("INPUT", "STRU")):
        return True
    return any(path.is_file()
               for pattern in ("running_*.log", "OUT.*/STRU_FINAL",
                               "OUT.*/running_*.log", "OUT.*/KPT.info")
               for path in directory.glob(pattern))


def read_input_settings(directory):
    """Read settings from an optional ABACUS INPUT file."""
    settings = {}
    input_file = directory / "INPUT"
    if input_file.is_file():
        for raw in input_file.read_text().splitlines():
            fields = re.split(r"#|//", raw, maxsplit=1)[0].split()
            if len(fields) >= 2:
                settings[fields[0].lower()] = fields[1]
    return settings


def read_kpt_info(text):
    """Read a uniform full grid and its reduction from ABACUS KPT.info.

    Return mesh, reduced fractional points, normalized weights, full fractional
    points, and one-based representative indices. Weights come from exact map
    multiplicities because KPT.info prints weights to limited precision.
    """
    lines = [line.strip() for line in text.splitlines() if line.strip()]
    counts = re.findall(r"\bnkstot now\s*=\s*(\d+)", text)
    if len(counts) != 1 or int(counts[0]) <= 0:
        raise ValueError("Expected one positive nkstot now in KPT.info")
    n_reduced = int(counts[0])

    def table(header, count, width):
        indices = [i for i, line in enumerate(lines) if line == header]
        if len(indices) != 1:
            raise ValueError(f"Expected one {header} table in KPT.info")
        start = indices[0] + 2  # Skip column headings.
        rows = [line.split() for line in lines[start:start + count]]
        if len(rows) != count or any(len(row) != width for row in rows):
            raise ValueError(f"Incomplete {header} table in KPT.info")
        if any(int(row[0]) != i for i, row in enumerate(rows, 1)):
            raise ValueError(f"Nonsequential indices in {header}")
        return rows

    rows = table("K-POINTS DIRECT COORDINATES", n_reduced, 5)
    reduced = [numbers(" ".join(row[1:4])) for row in rows]
    printed_weights = [numbers(row[4], 1)[0] for row in rows]
    # A zero printed weight can result from rounding on a large grid.
    if any(weight < 0 for weight in printed_weights) or sum(printed_weights) <= 0:
        raise ValueError("Invalid k-point weights in KPT.info")
    rounding = [0.5 * 10.0 ** Decimal(row[4].replace("D", "E").replace("d", "e")).as_tuple().exponent
                for row in rows]
    full_counts = re.findall(r"\bnkstot\s*=\s*(\d+)", text)
    if full_counts:
        if len(full_counts) != 1 or int(full_counts[0]) < n_reduced:
            raise ValueError("Invalid full-grid count in KPT.info")
        full_rows = table("K-POINTS REDUCTION ACCORDING TO SYMMETRY", int(full_counts[0]), 8)
        full = [numbers(" ".join(row[1:4])) for row in full_rows]
        mapping = [int(row[4]) for row in full_rows]
        for row, representative in zip(full_rows, mapping):
            if not 1 <= representative <= n_reduced:
                raise ValueError("KPT.info representative is outside the reduced grid")
            expected = reduced[representative - 1]
            reported = numbers(" ".join(row[5:8]))
            if any(abs(a - b - round(a - b)) > 2e-8 for a, b in zip(expected, reported)):
                raise ValueError("KPT.info representative coordinates disagree with the reduced table")
    else:
        if any("K-POINTS REDUCTION" in line for line in lines):
            raise ValueError("Missing full-grid count in KPT.info")
        full = reduced
        mapping = list(range(1, n_reduced + 1))

    # Infer divisions from the complete fractional grid, preserving its offsets.
    wrapped = [tuple(round(x % 1.0, 8) % 1.0 for x in point) for point in full]
    axes = [sorted({point[i] for point in wrapped}) for i in range(3)]
    mesh = [len(axis) for axis in axes]
    if math.prod(mesh) != len(full) or len(set(wrapped)) != len(full):
        raise ValueError("KPT.info does not contain a complete uniform three-dimensional grid")
    for axis in axes:
        spacing = 1.0 / len(axis)
        if any(abs((axis[(i + 1) % len(axis)] - value) % 1.0 - spacing) > 2e-8
               for i, value in enumerate(axis)) and len(axis) > 1:
            raise ValueError("KPT.info grid spacing is not uniform")
    grid = set(wrapped)
    if any(tuple(round(x % 1.0, 8) % 1.0 for x in point) not in grid for point in reduced):
        raise ValueError("KPT.info reduced point is absent from the full grid")
    multiplicity = Counter(mapping)
    if len(multiplicity) != n_reduced:
        raise ValueError("KPT.info full mapping omits a reduced point")
    weights = [multiplicity[i] / len(full) for i in range(1, n_reduced + 1)]
    total = sum(printed_weights)
    total_rounding = sum(rounding)
    for printed, expected, error in zip(printed_weights, weights, rounding):
        if abs(printed - expected * total) > error + expected * total_rounding + 1e-12:
            raise ValueError("KPT.info weight disagrees with full-grid multiplicity")
    return {"mesh": mesh, "direct": reduced, "weights": weights,
            "full_direct": wrapped, "full_map": mapping}
