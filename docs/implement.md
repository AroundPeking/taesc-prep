# Implementation details

## ABACUS

Output filenames follow LibRPA's `input_preset=abacus`: `stru_out.txt` and
`bz_sample.txt`.

### Structure

The program writes `/path/to/calculation/stru_out.txt`. It selects geometry in this
order:

1. `OUT.<suffix>/STRU_FINAL`.
2. The input structure specified by `stru_file`, or `STRU` by default.
3. A single-geometry running log (such as `running_scf.log`) in `OUT.<suffix>`,
   then in the supplied directory.

The `suffix`, `calculation`, and `stru_file` settings come from `INPUT` when
available. The suffix defaults to `ABACUS`. Structure files take priority over
logs to retain their numerical precision. Logs are used only when neither
structure file exists; ambiguous logs cause an error.

The structure reader supports explicit `LATTICE_VECTORS` and the `Direct`,
`Cartesian`, `Cartesian_au`, `Cartesian_angstrom`, and
`Cartesian_angstrom_center_*` coordinate modes. It preserves species and atom
order. For both `STRU_FINAL` and the input structure, lattice vectors in Bohr
are `LATTICE_VECTORS` multiplied by `LATTICE_CONSTANT`. Before converting
fractional positions to Cartesian coordinates, components outside `[0, 1)`
are translated by integer lattice periods into that interval (for example,
`-0.25` becomes `0.75`, and `1.0` becomes `0.0`). Components already within
the interval are unchanged.

Only Angstrom coordinate modes use `BOHR_TO_A = 0.5291770` to convert atomic
coordinates to Bohr. For centered coordinates, the reader selects the requested
Cartesian components of the cell center, as ABACUS does.

For relaxation or dynamics, the selected structure must correspond to the
exported electronic data. Structures specified through `latname` need an
explicit-lattice snapshot or the log route.

The program includes spatial symmetry operations when available in the
ABACUS output. These operations use the
fractional `row` convention. Rotations and translations retain the precision
of the ABACUS export. A symmetry count in `running_*.log` is insufficient
to reconstruct the operations. All input files must belong to the same geometry.
No k-point block is written.

`symrot_R.txt` exports the cell's full point group, while the retired
`stru_out.txt` kept only the operations that map the sampled k-mesh onto
itself. The program applies the same filter, using the mesh it reads from
`KPT.info`; operations that a non-cubic mesh breaks are dropped, because
LibRPA reconstructs k-stars from this set. When `KPT.info` cannot be read the
filter is skipped and every exported operation is written.

For `nspin 4` runs the program also writes the `spin_symmetry` block, which
carries the SU(2) representative of every operation. The convention is
`U sigma_i U^dagger = sum_j W_ij sigma_j` with `W = R_proper^T`, matching
ABACUS's `ModuleSymmetry::SpinRotation::so3_to_su2`. Magnetic runs add the
antiunitary elements `Theta*g` from the `Antiunitary operations:` trailer
that ABACUS appends to `symrot_R.txt`; when that trailer is absent the block
holds the unitary operations only.

An existing `stru_out.txt` is replaced after parsing succeeds. The program uses
only the Python standard library. The current reader supports ABACUS;
FHI-aims directory support is reserved for a later extension.

The ABACUS conventions follow its
[coordinate conversion source](https://github.com/deepmodeling/abacus-develop/blob/develop/source/source_cell/read_atoms_helper.cpp)
and [symmetry export source](https://github.com/deepmodeling/abacus-develop/blob/develop/source/source_lcao/module_ri/module_exx_symmetry/symm_rot_out.cpp).

### Brillouin-zone sampling

ABACUS conversion requires `OUT.<suffix>/KPT.info` and writes `bz_sample.txt`.
The output suffix follows `INPUT`, with `ABACUS` as the default.

The parser reads reduced fractional points and the full-grid reduction table.
It infers mesh divisions from the distinct fractional coordinates of the full
grid, then checks completeness and uniform spacing. Reduced coordinates retain
their printed values. Full-grid coordinates are placed in `[0, 1)` before
Cartesian conversion. Cartesian vectors are fractional row vectors multiplied
by the structure's reciprocal lattice, including `2*pi`, in inverse Bohr.
This uses the same selected structure as `stru_out.txt`.

The format follows ABACUS commit `52bc7978961c6fb4be330205cf71f00438c68497`:
[`rpa_lri_bz.cpp`](https://github.com/deepmodeling/abacus-develop/blob/52bc7978961c6fb4be330205cf71f00438c68497/source/source_lcao/module_rpa/rpa_lri_bz.cpp)
and [`rpa_lri_bz_io.cpp`](https://github.com/deepmodeling/abacus-develop/blob/52bc7978961c6fb4be330205cf71f00438c68497/source/source_lcao/module_rpa/rpa_lri_bz_io.cpp).

The output contains:

1. Three mesh divisions.
2. SCF and irreducible point counts, both equal to the reduced count.
3. Reduced-point rows: index, normalized weight, three fractional coordinates,
   three Cartesian coordinates, irreducible index, and SCF representative index.
4. If the grid is reduced, `full_kmap N` followed by the full-grid rows:
   index, three Cartesian coordinates, and SCF representative index.

All indices are one-based. Reduced and full-grid row orders follow `KPT.info`.
Weights equal the number of full-grid points assigned to each representative,
divided by the full-grid size. Printed weights are checked against these
multiplicities, allowing their decimal rounding error. This avoids carrying
four-decimal weight rounding into the output. Coordinates remain limited by
`KPT.info` precision. If no reduction table is present, the listed points must
form a complete uniform grid and receive identity mappings. Band paths and
incomplete grids raise errors.

## Output conventions

- Lattice vectors and Cartesian positions are in Bohr.
- Reciprocal lattice vectors are in inverse Bohr and include `2*pi`:
  `A @ B.T = 2*pi*I`, with lattice vectors stored as rows.
- Atom types start at 1 and follow the ABACUS species order.
- The optional spatial symmetry block uses fractional row-vector operations:
  `x_new = x @ R + t`.

See the [LibRPA dataset format](https://minyez.github.io/LibRPA/develop/dataset_format.html#stru-out)
for the file specification.

## Python interface

`detect_source(directory)` identifies the calculation and returns a `Source`
instance. The CLI uses this function before conversion. ABACUS detection
uses the boolean result of `discover_abacus`, which checks file markers only.
`AbacusSource` selects the structure file and parser during extraction, using
`stru_file`, `suffix`, and `calculation` from `INPUT`. It locates `KPT.info`
when BZ sampling is extracted. Initialization only stores the directory and
enables the supported conversions; it does not read configuration or cache
input file paths. Unrecognized directories produce an error.

`Source` implements the common conversion workflow. `AbacusSource` supplies
ABACUS data through `_extract_stru`:

```python
from taesc_prep import Source, detect_source

source: Source = detect_source("/path/to/calculation")
result = source.convert()
print(result.output)
```

`convert(output=None, *, include_symmetry=True)` runs the enabled conversions
and returns a `ConversionResult` with `output`, `n_atoms`, `n_symops`, `n_kpoints`,
`n_full_kpoints`, and `notes`. K-point counts are zero when BZ conversion is disabled.
`output` is a destination directory and defaults to the input directory.
The returned `result.output` is that directory; the structure file is
`result.output / source.stru_filename`. Set `include_symmetry=False` to omit spatial
symmetry operations. If no conversion is enabled, the result has `output=None`,
zero counts, and no notes; no files are read or written.

The shared conversion methods use `stru_filename` and `bz_sampling_filename`
from the concrete source.

`Source` initializes `process_stru=False` and `process_bz_sampling=False`. A concrete source enables structure
processing in its initializer and implements `_extract_stru`:

```python
class MySource(Source):
    def __init__(self, directory):
        super().__init__(directory)
        self.process_stru = True

    def _extract_stru(self, *, include_symmetry=True) -> dict:
        # Read the calculation files and return the dictionary described below.
        ...
```

Concrete sources inherit `convert` and `convert_stru` without overriding them.
When `process_stru` is True, `convert` calls the shared `convert_stru` method.
That method calls `_extract_stru`, passes its dictionary to `write_taesc_stru`
in `output.py`, and returns the extracted dictionary. Only `convert` constructs
and returns `ConversionResult`.

`_extract_stru` returns a dictionary with the following keys:

| Key | Required | Value |
| --- | --- | --- |
| `lattice` | Yes | Three row vectors, each with three components, in Bohr. |
| `positions` | Yes | One Cartesian coordinate row per atom, in Bohr. |
| `types` | Yes | One integer type per atom, numbered from 1. |
| `operations` | No | `(rotation, translation)` pairs in the fractional row convention. Each rotation is a flat list of nine integers in row-major order; each translation has three components. |
| `notes` | No | Conversion notes as a sequence of strings. |
| `source_paths` | No | Input paths as strings or `Path` objects, protected from overwriting. |

Missing optional keys default to empty sequences. The writer computes reciprocal
vectors from `lattice`. It requires the three geometry keys and uses `operations`
when present. The shared conversion method handles `notes` and `source_paths`.

For BZ conversion, a concrete source sets `process_bz_sampling=True` and
implements `_extract_bz_sampling`. `Source.convert_bz_sampling` adds the lattice
from the extracted structure and calls `write_taesc_bz_sampling` in `output.py`.
It returns the sampling dictionary. The ABACUS constructor always sets
`process_bz_sampling=True`.

The sampling dictionary has these required keys:

| Key | Value |
| --- | --- |
| `mesh` | Three positive mesh divisions. |
| `direct` | Reduced fractional k vectors, in SCF order. |
| `weights` | Positive weights normalized to sum to one. |
| `full_direct` | Full-grid fractional vectors, in mapping-table order. |
| `full_map` | One-based indices into `direct`, one per full-grid point. |

`notes` and `source_paths` are optional, as for structure extraction.
The shared conversion method supplies `lattice` before writing. The writer
normalizes weights and validates dimensions, representative indices, and
multiplicities before replacing the destination file.
