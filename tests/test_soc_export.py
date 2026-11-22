"""Reading ABACUS's spatial-symmetry export, including the SOC trailer."""
import unittest

from taesc_prep.abacus import read_symrot, read_symrot_antiunitary

HEADER = """Number of irreducible sector: 1
Lmax of AOs: 1
Lmax of ABFs: 2
Format:
The index of the symmetry operation
"""

ORBITAL_ROWS = "\n".join(["1 0"] * 9)


def export(operations, antiunitary=()):
    """Assemble a symrot_R.txt body for the given operations."""
    lines = [HEADER]
    for index, rotation in enumerate(operations):
        lines.append(str(index))
        lines.extend(" ".join(str(x) for x in rotation[3 * i:3 * i + 3]) for i in range(3))
        lines.append("(0 0 0)")
        lines.append(ORBITAL_ROWS)
    if antiunitary is not None:
        lines.append("Antiunitary operations: {}".format(len(antiunitary)))
        for index, rotation in enumerate(antiunitary):
            lines.append(str(index))
            lines.extend(" ".join(str(x) for x in rotation[3 * i:3 * i + 3]) for i in range(3))
            lines.append("(0.5 0.5 0.5)")
    return "\n".join(lines) + "\n"


IDENTITY = [1, 0, 0, 0, 1, 0, 0, 0, 1]
CZ = [0, -1, 0, 1, 0, 0, 0, 0, 1]


class SymrotExportTests(unittest.TestCase):
    def test_only_unitary_operations(self):
        unitary, antiunitary = read_symrot(export([IDENTITY, CZ]))
        self.assertEqual(unitary, [(IDENTITY, [0.0, 0.0, 0.0]), (CZ, [0.0, 0.0, 0.0])])
        self.assertEqual(antiunitary, [])

    def test_magnetic_trailer_is_read(self):
        unitary, antiunitary = read_symrot(export([IDENTITY], [CZ]))
        self.assertEqual(len(unitary), 1)
        self.assertEqual(antiunitary, [(CZ, [0.5, 0.5, 0.5])])

    def test_zero_antiunitary_trailer(self):
        unitary, antiunitary = read_symrot(export([IDENTITY], []))
        self.assertEqual(len(unitary), 1)
        self.assertEqual(antiunitary, [])

    def test_legacy_export_without_trailer(self):
        text = export([IDENTITY], antiunitary=None)
        self.assertEqual(read_symrot_antiunitary(text), [])
        self.assertEqual(len(read_symrot(text)[0]), 1)


if __name__ == "__main__":
    unittest.main()


class SpatialOnlyExportTests(unittest.TestCase):
    """nspin=4 exports carry no orbital rotation rows."""

    TEXT = """Number of irreducible sector: 1
Lmax of AOs: 1
Lmax of ABFs: 2
Format:
The index of the symmetry operation
Spatial-only export: no orbital rotation matrices are written.

0
1.000000 0.000000 0.000000
0.000000 1.000000 0.000000
0.000000 0.000000 1.000000
(0.000000 0.000000 0.000000)
1
-1.000000 0.000000 0.000000
0.000000 -1.000000 0.000000
0.000000 0.000000 -1.000000
(0.500000 0.500000 0.500000)
Antiunitary operations: 1
0
0.000000 -1.000000 0.000000
1.000000 0.000000 0.000000
0.000000 0.000000 1.000000
(0.500000 0.500000 0.000000)
"""

    def test_two_unitary_and_one_antiunitary(self):
        unitary, antiunitary = read_symrot(self.TEXT)
        self.assertEqual([r for r, _ in unitary],
                         [[1, 0, 0, 0, 1, 0, 0, 0, 1],
                          [-1, 0, 0, 0, -1, 0, 0, 0, -1]])
        self.assertEqual(len(antiunitary), 1)
        self.assertEqual(antiunitary[0][0], [0, -1, 0, 1, 0, 0, 0, 0, 1])

    def test_translation_survives(self):
        _, antiunitary = read_symrot(self.TEXT)
        self.assertEqual(antiunitary[0][1], [0.5, 0.5, 0.0])
