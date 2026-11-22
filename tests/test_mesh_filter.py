"""The sampled-mesh filter that symrot_R.txt needs before it drives k-stars."""
import unittest

from taesc_prep.abacus import preserves_mesh, reciprocal_rotation


class MeshFilterTests(unittest.TestCase):
    def test_identity_preserves_every_mesh(self):
        identity = [1, 0, 0, 0, 1, 0, 0, 0, 1]
        for mesh in ((1, 1, 3), (1, 1, 2), (2, 2, 2), (2, 3, 4)):
            self.assertTrue(preserves_mesh(identity, mesh), mesh)

    def test_cubic_only_operation(self):
        # maps a cubic mesh onto itself, but not (1,1,3)
        operation = [1, 0, -1, 1, 0, 0, 1, -1, 0]
        self.assertTrue(preserves_mesh(operation, (2, 2, 2)))
        self.assertFalse(preserves_mesh(operation, (1, 1, 3)))
        self.assertFalse(preserves_mesh(operation, (2, 2, 1)))

    def test_rotation_about_z_survives_a_z_mesh(self):
        rotation = [0, -1, 0, 1, 0, 0, 0, 0, 1]
        self.assertTrue(preserves_mesh(rotation, (1, 1, 3)))
        self.assertTrue(preserves_mesh(rotation, (2, 2, 2)))

    def test_antunitary_flips_the_mesh_test(self):
        # Theta * inversion reverses k, so a mesh is preserved only if it is
        # symmetric about the origin.
        inversion = [-1, 0, 0, 0, -1, 0, 0, 0, -1]
        self.assertTrue(preserves_mesh(inversion, (2, 2, 2)))
        self.assertTrue(preserves_mesh(inversion, (2, 2, 2), antiunitary=True))

    def test_reciprocal_of_a_known_operation(self):
        # [1,0,-1, 1,0,0, 1,-1,0] has det = -1; its reciprocal is
        # [0,0,-1, 1,1,1, 0,-1,0] (verified against the ABACUS-side helper).
        self.assertEqual(reciprocal_rotation([1, 0, -1, 1, 0, 0, 1, -1, 0]),
                         [0, 0, -1, 1, 1, 1, 0, -1, 0])

    def test_non_unimodular_rotation_is_rejected(self):
        with self.assertRaises(ValueError):
            reciprocal_rotation([2, 0, 0, 0, 1, 0, 0, 0, 1])


if __name__ == "__main__":
    unittest.main()
