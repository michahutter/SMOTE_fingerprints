import unittest

from rdkit import DataStructs
from rdkit.Chem import ChemicalFeatures
from rdkit.Chem.Pharm2D import Generate
from rdkit.Chem.Pharm2D.SigFactory import SigFactory
from rdkit.DataStructs import ExplicitBitVect

from source.code import fingerprints, matching, evaluation, compare
from source.code import combine
from source.code import distances
from source.code.ranks import *
from source.code import statistics
from source.code.globals import *


# fingerprints.py: finished
# distances.py: finished
# ranks.py: finished
# combine.py: finished
# matching.py: in progress
# compare.py: in progress
# evaluate.py: finished


class TestFunctions(unittest.TestCase):

    # -------------------------- fingerprints -----git ---------------------

    def test_read_compound_file_fail(self):

        self.assertIsNone(fingerprints.read_compound_file('thisfileshouldnotexistadsacsacsa'))
        self.assertIsNone(fingerprints.read_compound_file(''))
        self.assertIsNone(fingerprints.read_compound_file(1243543))
        self.assertIsNone(fingerprints.read_compound_file(1.3434))
        self.assertIsNone(fingerprints.read_compound_file(['common', 'fail']))

    def test_regenerate_fingerprint(self):
        m1 = combine.Chem.MolFromSmiles('C[C@H](N[C@@H](CCC1=CC=CC=C1)C(O)=O)C(=O)N2CCC[C@H]2C(O)=O')
        fdef_name = "source/general/BaseFeatures.fdef"

        fp_type1 = 'MACCS'
        fp_type2 = 'Morgan'
        fp_type3 = 'Atom Pairs'
        fp_type4 = 'Topological Torsion'
        fp_type5 = '2D Pharmacophore'

        fp1 = fingerprints.MACCSkeys.GenMACCSKeys(m1)
        fp2 = fingerprints.AllChem.GetMorganFingerprintAsBitVect(m1, 2)
        fp3 = fingerprints.Pairs.GetHashedAtomPairFingerprint(m1, nBits=2048)
        fp4 = fingerprints.Torsions.GetHashedTopologicalTorsionFingerprint(m1, nBits=2048)

        # Setting this up the same way I would in my personal fingerprint generation
        factory = fingerprints.ChemicalFeatures.BuildFeatureFactory(fdef_name)
        sig_factory = fingerprints.SigFactory(factory, trianglePruneBins=False)
        sig_factory.SetBins([(0, 2), (2, 5), (5, 8)])
        sig_factory.Init()
        sig_factory.GetSigSize()

        fp5 = fingerprints.Generate.Gen2DFingerprint(m1, sig_factory)

        fp1_e = fp1.ToBinary()
        fp2_e = fp2.ToBinary()
        fp3_e = fp3.ToBinary()
        fp4_e = fp4.ToBinary()
        fp5_e = fp5.ToBinary()

        res1 = fingerprints.regenerate_fingerprint(fp1_e, fp_type1)
        res2 = fingerprints.regenerate_fingerprint(fp2_e, fp_type2)
        res3 = fingerprints.regenerate_fingerprint(fp3_e, fp_type3)
        res4 = fingerprints.regenerate_fingerprint(fp4_e, fp_type4)
        res5 = fingerprints.regenerate_fingerprint(fp5_e, fp_type5)

        self.assertEqual(res1, fp1)
        self.assertEqual(res2, fp2)
        self.assertEqual(res3, fp3)
        self.assertEqual(res4, fp4)
        self.assertEqual(res5, fp5)

    def test_read_byte_fingerprint_file_fail(self):

        self.assertIsNone(fingerprints.read_byte_fingerprint_file(''))
        self.assertIsNone(fingerprints.read_byte_fingerprint_file(TEST_DATA_EMPTY_ASC))
        self.assertIsNone(fingerprints.read_byte_fingerprint_file(TEST_DATA_FAIL1_ASC))
        self.assertIsNone(fingerprints.read_byte_fingerprint_file(TEST_DATA_FAIL2))
        self.assertIsNone(fingerprints.read_byte_fingerprint_file(TEST_DATA_FAIL3))
        self.assertIsNone(fingerprints.read_byte_fingerprint_file(1234567))

    def test_read_byte_header_type(self):

        res1 = fingerprints.read_byte_header_type(TEST_DATA_ACTIVES1_FPS)
        res2 = fingerprints.read_byte_header_type(TEST_DATA_ACTIVES2_FPS)
        res22 = fingerprints.read_byte_header_type(TEST_DATA_ACTIVES22_FPS)
        res222 = fingerprints.read_byte_header_type(TEST_DATA_ACTIVES222_FPS)
        res3 = fingerprints.read_byte_header_type(TEST_DATA_ACTIVES3_FPS)
        res4 = fingerprints.read_byte_header_type(TEST_DATA_ACTIVES4_FPS)
        res5 = fingerprints.read_byte_header_type(TEST_DATA_ACTIVES5_FPS)
        res55 = fingerprints.read_byte_header_type(TEST_DATA_ACTIVES55_FPS)

        self.assertEqual(res1.lower(), 'MACCS'.lower())
        self.assertEqual(res2.lower(), 'Morgan'.lower())
        self.assertEqual(res22.lower(), 'Morgan'.lower())
        self.assertEqual(res222.lower(), 'Morgan'.lower())
        self.assertEqual(res3.lower(), 'Atom Pairs'.lower())
        self.assertEqual(res4.lower(), 'Topological Torsion'.lower())
        self.assertEqual(res5.lower(), '2D Pharmacophore'.lower())
        self.assertEqual(res55.lower(), '2D Pharmacophore'.lower())

    def test_check_fingerprint_byte_file(self):
        self.assertTrue(fingerprints.check_fingerprint_byte_file('compound name size MACCS\n'.encode()))
        self.assertTrue(fingerprints.check_fingerprint_byte_file('compound name size Morgan\n'.encode()))
        self.assertTrue(fingerprints.check_fingerprint_byte_file('compound name size Atom Pairs\n'.encode()))
        self.assertTrue(fingerprints.check_fingerprint_byte_file('compound name size Topological Torsion\n'.encode()))
        self.assertTrue(fingerprints.check_fingerprint_byte_file('compound name size 2D Pharmacophore\n'.encode()))

        self.assertFalse(fingerprints.check_fingerprint_byte_file('compound nam size MACCS\n'.encode()))
        self.assertFalse(fingerprints.check_fingerprint_byte_file('compound name size MAC\n'.encode()))
        self.assertFalse(fingerprints.check_fingerprint_byte_file('compound name'.encode()))

    def test_writing_and_reading_actives(self):

        morgan = '{"radius":2,"useFeatures":True}'
        morgan_incomplete = '{"useFeatures":True}'
        fdef_file = 'source/general/BaseFeatures.fdef'

        actives1 = fingerprints.create_fingerprints_output(TEST_DATA_ACTIVES,
                                                           'tests_output/actives1.asc', 'MACCS')

        actives2 = fingerprints.create_fingerprints_output(TEST_DATA_ACTIVES,
                                                           'tests_output/actives2.asc', 'Morgan')

        actives22 = fingerprints.create_fingerprints_output(TEST_DATA_ACTIVES,
                                                            'tests_output/actives22.asc', 'Morgan', morgan)

        actives222 = fingerprints.create_fingerprints_output(TEST_DATA_ACTIVES,
                                                             'tests_output/actives222.asc', 'Morgan', morgan_incomplete)

        actives3 = fingerprints.create_fingerprints_output(TEST_DATA_ACTIVES,
                                                           'tests_output/actives3.asc', 'Atom Pairs')

        actives4 = fingerprints.create_fingerprints_output(TEST_DATA_ACTIVES,
                                                           'tests_output/actives4.asc', 'Topological Torsion')

        actives5 = fingerprints.create_fingerprints_output(TEST_DATA_ACTIVES,
                                                           'tests_output/actives5.asc', '2D Pharmacophore')

        actives55 = fingerprints.create_fingerprints_output(TEST_DATA_ACTIVES,
                                                            'tests_output/actives55.asc', '2D Pharmacophore', fdef_file)

        actives_help1 = fingerprints.read_byte_fingerprint_file('tests_output/actives1.asc')
        actives_help2 = fingerprints.read_byte_fingerprint_file('tests_output/actives2.asc')
        actives_help22 = fingerprints.read_byte_fingerprint_file('tests_output/actives22.asc')
        actives_help222 = fingerprints.read_byte_fingerprint_file('tests_output/actives222.asc')
        actives_help3 = fingerprints.read_byte_fingerprint_file('tests_output/actives3.asc')
        actives_help4 = fingerprints.read_byte_fingerprint_file('tests_output/actives4.asc')
        actives_help5 = fingerprints.read_byte_fingerprint_file('tests_output/actives5.asc')
        actives_help55 = fingerprints.read_byte_fingerprint_file('tests_output/actives55.asc')

        for i in range(len(actives_help1[3])):
            self.assertEqual(actives_help1[3][i], actives1[2][i])

        for i in range(len(actives_help2[3])):
            self.assertEqual(actives_help2[3][i], actives2[2][i])

        for i in range(len(actives_help22[3])):
            self.assertEqual(actives_help22[3][i], actives22[2][i])

        for i in range(len(actives_help222[3])):
            self.assertEqual(actives_help222[3][i], actives222[2][i])

        # Should be same as well since same arguments given
        for i in range(len(actives_help222[3])):
            self.assertEqual(actives22[2][i], actives222[2][i])

        for i in range(len(actives_help3[3])):
            self.assertEqual(actives_help3[3][i], actives3[2][i])

        for i in range(len(actives_help4[3])):
            self.assertEqual(actives_help4[3][i], actives4[2][i])

        for i in range(len(actives_help5[3])):
            self.assertEqual(actives_help5[3][i], actives5[2][i])

        for i in range(len(actives_help55[3])):
            self.assertEqual(actives_help55[3][i], actives55[2][i])

    def test_writing_and_reading_inactives(self):
        morgan = '{"radius":2,"useFeatures":True}'
        morgan_incomplete = '{"useFeatures":True}'
        fdef_file = 'source/general/BaseFeatures.fdef'

        inactives1 = fingerprints.create_fingerprints_output(TEST_DATA_INACTIVES,
                                                             'tests_output/inactives1.asc', 'MACCS')

        inactives2 = fingerprints.create_fingerprints_output(TEST_DATA_INACTIVES,
                                                             'tests_output/inactives2.asc', 'Morgan')

        inactives22 = fingerprints.create_fingerprints_output(TEST_DATA_INACTIVES,
                                                              'tests_output/inactives22.asc', 'Morgan', morgan)

        inactives222 = fingerprints.create_fingerprints_output(TEST_DATA_INACTIVES,
                                                               'tests_output/inactives222.asc', 'Morgan',
                                                               morgan_incomplete)

        inactives3 = fingerprints.create_fingerprints_output(TEST_DATA_INACTIVES,
                                                             'tests_output/inactives3.asc', 'Atom Pairs')

        inactives4 = fingerprints.create_fingerprints_output(TEST_DATA_INACTIVES,
                                                             'tests_output/inactives4.asc', 'Topological Torsion')

        inactives5 = fingerprints.create_fingerprints_output(TEST_DATA_INACTIVES,
                                                             'tests_output/inactives5.asc',
                                                             '2D Pharmacophore')

        inactives55 = fingerprints.create_fingerprints_output(TEST_DATA_INACTIVES,
                                                              'tests_output/inactives55.asc', '2D Pharmacophore',
                                                              fdef_file)

        inactives_help1 = fingerprints.read_byte_fingerprint_file('tests_output/inactives1.asc')
        inactives_help2 = fingerprints.read_byte_fingerprint_file('tests_output/inactives2.asc')
        inactives_help22 = fingerprints.read_byte_fingerprint_file('tests_output/inactives22.asc')
        inactives_help222 = fingerprints.read_byte_fingerprint_file('tests_output/inactives222.asc')
        inactives_help3 = fingerprints.read_byte_fingerprint_file('tests_output/inactives3.asc')
        inactives_help4 = fingerprints.read_byte_fingerprint_file('tests_output/inactives4.asc')
        inactives_help5 = fingerprints.read_byte_fingerprint_file('tests_output/inactives5.asc')
        inactives_help55 = fingerprints.read_byte_fingerprint_file('tests_output/inactives55.asc')

        for i in range(len(inactives_help1[3])):
            self.assertEqual(inactives_help1[3][i], inactives1[2][i])

        for i in range(len(inactives_help2[3])):
            self.assertEqual(inactives_help2[3][i], inactives2[2][i])

        for i in range(len(inactives_help22[3])):
            self.assertEqual(inactives_help22[3][i], inactives22[2][i])

        for i in range(len(inactives_help222[3])):
            self.assertEqual(inactives_help222[3][i], inactives222[2][i])

        for i in range(len(inactives_help222[3])):
            self.assertEqual(inactives22[2][i], inactives222[2][i])

        for i in range(len(inactives_help3[3])):
            self.assertEqual(inactives_help3[3][i], inactives3[2][i])

        for i in range(len(inactives_help4[3])):
            self.assertEqual(inactives_help4[3][i], inactives4[2][i])

        for i in range(len(inactives_help5[3])):
            self.assertEqual(inactives_help5[3][i], inactives5[2][i])

        for i in range(len(inactives_help55[3])):
            self.assertEqual(inactives_help55[3][i], inactives55[2][i])

    def test_writing_and_reading2(self):
        morgan = '{"radius":2,"useFeatures":True}'
        morgan_incomplete = '{"useFeatures":True}'
        fdef_file = 'source/general/BaseFeatures.fdef'

        actives1 = fingerprints.create_fingerprints_output(TEST_DATA_ACTIVES,
                                                           'tests_output/actives1.asc', 'MACCS')

        actives2 = fingerprints.create_fingerprints_output(TEST_DATA_ACTIVES,
                                                           'tests_output/actives2.asc', 'Morgan')

        actives22 = fingerprints.create_fingerprints_output(TEST_DATA_ACTIVES,
                                                            'tests_output/actives22.asc', 'Morgan', morgan)

        actives222 = fingerprints.create_fingerprints_output(TEST_DATA_ACTIVES,
                                                             'tests_output/actives222.asc', 'Morgan', morgan_incomplete)

        actives3 = fingerprints.create_fingerprints_output(TEST_DATA_ACTIVES,
                                                           'tests_output/actives3.asc', 'Atom Pairs')

        actives4 = fingerprints.create_fingerprints_output(TEST_DATA_ACTIVES,
                                                           'tests_output/actives4.asc', 'Topological Torsion')

        actives5 = fingerprints.create_fingerprints_output(TEST_DATA_ACTIVES,
                                                           'tests_output/actives5.asc', '2D Pharmacophore')

        actives55 = fingerprints.create_fingerprints_output(TEST_DATA_ACTIVES,
                                                            'tests_output/actives55.asc', '2D Pharmacophore', fdef_file)

        actives_help1 = fingerprints.read_byte_fingerprint_file('tests_output/actives1.asc')
        actives_help2 = fingerprints.read_byte_fingerprint_file('tests_output/actives2.asc')
        actives_help22 = fingerprints.read_byte_fingerprint_file('tests_output/actives22.asc')
        actives_help222 = fingerprints.read_byte_fingerprint_file('tests_output/actives222.asc')
        actives_help3 = fingerprints.read_byte_fingerprint_file('tests_output/actives3.asc')
        actives_help4 = fingerprints.read_byte_fingerprint_file('tests_output/actives4.asc')
        actives_help5 = fingerprints.read_byte_fingerprint_file('tests_output/actives5.asc')
        actives_help55 = fingerprints.read_byte_fingerprint_file('tests_output/actives55.asc')

        # checking that they are actually unequal when they are supposed to
        for i in range(len(actives_help1[3]) - 1):
            self.assertNotEqual(actives_help1[3][i], actives1[2][i + 1])

        for i in range(len(actives_help2[3]) - 1):
            self.assertNotEqual(actives_help2[3][i], actives2[2][i + 1])

        for i in range(len(actives_help22[3]) - 1):
            self.assertNotEqual(actives_help22[3][i], actives22[2][i + 1])

        for i in range(len(actives_help222[3]) - 1):
            self.assertNotEqual(actives_help222[3][i], actives222[2][i + 1])

        for i in range(len(actives_help3[3]) - 1):
            self.assertNotEqual(actives_help3[3][i], actives3[2][i + 1])

        for i in range(len(actives_help4[3]) - 1):
            self.assertNotEqual(actives_help4[3][i], actives4[2][i + 1])

        for i in range(len(actives_help5[3]) - 1):
            self.assertNotEqual(actives_help5[3][i], actives5[2][i + 1])

        for i in range(len(actives_help55[3]) - 1):
            self.assertNotEqual(actives_help55[3][i], actives55[2][i + 1])

    def test_writing_and_reading_fail(self):
        actives1 = combine.create_fingerprints_output('notafileqwertzuio', 'tests_output/actives0.asc', 'MACCS')
        actives2 = combine.create_fingerprints_output(TEST_DATA_ACTIVES, '', 'MACCS')
        actives3 = combine.create_fingerprints_output(TEST_DATA_ACTIVES, 'tests_output/act1.asc', 'n')

        self.assertIsNone(actives1)
        self.assertIsNone(actives2)
        self.assertIsNone(actives3)

    # -------------------------- distances --------------------------

    def test_manhatten_simple(self):
        m1 = fingerprints.Chem.MolFromSmiles('C[C@H](N[C@@H](CCC1=CC=CC=C1)C(O)=O)C(=O)N2CCC[C@H]2C(O)=O')
        fdef_name = "source/general/BaseFeatures.fdef"

        fp1 = fingerprints.MACCSkeys.GenMACCSKeys(m1)
        fp2 = fingerprints.AllChem.GetMorganFingerprintAsBitVect(m1, 2)
        fp3 = fingerprints.Pairs.GetAtomPairFingerprintAsBitVect(m1)
        fp4 = fingerprints.Torsions.GetTopologicalTorsionFingerprintAsIntVect(m1)

        factory = fingerprints.ChemicalFeatures.BuildFeatureFactory(fdef_name)
        sig_factory = fingerprints.SigFactory(factory, trianglePruneBins=False)
        sig_factory.SetBins([(0, 2), (2, 5), (5, 8)])
        sig_factory.Init()
        sig_factory.GetSigSize()

        fp5 = fingerprints.Generate.Gen2DFingerprint(m1, sig_factory)

        self.assertEqual(distances.manhatten_distance(fp1, fp1), 0)
        self.assertEqual(distances.manhatten_distance(fp2, fp2), 0)
        self.assertEqual(distances.manhatten_distance(fp3, fp3), 0)
        self.assertEqual(distances.manhatten_distance(fp4, fp4), 0)
        self.assertEqual(distances.manhatten_distance(fp5, fp5), 0)

    def test_manhatten_check_MACCS_Morgan(self):
        fp1 = statistics.ExplicitBitVect(8)
        fp2 = statistics.ExplicitBitVect(8)
        self.assertEqual(distances.manhatten_distance(fp1, fp2), 0)

        fp1.SetBit(0)
        self.assertEqual(distances.manhatten_distance(fp1, fp2), 1)

        fp2.SetBit(1)
        self.assertEqual(distances.manhatten_distance(fp1, fp2), 2)

        fp1.SetBit(2)
        self.assertEqual(distances.manhatten_distance(fp1, fp2), 3)

        fp2.SetBit(2)
        self.assertEqual(distances.manhatten_distance(fp1, fp2), 2)

    def test_manhatten_check_values_Atom_pairs(self):
        fp1 = statistics.IntSparseIntVect(8)
        fp2 = statistics.IntSparseIntVect(8)
        self.assertEqual(distances.manhatten_distance(fp1, fp2), 0)

        fp1[0] = 2
        self.assertEqual(distances.manhatten_distance(fp1, fp2), 1)

        fp2[0] = 1
        self.assertEqual(distances.manhatten_distance(fp1, fp2), 0)

        fp2[1] = 1
        self.assertEqual(distances.manhatten_distance(fp1, fp2), 1)

    def test_manhatten_check_values_Topological_Torsion(self):
        fp1 = statistics.LongSparseIntVect(8)
        fp2 = statistics.LongSparseIntVect(8)
        self.assertEqual(distances.manhatten_distance(fp1, fp2), 0)

        fp1[0] = 2
        self.assertEqual(distances.manhatten_distance(fp1, fp2), 1)

        fp2[0] = 1
        self.assertEqual(distances.manhatten_distance(fp1, fp2), 0)

        fp2[1] = 1
        self.assertEqual(distances.manhatten_distance(fp1, fp2), 1)

    def test_manhatten_check_values_2D_Pharmacophore(self):
        fp1 = statistics.SparseBitVect(8)
        fp2 = statistics.SparseBitVect(8)
        self.assertEqual(distances.manhatten_distance(fp1, fp2), 0)

        fp1[0] = 2
        self.assertEqual(distances.manhatten_distance(fp1, fp2), 1)

        fp2[0] = 1
        self.assertEqual(distances.manhatten_distance(fp1, fp2), 0)

        fp2[1] = 1
        self.assertEqual(distances.manhatten_distance(fp1, fp2), 1)

    def test_manhatten_distance_fail(self):
        fp1 = statistics.ExplicitBitVect(7)
        fp11 = statistics.ExplicitBitVect(9)
        fp2 = statistics.LongSparseIntVect(7)
        fp22 = statistics.LongSparseIntVect(9)

        self.assertIsNone(distances.manhatten_distance(1, 2))
        self.assertIsNone(distances.manhatten_distance(fp1, fp11))
        self.assertIsNone(distances.manhatten_distance(fp2, fp22))

    def test_group_manhatten_simple(self):
        test1 = [['aa', 'ba', 20], ['aa', 'bb', 2], ['aa', 'bc', 1],
                 ['ab', 'ba', 33], ['ab', 'bb', 2], ['ac', 'ba', 42]]

        out1 = [['aa', [20, 2, 1]], ['ab', [33, 2]], ['ac', [42]]]

        self.assertEqual(group_manhatten(test1), out1)

    def test_group_manhatten_separated(self):
        test1 = [['aa', 'ba', 20], ['aa', 'bb', 2], ['aa', 'bc', 1],
                 ['ab', 'ba', 33], ['aa', 'bb', 2], ['ac', 'ba', 42]]

        out1 = [['aa', [20, 2, 1, 2]], ['ab', [33]], ['ac', [42]]]
        out2 = [['aa', [20, 2, 1]], ['ab', [33]], ['aa', [2]], ['ac', [42]]]

        self.assertNotEqual(group_manhatten(test1), out1)
        self.assertEqual(group_manhatten(test1), out2)

    def test_reading_writing_distances_simple(self):
        actives1 = fingerprints.read_byte_fingerprint_file(TEST_DATA_ACTIVES1_FPS)
        inactives1 = fingerprints.read_byte_fingerprint_file(TEST_DATA_INACTIVES1_FPS)
        actives2 = fingerprints.read_byte_fingerprint_file(TEST_DATA_ACTIVES2_FPS)
        inactives2 = fingerprints.read_byte_fingerprint_file(TEST_DATA_INACTIVES2_FPS)
        actives3 = fingerprints.read_byte_fingerprint_file(TEST_DATA_ACTIVES3_FPS)
        inactives3 = fingerprints.read_byte_fingerprint_file(TEST_DATA_INACTIVES3_FPS)
        actives4 = fingerprints.read_byte_fingerprint_file(TEST_DATA_ACTIVES4_FPS)
        inactives4 = fingerprints.read_byte_fingerprint_file(TEST_DATA_INACTIVES4_FPS)
        actives5 = fingerprints.read_byte_fingerprint_file(TEST_DATA_ACTIVES5_FPS)
        inactives5 = fingerprints.read_byte_fingerprint_file(TEST_DATA_INACTIVES5_FPS)

        distances.create_distances_output(actives1[3], inactives1[3], actives1[1], inactives1[1],
                                          'tests_output/distances1.asc')
        distances.create_distances_output(actives2[3], inactives2[3], actives2[1], inactives2[1],
                                          'tests_output/distances2.asc')
        distances.create_distances_output(actives3[3], inactives3[3], actives3[1], inactives3[1],
                                          'tests_output/distances3.asc')
        distances.create_distances_output(actives4[3], inactives4[3], actives4[1], inactives4[1],
                                          'tests_output/distances4.asc')
        distances.create_distances_output(actives5[3], inactives5[3], actives5[1], inactives5[1],
                                          'tests_output/distances5.asc')

    # -------------------------- ranking --------------------------

    def test_ranking_manhatten_simple(self):
        test1 = [['aa', 'ba', 200], ['aa', 'bb', 2], ['aa', 'bc', 5],
                 ['ab', 'ba', 32], ['ab', 'bb', 2], ['ac', 'ba', 42]]

        out1 = [['aa', 207], ['ac', 42], ['ab', 34]]

        self.assertEqual(ranking_manhatten(test1), out1)

    def test_create_ranking_output_simple(self):

        distances1 = distances.read_distances(TEST_DATA_DISTS1)
        distances2 = distances.read_distances(TEST_DATA_DISTS2)
        distances3 = distances.read_distances(TEST_DATA_DISTS3)
        distances4 = distances.read_distances(TEST_DATA_DISTS4)
        distances5 = distances.read_distances(TEST_DATA_DISTS5)

        ranking1 = ranking_manhatten(distances1)
        ranking2 = ranking_manhatten(distances2)
        ranking3 = ranking_manhatten(distances3)
        ranking4 = ranking_manhatten(distances4)
        ranking5 = ranking_manhatten(distances5)

        create_ranking_output(ranking1, 'tests_output/ranks1.asc')
        create_ranking_output(ranking2, 'tests_output/ranks2.asc')
        create_ranking_output(ranking3, 'tests_output/ranks3.asc')
        create_ranking_output(ranking4, 'tests_output/ranks4.asc')
        create_ranking_output(ranking5, 'tests_output/ranks5.asc')

    def test_read_ranking_output_simple(self):
        ranks = read_ranking_output(TEST_DATA_ALTER)

        self.assertEqual(len(ranks), 5)
        self.assertEqual(ranks[0][0], 'TOP1')
        self.assertEqual(ranks[1][0], 'TOP2')
        self.assertEqual(ranks[2][0], 'TOP3')
        self.assertEqual(ranks[3][0], 'TOP4')
        self.assertEqual(ranks[4][0], 'TOP5')
        self.assertEqual(ranks[0][1], '50')
        self.assertEqual(ranks[1][1], '40')
        self.assertEqual(ranks[2][1], '30')
        self.assertEqual(ranks[3][1], '20')
        self.assertEqual(ranks[4][1], '10')

    # -------------------------- combine --------------------------

    def test_combine_fingerprints_simple(self):
        fp1 = statistics.ExplicitBitVect(8)
        fp2 = statistics.IntSparseIntVect(7)
        fp3 = statistics.LongSparseIntVect(6)
        fp4 = statistics.SparseBitVect(5)

        res1 = combine.combine_fingerprints([fp1])
        self.assertEqual(len(res1[0]), 8)

        res2 = combine.combine_fingerprints([fp2])
        self.assertEqual(res2[0].GetLength(), 7)

        res3 = combine.combine_fingerprints([fp3])
        self.assertEqual(res3[0].GetLength(), 6)

        res4 = combine.combine_fingerprints([fp4])
        self.assertEqual(len(res4[0]), 5)

    def test_combine_fingerprints_Explicit_Bit_Vect(self):
        fp1 = statistics.ExplicitBitVect(8)
        fp2 = statistics.ExplicitBitVect(8)

        fp1[1], fp1[2], fp1[3], fp2[2], fp2[4] = 1, 1, 1, 1, 1

        # fp1 = [0, 1, 1, 1, 0, 0, 0, 0]
        # fp2 = [0, 0, 1, 0, 1, 0, 0, 0]

        res3 = combine.combine_fingerprints([fp1, fp2])
        self.assertEqual(len(res3[0]), 8)
        self.assertEqual(res3[1], 8)
        self.assertEqual(res3[0].ToList(), [0, 1, 1, 1, 1, 0, 0, 0])

    def test_combine_fingerprints_Int_Sparse_Int_Vect(self):
        fp1 = statistics.IntSparseIntVect(8)
        fp2 = statistics.IntSparseIntVect(8)

        fp1[1], fp1[2], fp1[3], fp2[2], fp2[4] = 1, 1, 1, 1, 1

        # fp1 = [0, 1, 1, 1, 0, 0, 0, 0]
        # fp2 = [0, 0, 1, 0, 1, 0, 0, 0]

        res3 = combine.combine_fingerprints([fp1, fp2])
        self.assertEqual(res3[0].GetLength(), 8)
        self.assertEqual(res3[1], 8)
        self.assertEqual(res3[0].ToList(), [0, 1, 1, 1, 1, 0, 0, 0])

    def test_combine_fingerprints_Long_Sparse_Int_Vect(self):
        fp1 = statistics.LongSparseIntVect(8)
        fp2 = statistics.LongSparseIntVect(8)

        fp1[1], fp1[2], fp1[3], fp2[2], fp2[4] = 1, 1, 1, 1, 1

        # fp1 = [0, 1, 1, 1, 0, 0, 0, 0]
        # fp2 = [0, 0, 1, 0, 1, 0, 0, 0]

        res3 = combine.combine_fingerprints([fp1, fp2])
        self.assertEqual(res3[0].GetLength(), 8)
        self.assertEqual(res3[1], 8)
        self.assertEqual(res3[0].ToList(), [0, 1, 1, 1, 1, 0, 0, 0])

    def test_combine_fingerprints_Sparse_Bit_Vect(self):
        fp1 = statistics.SparseBitVect(8)
        fp2 = statistics.SparseBitVect(8)

        fp1[1], fp1[2], fp1[3], fp2[2], fp2[4] = 1, 1, 1, 1, 1

        # fp1 = [0, 1, 1, 1, 0, 0, 0, 0]
        # fp2 = [0, 0, 1, 0, 1, 0, 0, 0]

        res3 = combine.combine_fingerprints([fp1, fp2])
        self.assertEqual(len(res3[0]), 8)
        self.assertEqual(res3[1], 8)
        self.assertEqual(res3[0].ToList(), [0, 1, 1, 1, 1, 0, 0, 0])

    def test_get_top_fingerprints(self):
        ranks = read_ranking_output(TEST_DATA_RANKS_SHORT)
        fps = combine.read_byte_fingerprint_file(TEST_DATA_ACTIVES_SHORT)

        data = combine.get_top_fingerprints(ranks, fps, 4)

        self.assertEqual(len(ranks), 4)
        self.assertEqual(data[0][0], 'benazepril')
        self.assertEqual(data[1][0], 'enalapril')
        self.assertEqual(data[2][0], 'spiraprilat')
        self.assertEqual(data[3][0], 'enalaprilat')
        self.assertEqual(data[0][1], fps[3][3])
        self.assertEqual(data[1][1], fps[3][1])
        self.assertEqual(data[2][1], fps[3][2])
        self.assertEqual(data[3][1], fps[3][0])

    def test_get_k_nearest_fingerprints(self):
        ranks = read_ranking_output(TEST_DATA_RANKS_SHORT)
        fps = combine.read_byte_fingerprint_file(TEST_DATA_ACTIVES_SHORT)

        data = combine.get_k_nearest_fingerprints(ranks, fps, 3)

        self.assertEqual(len(ranks), 4)
        self.assertEqual(data[0][0], 'benazepril')
        self.assertEqual(data[1][0], 'enalapril')
        self.assertEqual(data[2][0], 'enalaprilat')
        self.assertEqual(data[3][0], 'spiraprilat')
        self.assertEqual(data[0][1], fps[3][3])
        self.assertEqual(data[1][1], fps[3][1])
        self.assertEqual(data[2][1], fps[3][0])
        self.assertEqual(data[3][1], fps[3][2])
        self.assertGreater(data[0][2], data[1][2])
        self.assertGreater(data[1][2], data[2][2])
        self.assertGreater(data[2][2], data[3][2])

    def test_combine_fingerprints_threshold(self):
        fp1 = statistics.ExplicitBitVect(8)
        fp2 = statistics.ExplicitBitVect(8)
        fp3 = statistics.ExplicitBitVect(8)
        fp4 = statistics.ExplicitBitVect(8)

        fp1[0], fp2[0], fp3[0], fp4[0] = 1, 1, 1, 1
        fp1[1], fp2[1], fp3[1] = 1, 1, 1
        fp1[2], fp2[2] = 1, 1
        fp1[3] = 1

        res_fp1 = combine.combine_fingerprints([fp1, fp2, fp3, fp4], 0.9)
        res_fp2 = combine.combine_fingerprints([fp1, fp2, fp3, fp4], 0.75)
        res_fp3 = combine.combine_fingerprints([fp1, fp2, fp3, fp4], 0.5)
        res_fp4 = combine.combine_fingerprints([fp1, fp2, fp3, fp4], 0.25)
        res_fp5 = combine.combine_fingerprints([fp1, fp2, fp3, fp4], 1)

        self.assertEqual(res_fp1[0].ToList(), [1, 0, 0, 0, 0, 0, 0, 0])
        self.assertEqual(res_fp2[0].ToList(), [1, 1, 0, 0, 0, 0, 0, 0])
        self.assertEqual(res_fp3[0].ToList(), [1, 1, 1, 0, 0, 0, 0, 0])
        self.assertEqual(res_fp4[0].ToList(), [1, 1, 1, 1, 0, 0, 0, 0])
        self.assertEqual(res_fp5[0].ToList(), [1, 0, 0, 0, 0, 0, 0, 0])


    def test_create_combine_output_simple(self):
        ranks = read_ranking_output(TEST_DATA_RANKS1)
        fps = combine.read_byte_fingerprint_file(TEST_DATA_ACTIVES1_FPS)
        fp_type = combine.read_byte_header_type(TEST_DATA_ACTIVES1_FPS)

        self.assertIsNotNone(combine.create_combine_output(ranks, fps, fp_type,
                                                           "tests_output/combine1.asc", "top", 2))
        self.assertIsNotNone(combine.create_combine_output(ranks, fps, fp_type,
                                                           "tests_output/combine2.asc", "nearest", 2))
        self.assertIsNotNone(combine.create_combine_output(ranks, fps, fp_type,
                                                           "tests_output/combine1.asc", "top", 2))
        self.assertIsNotNone(combine.create_combine_output(ranks, fps, fp_type,
                                                           "tests_output/combine2.asc", "nearest", 2))
        self.assertIsNone(combine.create_combine_output(ranks, fps, fp_type,
                                                        "tests_output/combine2.asc", "neart", 2))

    def test_write_read_combine(self):
        ranks = read_ranking_output(TEST_DATA_ALTER)
        fps = list()

        fp1 = ExplicitBitVect(8)
        fp2 = ExplicitBitVect(8)
        fp3 = ExplicitBitVect(8)
        fp4 = ExplicitBitVect(8)
        fp5 = ExplicitBitVect(8)

        fp1[0], fp1[1], fp1[2], fp1[3] = 1, 1, 1, 1
        fp2[0], fp2[1], fp2[2] = 1, 1, 1
        fp3[0], fp3[1] = 1, 1
        fp4[0] = 1

        fps.append(['C', 'TOP1', '8', fp1])
        fps.append(['CC', 'TOP2', '8', fp2])
        fps.append(['CCC', 'TOP3', '8', fp3])
        fps.append(['CCCC', 'TOP4', '8', fp4])
        fps.append(['CCCCC', 'TOP5', '8', fp5])

        fps = list(map(list, zip(*fps)))

        combine.create_combine_output(ranks, fps, 'maccs', "tests_output/combine_test.asc", "top", 5)
        data = combine.read_combine("tests_output/combine_test.asc")

        fp = combine.regenerate_fingerprint(data[0][1], data[0][0])

        self.assertEqual(data[1][0], 'TOP1')
        self.assertEqual(data[1][1], 'TOP2')
        self.assertEqual(data[1][2], 'TOP3')
        self.assertEqual(data[1][3], 'TOP4')
        self.assertEqual(data[1][4], 'TOP5')
        self.assertEqual(fp, fp3)

    # -------------------------- matching --------------------------

    def test_folding_simple(self):
        m1 = combine.Chem.MolFromSmiles('C[C@H](N[C@@H](CCC1=CC=CC=C1)C(O)=O)C(=O)N2CCC[C@H]2C(O)=O')

        fp_type1 = 'MACCS'
        fp_type2 = 'Morgan'
        fp_type3 = 'Atom Pairs'
        fp_type4 = 'Topological Torsion'
        fp_type5 = '2D Pharmacophore'

        fp1 = fingerprints.MACCSkeys.GenMACCSKeys(m1)
        fp2 = fingerprints.AllChem.GetMorganFingerprintAsBitVect(m1, 2)
        fp3 = fingerprints.Pairs.GetHashedAtomPairFingerprint(m1, nBits=2048)
        fp4 = fingerprints.Torsions.GetHashedTopologicalTorsionFingerprint(m1, nBits=2048)

        factory = ChemicalFeatures.BuildFeatureFactory(PHARMACOPHORE_FEATURES)
        sig_factory = SigFactory(factory, trianglePruneBins=False)
        sig_factory.SetBins([(0, 2), (2, 5), (5, 8)])
        sig_factory.Init()
        sig_factory.GetSigSize()

        fp5 = Generate.Gen2DFingerprint(m1, sig_factory)

        test = ExplicitBitVect(fp3.GetLength())
        test.SetBitsFromList(list(fp3.GetNonzeroElements().keys()))
        test2 = ExplicitBitVect(fp4.GetLength())
        test2.SetBitsFromList(list(fp4.GetNonzeroElements().keys()))

        fp3 = test
        fp4 = test2

        f_fp1 = DataStructs.FoldFingerprint(fp1, 8)
        f_fp2 = DataStructs.FoldFingerprint(fp2, 8)
        f_fp3 = DataStructs.FoldFingerprint(fp3, 8)
        f_fp4 = DataStructs.FoldFingerprint(fp4, 8)
        f_fp5 = DataStructs.FoldFingerprint(fp5, 8)

        self.assertEqual(len(f_fp1), 20)
        self.assertEqual(len(f_fp2), 256)
        self.assertEqual(len(f_fp3), 256)
        self.assertEqual(len(f_fp4), 256)
        self.assertEqual(len(f_fp5), 418)

    # -------------------------- compare --------------------------

    # -------------------------- evaluation --------------------------

    def test_barchart_sanity(self):

        data = matching.read_all_scores(TEST_DATA_BARCHART)
        evaluation.barchart(data, 'tests_output/barchart.png')

    def test_create_evaluation_output(self):

        group1 = compare.read_compare_file(TEST_DATA_COMPARE1)
        group2 = compare.read_compare_file(TEST_DATA_COMPARE2)
        group3 = compare.read_compare_file(TEST_DATA_COMPARE3)

        self.assertIsNotNone(group1)
        self.assertIsNotNone(group2)
        self.assertIsNotNone(group3)

        res1 = evaluation.create_evaluation_output(group1, group2, TEST_DATA_COMPARE1,
                                                   TEST_DATA_COMPARE2, 'tests_output/eval1.asc')
        res2 = evaluation.create_evaluation_output(group3, group1, TEST_DATA_COMPARE1,
                                                   TEST_DATA_COMPARE3, 'tests_output/eval2.asc')

        # no difference
        self.assertEqual(res1, 0)
        # difference
        self.assertEqual(res2, 1)

if __name__ == '__main__':
    # Suppressing output of tested functions, using buffer.
    unittest.main(buffer=True)


