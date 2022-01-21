"""
Unit Tests for the scientific_data_challenge.py.

There are a total of eight unit tests covering almost all the methods in the subsequent code.
Prerequisites to run this code:
a. UnitTest - The unittest unit testing framework, downloaded directly from Pycharm. The pattern of the UnitTest will
follow the AAA (ARRANGE-ACT-ASSERT) format.
b. Pandas - open source data analysis and manipulation tool. Installation through general Pycharm --> Preferences -->
Project Interpreter --> install pandas.
"""
import unittest
from pandas.testing import assert_series_equal
import pandas as pd
from scientific_data_challenge import ScientificData
import filecmp


class TestDemo(unittest.TestCase):

    def test_dataframe(self):
        """
        UnitTest case to check the shape (rows and columns) of an sdf imported as a dataframe.
        """
        # arrange
        sdf_obj = ScientificData("Drugs.sdf")
        # act
        actual_count_row = sdf_obj.sdf_df.shape[0]
        actual_count_columns = sdf_obj.sdf_df.shape[1]
        # assert
        assert (actual_count_row == 10) and (actual_count_columns == 34),\
            "The shape (rows and columns) of the dataframe is not as expected"

    def test_sdf_sanity(self):
        """
        UnitTest case to check the sanity of an uploaded sdf.
        """
        with self.assertRaises(Exception) as ex:
            ScientificData("empty.txt")

        print(ex.exception)

    def test_smiles_in_sdf(self):
        """
        UnitTest case to check the SMILES returned for an sdf matches the expected SMILES.
        """
        # arrange
        sdf_obj = ScientificData("Drugs.sdf")
        # act
        actual_smiles = sdf_obj.sdf_df["SMILES"].tolist()
        expected_smiles_list = ["Nc1nc2[nH]nnc2c(=O)[nH]1", "N=C(N)NS(=O)(=O)c1ccc(N)cc1",
                                "COC(=O)Nc1nc2cc(C(=O)c3cccs3)ccc2[nH]1",
                                "CN(C)[C@H]1C(O)=C(C(=O)NCNC(CCCCN)C(=O)O)C(=O)[C@]2(O)C(O)=C3C(=O)c4c(O)cccc4[C@](C)(O)[C@@H]3C[C@H]12",
                                "CC(=O)OCC(=O)[C@H]1CC[C@H]2[C@@H]3CC[C@H]4C[C@H](O)CC[C@]4(C)[C@H]3C(=O)C[C@]12C",
                                "CC(=O)[C@H]1CC[C@H]2[C@@H]3CC[C@H]4C[C@H](O)CC[C@]4(C)[C@H]3C(=O)C[C@]12C",
                                "CCC[C@@H]1C(=O)N2C(N(C)C)=Nc3ccc(C)cc3N2C1=O",
                                "CC[C@]1(c2cccc(O)c2)CCCCN(C)C1.Cl",
                                "CN[C@@H]1C(O[C@H]2O[C@H](CO)[C@@H](N)[C@H](O)[C@H]2O)O[C@H]2C[C@@H](N)[C@@H](O[C@@H]3[C@@H](N)C[C@@H](N)[C@H](O)[C@H]3O)O[C@@H]2[C@@H]1O",
                                "C/C(=C(/CCO)SSC[C@H]1CCCO1)N(C=O)Cc1cnc(C)nc1N.Cl"]
        # assert
        assert actual_smiles == expected_smiles_list, \
            "The expected SMILES strings in the sdf {} does not match the actual SMILES string {}". \
            format(expected_smiles_list, actual_smiles)

    def test_molecular_weight(self):
        """
        UnitTest case to check the Molecular weight returned for an sdf matches the expected Molecular weight.
        """
        # arrange
        test = ScientificData("test.sdf")
        test.calc_descriptors()
        # act
        expected_mol_wt_from_empirical_formula = round(((12.0109 * 9) + (1.00784 * 13) + 14.0067), 2)
        # assert
        assert round(test.sdf_df.at[0, "Molwt"], 2) == expected_mol_wt_from_empirical_formula, \
            "The expected Molecular weights in the sdf {} does not match the actual actual string {}". \
            format(expected_mol_wt_from_empirical_formula, test.sdf_df["Molwt"])

    def test_logp(self):
        """
        UnitTest case to check if the logP value returned for an sdf matches the expected SMILES.
        """
        # arrange
        test = ScientificData("test.sdf")
        test.calc_descriptors()

        # assert
        assert round(test.sdf_df.at[0, "MolLogP"], 2) == 2.1, \
            "The expected LogP values in the sdf {} does not match the actual actual LogP values {}". \
            format(2.1, test.sdf_df["MolLogP"])

    def test_iupac_name_using_pubchempy(self):
        """
        UnitTest case to check if the IUPAC name returned for the first compound in an sdf (using pubchempy) matches the
        expected iupac name.
        """
        # arrange
        test = ScientificData("test.sdf")
        actual_smiles_iupac_dict = test.iupac_name_generator()
        # act
        expected_smiles_iupac_dictionary = {'CCC(N)c1ccccc1': '1-phenylpropan-1-amine'}
        # assert
        assert expected_smiles_iupac_dictionary == actual_smiles_iupac_dict, \
            "The expected IUPAC names for the first compound in the sdf {} does not match the actual IUPAC string {}". \
            format(expected_smiles_iupac_dictionary, actual_smiles_iupac_dict)

    def test_mol_images(self):
        """
        UnitTest case to check if the images returned from a file is as expected.
        """
        # arrange
        test = ScientificData("Drugs.sdf")
        test.mol_image()
        # assert
        self.assertTrue(filecmp.cmp("./output_image.png", "./output.png", shallow=False))

    def test_substructure_search(self):
        """
        UnitTest case to check if the substructure search for a compound returns expected hits for an sdf.
        """
        # arrange
        test = ScientificData('Drugs.sdf')
        subs_df = test.substructure_search("C[C@H]1CC[C@H]2[C@@H]3CC[C@H]4C[C@H](O)CC[C@]4(C)[C@H]3CC[C@]12C")

        # act
        expected_subs_match = ["CC(=O)OCC(=O)[C@H]1CC[C@H]2[C@@H]3CC[C@H]4C[C@H](O)CC[C@]4(C)[C@H]3C(=O)C[C@]12C",
             "CC(=O)[C@H]1CC[C@H]2[C@@H]3CC[C@H]4C[C@H](O)CC[C@]4(C)[C@H]3C(=O)C[C@]12C"]
        # assert
        assert subs_df["SMILES"].tolist() == expected_subs_match, \
            "The expected substructures for the sdf does not match the actual substructures {}". \
            format(expected_subs_match, subs_df["SMILES"].tolist())

    def test_lipinski_rule(self):
        """
        UnitTest case to test if a molecule passes the Lipinski filters and returns True or False accordingly.
        """
        # arrange
        test = ScientificData('Drugs.sdf')
        df2 = test.lipinski_rule()
        # act
        expected_lipinski_results = pd.Series([True, True, True, False, True, True, True, True, False, True],
                                              dtype=bool)
        # assert
        assert_series_equal(df2, expected_lipinski_results), \
            "The expected booleans for Lipinski Rule for the sdf does not match the actual defined booleans {}". \
                format(expected_lipinski_results, df2)

    def test_stereochemistry_checker(self):
        """
        UnitTest case to test the count of number of molecules which are stereo-active and stereo-inactive in an sdf.
        """
        # arrange
        test = ScientificData('Drugs.sdf')
        actual_dict_with_stereoinfo, actual_dict_with_stereo_number = test.stereochemistry_checker()
        # act
        expected_dict_with_stereo_number = {'StereoActive': 7, 'Stereo_Inactive': 3}
        expected_stereo_config_for_one_compound = [(2, 'R')]
        # asserts
        assert actual_dict_with_stereoinfo["CC[C@]1(c2cccc(O)c2)CCCCN(C)C1.Cl"][1] == expected_stereo_config_for_one_compound
        assert actual_dict_with_stereo_number == expected_dict_with_stereo_number, \
            "The expected counts for stereo_active and stereo_inactive molecules in the sdf does not match the actual " \
            "counts"

    def test_sorting_problems(self):
        """
        UnitTest case to test two things:
        1. Returns the exact number of rows for a particular therapeutic class.
        2. The Molwt. is sorted for that particular therapeutic class.
        """
        # arrange
        test = ScientificData("Drugs.sdf")
        filtered_df = test.sorting_problems("Metabolism")
        row_count = filtered_df.shape[0]
        # act
        expected_mol_weight = pd.Series([214.250, 300.362, 435.015, 539.583, 602.641], index=[1, 6, 9, 8, 3],
                                        name="Molwt")
        # assert
        assert (row_count, 5)
        assert_series_equal(filtered_df["Molwt"], expected_mol_weight), \
            "The expected molecular weights of a particular therapeutic class {} does not match the actual molecular " \
            "weights {}". format(expected_mol_weight, filtered_df["Molwt"])


if __name__ == '__main__':
    unittest.main()
