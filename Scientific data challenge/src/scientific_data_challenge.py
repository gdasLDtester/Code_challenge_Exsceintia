"""
Code challenge round for Scientific Data Engineer position in Exscientia. Refer to the README file on instructions on
how to execute this script.
"""
from rdkit import Chem
import pubchempy
from rdkit.Chem import Draw, PandasTools
from rdkit.Chem.Descriptors import MolWt, MolLogP, NumAromaticRings, NumHDonors, NumHAcceptors
import argparse
import sys


def get_arguments(argv):
    """
    :param argv: list, the list of arguments to be parsed
    :return: list, the parsed arguments
    """
    parser = argparse.ArgumentParser(usage='python {} --sdf <SDF Path> [optional arguments]'.format(sys.argv[0]))
    required = parser.add_argument_group('required arguments')
    required.add_argument('--sdf', type=str,  help='Path of the sdf file in the directory', required=True)
    parser.add_argument('--calc_descriptors', action='store_true', help='gives the user the list of descriptors like '
                                                                        'Molwt, MolLogP, Num Aromatic Rings, HBA, HBD')
    parser.add_argument('--image', action='store_true', help='Gives the images of a set of molecules in an sdf')
    parser.add_argument('--iupac_name', action='store_true', help='Returns the IUPAC names of all compounds in an sdf '
                                                                  'using PubChempy')
    parser.add_argument('--substructure_search', type=str, help='Gives the list of SMILES')
    parser.add_argument('--lipinski_rule', action='store_true', help='Returns the IUPAC names of all compounds in an '
                                                                     'sdf using PubChempy')
    parser.add_argument('--sorted_sdf', type=str, help='Gives a sorted list of Molwt based on a particular therapeutic class')
    parser.add_argument('--stereochemistry_checker', action='store_true', help='Returns the IUPAC names of all '
                                                                               'compounds in an sdf using PubChempy')

    return parser.parse_args(argv)


class ScientificData:

    # instance variable for sdf dataframe
    sdf_df = None

    def __init__(self, sdf_file_path):
        # converting an sdf into a dataframe
        self.sdf_df = self.loading_sdf(sdf_file_path)

    @staticmethod
    def loading_sdf(sdf_file_path):
        """
        Method to invoke an sdf file given its path and returns a dataframe of rows and columns.

        :param sdf_file_path:str, path of the file in your directory
        :return df: DataFrame, 2-D table of rows and columns containing data
        """
        # Using PandasTools library to read file in SDF format and load into a dataframe
        df = PandasTools.LoadSDF(sdf_file_path, smilesName='SMILES', molColName='Molecule')

        # Check if the SDF is invalid and an empty dataframe is returned
        if df.empty:
            raise Exception("Invalid SDF")
        return df

    def calc_descriptors(self):
        """
        Method to calculate the descriptors of compounds in an sdf using RDKit. The descriptors calculated are 'Molwt',
        'MolLogP', 'Num Aromatic Rings', 'HBA', 'HBD"""
        # Calculating descriptors using RDKit descriptor calculation by iterating over the molecule object
        descriptors = ['Molwt', 'MolLogP', 'Num Aromatic Rings', 'HBA', 'HBD']
        mol_wt, logp, num_arom_rings, num_h_acceptors, num_h_donors = [], [], [], [], []
        for mol in self.sdf_df["Molecule"]:
            mol_wt.append(MolWt(mol))
            logp.append(MolLogP(mol))
            num_arom_rings.append(NumAromaticRings(mol))
            num_h_acceptors.append(NumHAcceptors(mol))
            num_h_donors.append(NumHDonors(mol))

        # Assigning values of these descriptors wrt to the column names defined
        self.sdf_df[descriptors[0]] = mol_wt
        self.sdf_df[descriptors[1]] = logp
        self.sdf_df[descriptors[2]] = num_arom_rings
        self.sdf_df[descriptors[3]] = num_h_acceptors
        self.sdf_df[descriptors[4]] = num_h_donors

    def mol_image(self):
        """
        Method to give the image of an molecule with stereo-annotations included - both absolute and enhanced
        stereo-labels.
        """
        # Calling RDKit.Draw.MolDrawOptions to add Stereo-annotation to compound images.
        dos = Draw.MolDrawOptions()
        dos.addStereoAnnotation = True

        # Generating images in a grid and saving them in a file
        img = Draw.MolsToGridImage(self.sdf_df["Molecule"], molsPerRow=4, subImgSize=(800, 800), drawOptions=dos)
        img.save("output_image.png")

    def iupac_name_generator(self):
        """
        Method to find out the IUPAC name of compounds in an sdf from the SMILES using the PubChemPy package.

        :return struct_iupac_name: dict, dictionary containing the SMILES and their corresponding IUPAC name
        """
        struct_iupac_names = {}
        # Iterating over the smiles and using the PubChemPy library to generate IUPAC names.
        for smiles in self.sdf_df["SMILES"]:
            compounds = pubchempy.get_compounds(smiles, namespace='smiles')
            match = compounds[0]

            # Adding SMILES as keys and IUPAC names as values in the dictionary
            struct_iupac_names[smiles] = match.iupac_name
        return struct_iupac_names

    def substructure_search(self, substructure_smiles):
        """
        Method to find out the list of SMILES that match a given substructure SMILES in an sdf.

        :param substructure_smiles: str, the SMILES string on which a substructure search will be performed.
        :return subs_match: list of str, the SMILES that match the substructure.
        """
        # Converting the smiles passed as a mol object
        steroid_scaffold = Chem.MolFromSmiles(substructure_smiles)

        # Using the >= operator to perform a substructure search against the molecule object
        steroid_scaffold_comps = self.sdf_df[self.sdf_df['Molecule'] >= steroid_scaffold]

        # Storing the matched smiles as a dataframe with IUPAC names and SMILES
        subs_match = steroid_scaffold_comps[["name IUPAC", "SMILES"]]
        return subs_match

    def lipinski_rule(self):
        """
        Method to check if a compound passes the Lipinski filters in an sdf, returns True if so or else false.

        :return lipinski_filter: list of bool, return True if it passes a filtering criteria.
        """
        # Calculating the descriptors which will be required for the Lipinski rule calculation
        self.calc_descriptors()
        # Lipinski filter from the dataframe which returns a boolean values for every entry
        lipinski_filter = ((self.sdf_df["Molwt"] < 500) & (self.sdf_df["MolLogP"] < 5) & (self.sdf_df["HBD"] < 5) &
                           (self.sdf_df["HBA"] < 10))
        return lipinski_filter

    def sorting_problems(self, therapeutic_class):
        """
        Method to address the sorting problems listed in the code challenge.
        a. sort compounds by launched date
        b. sort compounds by therapeutic class, and within each therapeutic class by molecular weight

        :param therapeutic_class: str, the therapeutic class to be used for sorting.
        :return filtered_df: DataFrame, sorted table by Molwt for a particular therapeutic class.
        """
        # Calculating the descriptors which will be required for the sorting operation
        self.calc_descriptors()
        # Sort by launched date
        self.sdf_df.sort_values(by="Launched date", ascending=True)
        # Sort by therapeutic class and secondary sort by Molecular weight
        self.sdf_df.sort_values(["therapeutic class name", "Molwt"], ascending=[True, True])
        # Filter to select only rows that containing a therapeutic class name and sorting the filtered rows by Molwt.
        filt = self.sdf_df.loc[self.sdf_df['therapeutic class name'].str.contains(therapeutic_class, case=False)]
        filtered_df = filt.sort_values(by="Molwt", ascending=True)
        return filtered_df

    def stereochemistry_checker(self):
        """
        Method for finding the stereochemistry (absolute configurations) of molecules in an sdf.

        :return dict_smi_stereo_status: dict, returns a dictionary containing smiles and stereoinformation.
        :return stereo_count: dict, returns a dictionary of stereoactive and non-stereoactive compounds,
        """
        stereo_active = 0
        stereo_inactive = 0
        dict_smi_stereo_status = {}
        # Iterating over the smiles to using the RDKit FindMolChiralCenters method to calculate absolute stereo-
        # configuration.
        for smi in self.sdf_df["SMILES"]:
            chirality = Chem.FindMolChiralCenters(Chem.MolFromSmiles(smi), includeUnassigned=True)
            if not chirality:
                stereo_status = "Not Stereoactive"
                absolute_stereo_config = "N/A"
                stereo_inactive = stereo_inactive + 1
            else:
                stereo_status = "Stereoactive"
                absolute_stereo_config = chirality
                stereo_active = stereo_active + 1
            # Storing the results in a dictionary with the smiles as a key and the values as a tuple of stereostatus
            # and absolute stereoconfig.
            dict_smi_stereo_status[smi] = (stereo_status, absolute_stereo_config)
        # Storing the count of molecules that are stereoactive and stereo-inactive in a dictionary.
        stereo_count = {"StereoActive": stereo_active, "Stereo_Inactive": stereo_inactive}
        return dict_smi_stereo_status, stereo_count


def main(argv=sys.argv[1:]):

    args = get_arguments(argv)
    original_sdf = args.sdf

    get_sci_data_obj = ScientificData(original_sdf)

    # Checking for the arguments passed and calling only the required functions
    if args.calc_descriptors:
        get_sci_data_obj.calc_descriptors()
        print(get_sci_data_obj.sdf_df)

    if args.image:
        print(get_sci_data_obj.mol_image())

    if args.iupac_name:
        print(get_sci_data_obj.iupac_name_generator())

    if args.substructure_search:
        print(get_sci_data_obj.substructure_search(args.substructure_search))

    if args.lipinski_rule:
        print(get_sci_data_obj.lipinski_rule())

    if args.sorted_sdf:
        print(get_sci_data_obj.sorting_problems(args.sorted_sdf))

    if args.stereochemistry_checker:
        print(get_sci_data_obj.stereochemistry_checker())


if __name__ == '__main__':
    main()
