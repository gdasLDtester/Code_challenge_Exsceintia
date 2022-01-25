# Overview #

This directory contains code for the Code Challenge round of Scientific Data Engineer position at Exscientia.

#Pre-requisites
Use python3.6 or later and execute scientific_data_challenge.py
```bash
python --version
```
## Required Packages/Modules
Begin by installing the required packages listed in requirements.txt using the command:
```bash
pip3 install -r requirements.txt
```
## Alternative RDKit installation 
Alternatively RDKit can also be installed via Conda using the following steps
```bash
$ conda create -c conda-forge -n my-rdkit-env rdkit
$ conda activate my-rdkit-env
$ cd [anaconda folder]/bin
$ source activate my-rdkit-env
```
## Executing the script
Help information
```bash
python3 scientific_data_challenge.py --help
```
Actual running the script to calculate descriptors of an sdf
```bash
python3 scientific_data_challenge.py --sdf <sdf_file_path> --calc_descriptors
```
Running of the script to generate IUPAC names using PubChemPy
```bash
python3 scientific_data_challenge.py --sdf <sdf_file_path> --iupac_name
```
Simultaneously the script can be used to calculate substructure search, check for lipinski rule, print images, print a sorted dataframe or run a stereochemistry check on a set of molecules by, 
```bash
python3 scientific_data_challenge.py --sdf <sdf_file_path> --substructure_search
```
```bash
python3 scientific_data_challenge.py --sdf <sdf_file_path> --lipinski_rule
```
```bash
python3 scientific_data_challenge.py --sdf <sdf_file_path> --image
```
```bash
python3 scientific_data_challenge.py --sdf <sdf_file_path> --sorted_df
```
```bash
python3 scientific_data_challenge.py --sdf <sdf_file_path> --stereochemistry_checker
```

