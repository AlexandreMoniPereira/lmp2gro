import re
import argparse
import os


def extract_numbers(text):
    # Extract all numbers including decimals and negatives
    numbers = re.findall(r'-?\d+\.?\d*', text)
    return [float(num) if '.' in num else int(num) for num in numbers]

def atom_id_to_type(atom_id, atom_data):
    atom_nr_2_type = dict(zip(atom_data['atom_id'], atom_data['gro_type']))
    return atom_nr_2_type.get(atom_id)

def check_duplicate_bond_types(bond_data):
    bond_type_count = bond_data[['bond_type', 'element_i', 'element_j']].drop_duplicates()
    pair_bond_type_count = bond_data[['element_i', 'element_j']].drop_duplicates()

    return len(bond_type_count)!=len(pair_bond_type_count)

def setup_argparser():

    """Set up the argument parser for command-line options."""

    parser = argparse.ArgumentParser(description="Convert molecular data files into GROMACS format.",
                                     formatter_class=argparse.RawDescriptionHelpFormatter,
                                    epilog='''
    Example usage:
        python lmp2gro.py datafile.data
        python lmp2gro.py datafile.data -r RESIDUE_NAME
    ''')

    parser.add_argument("datafile", type=str, help="Path to the input data file.")

    parser.add_argument('-r', '--resname', type=str, default='UNL',
                        help="Residue name to use in the output files (default: UNL).")
    
    return parser
