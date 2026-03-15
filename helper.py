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
        python lmp2gro.py datafile.data -r RESIDUE_NAME --folder output_folder
        python lmp2gro.py datafile.data --clean -b "1 2" -a "1 2 3" -d "1" -i "1 2"
        python lmp2gro.py datafile.data --clean (removes all types with zero coefficients)
    ''')

    parser.add_argument('datafile', type=str, help="Path to the input data file.")

    parser.add_argument('-r', '--resname', type=str, default='UNL',
                        help="Residue name to use in the output files (default: UNL).")
    
    parser.add_argument('--clean', action='store_true',
                        help="Clean data by removing specified types in (-b for bonds, -a for angles, -d for dihedrals, and -i for impropers). If no specific types are provided, it will remove all types with zero coefficients.")

    parser.add_argument('-b', '--bond_rm', type=str, default='',
                        help="Requires --clean. String of bond type indices to remove (e.g., '1 2 3'). If empty, it will remove all bond types with zero coefficients.")

    parser.add_argument('-a', '--angle_rm', type=str, default='',
                        help="Requires --clean. String of angle type indices to remove (e.g., '1 2 3'). If empty, it will remove all angle types with zero coefficients.")
    
    parser.add_argument('-d', '--dihedral_rm', type=str, default='',
                        help="Requires --clean. String of dihedral type indices to remove (e.g., '1 2 3'). If empty, it will remove all dihedral types with zero coefficients.")
    
    parser.add_argument('-i', '--improper_rm', type=str, default='',
                        help="Requires --clean. String of improper type indices to remove (e.g., '1 2 3'). If empty, it will remove all improper types with zero coefficients.")                        

    parser.add_argument('--folder', type=str, default=None,
                        help="Optional name for the output folder.") 
    
    return parser
