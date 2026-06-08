import os
import re
import pandas as pd
import numpy as np
import math


import helper
import data_colector as dc
import data_writer as dw
import data_cleaner as dcl


def create_output_folder(input_file):
    base_name = 'output_' + input_file.replace('data','').replace('.','')
    folder_name = base_name
    counter = 1
    
    while os.path.exists(folder_name):
        folder_name = f"output_{counter}_{input_file.replace('data','').replace('.','')}"
        counter += 1
    
    os.makedirs(folder_name)
    return folder_name

def gen_types_dict(lines):
    types_dict={}
    for each in lines: 
        if 'types' in each:
            line=each.split()
            variable=line[1]+'_'+line[2]
            value=line[0]
    
            types_dict[variable] = int(value) 
    return types_dict


def gen_headers_dict(lines, default=-1):
    headers=['Masses','Pair Coeffs','Bond Coeffs','Angle Coeffs','Dihedral Coeffs'
             ,'Improper Coeffs','Atoms','Bonds','Angles','Dihedrals','Impropers']
    headers_dict = {}

    for i, line in enumerate(lines):
        for header in headers:
            if header not in headers_dict and header in line:
                headers_dict[header] = i
                if len(headers_dict) == len(headers):
                    return headers_dict

    # Add default values for missing headers
    for header in headers:
        if header not in headers_dict:
            headers_dict[header] = default

    return headers_dict


if __name__ == "__main__":

    #=================== Input name and creation of output folder  ===================
    config=helper.setup_argparser()
    args=config.parse_args()

    input_file=args.datafile
    resname=args.resname
    clean=args.clean
    bond_rm_str=args.bond_rm
    angle_rm_str=args.angle_rm
    dihedral_rm_str=args.dihedral_rm
    improper_rm_str=args.improper_rm
    output_folder = args.folder

    with open(input_file, 'r') as file:
        lines = file.readlines()

    if output_folder is None:
        output_folder = create_output_folder(input_file)
    else:
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)
        else:
            print(f"Warning: Output folder '{output_folder}' already exists. Files may be overwritten.")

    #=================== Extracting General Information ===================
    comment=lines[0]

    count_titles={"atoms", "bonds", "angles", "dihedrals", "impropers"}
    counts = {key: 0 for key in count_titles}

    for i, line in enumerate(lines):
        stripped = line.strip()
        
        for key in count_titles:
            if stripped.endswith(key):
                numbers = helper.extract_numbers(line)
                if numbers:
                    counts[key] = numbers[0]
                break

    atom_count=counts["atoms"]
    bond_count=counts["bonds"]
    angle_count=counts["angles"]
    dihedral_count=counts["dihedrals"]
    improper_count=counts["impropers"]

    types_dict=gen_types_dict(lines)
    headers_dict=gen_headers_dict(lines)

    element_data=pd.read_csv( 'element_data.csv')

    #=================== Extracting Data ===================
    print(f"Extracting and converting atom data for {atom_count} atoms.")
    atom_data,ff_types=dc.atom_data(input_file,headers_dict,atom_count,element_data,types_dict)

    with helper.TqdmSpinner(f"Extracting and converting bond data for {bond_count} bonds."):
        if bond_count>0:
            bond_data=dc.bond_data(input_file,headers_dict,bond_count,atom_data)
            bond_types_coeffs=dc.bond_coeffs(input_file, bond_data, headers_dict, types_dict['bond_types'])

            if helper.check_duplicate_bond_types(bond_data):
                bond_data=dc.complete_bond_data(bond_data,bond_types_coeffs)

            if clean ==True:
                bond_types_coeffs, bond_coeffs, bond_data, bond_count, types_dict['bond_types'] = dcl.clean_bond_data(bond_types_coeffs, bond_types_coeffs, bond_data, bond_rm_str, bond_count, types_dict['bond_types'])
        else:
            bond_data=None
            bond_types_coeffs=None
    with helper.TqdmSpinner(f"Extracting and converting angle data for {angle_count} angles."):
        if angle_count>0:
            angle_data=dc.angle_data(input_file,headers_dict,angle_count,atom_data)
            angle_coeffs,angle_types_coeffs,angle_style=dc.angle_coeffs(input_file, angle_data, headers_dict, types_dict['angle_types'])

            if clean ==True:
                angle_types_coeffs,angle_coeffs, angle_data, angle_count, types_dict['angle_types'] = dcl.clean_angle_data(angle_types_coeffs, angle_coeffs, angle_data, angle_rm_str, angle_count, types_dict['angle_types'])

        else:
            angle_data=None
            angle_coeffs=None
            angle_style=None
    with helper.TqdmSpinner(f"Extracting and converting dihedral data for {dihedral_count} dihedrals."):
        if dihedral_count>0:
            dihedral_data=dc.dihedral_improper_data(input_file, headers_dict['Dihedrals'], dihedral_count, atom_data)
            dihedral_types_coeffs,dihedral_style=dc.dihedral_coeffs(input_file, dihedral_data, headers_dict, types_dict['dihedral_types'])
            if clean ==True:
                dihedral_types_coeffs, dihedral_data, dihedral_count, types_dict['dihedral_types'] = dcl.clean_dihedral_data(dihedral_types_coeffs, dihedral_data, dihedral_rm_str, dihedral_count, types_dict['dihedral_types'])
        else:
            dihedral_data=None
            dihedral_types_coeffs=None
            dihedral_style=None
    
    with helper.TqdmSpinner(f"Extracting and converting improper data for {improper_count} impropers."):
        if improper_count>0:
            improper_data=dc.dihedral_improper_data(input_file, headers_dict['Impropers'], improper_count, atom_data)
            improper_types_coeffs,improper_style=dc.improper_coeffs(input_file, improper_data, headers_dict, types_dict['improper_types'])
            if clean ==True:
                improper_types_coeffs, improper_data, improper_count, types_dict['improper_types'] = dcl.clean_improper_data(improper_types_coeffs,  improper_data, improper_rm_str, improper_count, types_dict['improper_types'])
        
        else:
            improper_data=None
            improper_types_coeffs=None
            improper_style=None

    #=================== Writing gromacs files ===================
    print(f"Writing GROMACS files on the {output_folder} folder.")
    dw.write_bonded_info(output_folder, types_dict, 
                         bond_count, bond_types_coeffs,
                         angle_count, angle_types_coeffs, angle_coeffs, angle_style,
                         dihedral_count, dihedral_types_coeffs, dihedral_style,
                         improper_count, improper_types_coeffs, improper_style)
    
    dw.write_molecule_itp(output_folder, resname, atom_count, atom_data,
                       bond_count, bond_data,
                       angle_count, angle_data, angle_style,
                       dihedral_count, dihedral_data, dihedral_style,
                       improper_count, improper_data,improper_style)
    
    
    dw.write_atomtypes(ff_types, types_dict['atom_types'], output_folder)

    gro_box_str=dc.extract_box_params(lines)
    dw.write_gro_file(output_folder, atom_count, atom_data, resname,gro_box_str)

    dw.write_topology(output_folder, resname)
