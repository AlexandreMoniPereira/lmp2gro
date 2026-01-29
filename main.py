import os
import re
import pandas as pd
import numpy as np
import math


import helper
import data_colector as dc
import data_writer as dw


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

def gen_gro_box(lines):
    xy, xz, yz = 0,0,0
    for line in lines:
        if "xlo xhi" in line:
            xlo, xhi = helper.extract_numbers(line)
            lx = xhi - xlo
        if "ylo yhi" in line:
            ylo, yhi = helper.extract_numbers(line)
            ly = yhi - ylo
        if "zlo zhi" in line:
            zlo, zhi = helper.extract_numbers(line)
            lz = zhi - zlo
        if "xy xz yz" in line:
            xy, xz, yz = helper.extract_numbers(line)
            break

    a=lx
    b=(ly**2+xy**2)**(1/2)
    c=(lz**2+xz**2+yz**2)**(1/2)

    alpha = (math.acos(((xy*xz)+(ly*yz))/(b*c)))*57.2958
    beta = (math.acos(xz/c))*57.2958
    gamma = (math.acos(xy/b))*57.2958

    gromacs_box_line=f'    {lx/10:.5f} {ly/10:.5f} {lz/10:.5f} {0:.5f} {0:.5f} {xy/10:.5f} {0/10:.5f} {xz/10:.5f} {yz/10:.5f}'

    return gromacs_box_line


if __name__ == "__main__":

    #=================== Input name and creation of output folder  ===================
    config=helper.setup_argparser()
    args=config.parse_args()

    input_file=args.datafile
    resname=args.resname

    with open(input_file, 'r') as file:
        lines = file.readlines()

    output_folder = create_output_folder(input_file)

    #=================== Extracting General Information ===================
    comment=lines[0]
    atom_count=helper.extract_numbers(lines[2])[0]
    bond_count=helper.extract_numbers(lines[3])[0]
    angle_count=helper.extract_numbers(lines[4])[0]
    dihedral_count=helper.extract_numbers(lines[5])[0]
    improper_count=helper.extract_numbers(lines[6])[0]

    types_dict=gen_types_dict(lines)
    headers_dict=gen_headers_dict(lines)
    gro_box_line=gen_gro_box(lines)

    element_data=pd.read_csv( 'element_data.csv')

    #=================== Extracting Data ===================

    atom_data,ff_types=dc.atom_data(input_file,headers_dict,atom_count,element_data,types_dict)

    if bond_count>0:
        bond_data=dc.bond_data(input_file,headers_dict,bond_count,atom_data)
        bond_types_coeffs=dc.bond_coeffs(input_file, bond_data, headers_dict, types_dict['bond_types'])

        if helper.check_duplicate_bond_types(bond_data):
            bond_data=dc.complete_bond_data(bond_data,bond_types_coeffs)
    else:
        bond_data=None
        bond_types_coeffs=None

    if angle_count>0:
        angle_data=dc.angle_data(input_file,headers_dict,angle_count,atom_data)
        angle_coeffs,angle_types_coeffs,angle_style=dc.angle_coeffs(input_file, angle_data, headers_dict, types_dict['angle_types'])
    else:
        angle_data=None
        angle_coeffs=None
        angle_style=None

    if dihedral_count>0:
        dihedral_data=dc.dihedral_improper_data(input_file, headers_dict['Dihedrals'], dihedral_count, atom_data)
        dihedral_types_coeffs=dc.dihedral_coeffs(input_file, dihedral_data, headers_dict, types_dict['dihedral_types'])
    else:
        dihedral_data=None
        dihedral_types_coeffs=None

    if improper_count>0:
        improper_data=dc.dihedral_improper_data(input_file, headers_dict['Impropers'], improper_count, atom_data)
        improper_types_coeffs=dc.improper_coeffs(input_file, improper_data, headers_dict, types_dict['improper_types'])
    else:
        improper_data=None
        improper_types_coeffs=None

    #=================== Writing gromacs files ===================
    
    dw.write_bonded_info(output_folder, types_dict, 
                         bond_count, bond_types_coeffs,
                         angle_count, angle_types_coeffs, angle_coeffs, angle_style,
                         dihedral_count, dihedral_types_coeffs,
                         improper_count, improper_types_coeffs)
    
    dw.write_molecule_itp(output_folder, resname, atom_count, atom_data,
                       bond_count, bond_data,
                       angle_count, angle_data, angle_style,
                       dihedral_count, dihedral_data,
                       improper_count, improper_data)
    
    
    dw.write_atomtypes(ff_types, types_dict['atom_types'], output_folder)

    gro_box_str=dc.extract_box_params(lines)
    dw.write_gro_file(output_folder, atom_count, atom_data, resname,gro_box_str)

    dw.write_topology(output_folder, resname)