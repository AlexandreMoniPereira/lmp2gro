import pandas as pd
import numpy as np


def clean_bond_data(bond_types_coeffs, bond_coeffs, bond_data, bond_rm_str, bond_count, bond_types):
    
    if bond_rm_str != '':
        bond_rm_list = [int(x) for x in bond_rm_str.split()]
    else:
        # select rows where every column except 'bond_type' and 'element_*' columns is zero
        element_cols = [col for col in bond_coeffs.columns if 'element_' in col]
        mask = (bond_coeffs.drop(columns=['bond_type'] + element_cols) == 0).all(axis=1)
        bond_rm_list = bond_coeffs.loc[mask, 'bond_type'].tolist()

    # Remove rows from bond_types_coeffs based on bond_rm_list
    bond_types_coeffs_cleaned = bond_types_coeffs[~bond_types_coeffs['bond_type'].isin(bond_rm_list)].reset_index(drop=True)
    
    # Remove rows from bond_data based on bond_rm_list
    bond_data_cleaned = bond_data[~bond_data['bond_type'].isin(bond_rm_list)].reset_index(drop=True)

    # Remove rows from bond_coeffs based on bond_rm_list
    bond_coeffs_cleaned = bond_coeffs[~bond_coeffs['bond_type'].isin(bond_rm_list)].reset_index(drop=True)

    # Create a dictionary to map old bond_type values to new ones
    bond_type_mapping = {old: new for new, old in enumerate(bond_types_coeffs_cleaned['bond_type'], start=1)}
    #update bond_type values in bond_data_cleaned using the mapping
    bond_data_cleaned['bond_type'] = bond_data_cleaned['bond_type'].map(bond_type_mapping)

    # If the "cleaned" is now empty, set the dataframes to None
    if bond_types_coeffs_cleaned.empty:
        bond_types_coeffs_cleaned = None
    if bond_data_cleaned.empty:
        bond_data_cleaned = None
    if bond_coeffs_cleaned.empty:
        bond_coeffs_cleaned = None

    # update bond_counts and bond_types in the main function after cleaning
    bond_count = len(bond_data_cleaned) if bond_data_cleaned is not None else 0
    bond_types = len(bond_types_coeffs_cleaned) if bond_types_coeffs_cleaned is not None else 0
    
    return bond_types_coeffs_cleaned, bond_coeffs_cleaned, bond_data_cleaned, bond_count, bond_types

def clean_angle_data(angle_types_coeffs,angle_coeffs, angle_data, angle_rm_str, angle_count, angle_types):

    if angle_rm_str != '':
        angle_rm_list = [int(x) for x in angle_rm_str.split()]
    else:
        # select rows where every column except 'angle_type' and 'element_*' columns is zero
        element_cols = [col for col in angle_coeffs.columns if 'element_' in col]
        mask = (angle_coeffs.drop(columns=['angle_type'] + element_cols) == 0).all(axis=1)
        angle_rm_list = angle_coeffs.loc[mask, 'angle_type'].tolist()

    # Remove rows from angle_types_coeffs based on angle_rm_list
    angle_types_coeffs_cleaned = angle_types_coeffs[~angle_types_coeffs['angle_type'].isin(angle_rm_list)].reset_index(drop=True)
    
    # Remove rows from angle_data based on angle_rm_list
    angle_data_cleaned = angle_data[~angle_data['angle_type'].isin(angle_rm_list)].reset_index(drop=True)

    # Remove rows from angle_coeffs based on angle_rm_list
    angle_coeffs_cleaned = angle_coeffs[~angle_coeffs['angle_type'].isin(angle_rm_list)].reset_index(drop=True)

    # Create a dictionary to map old angle_type values to new ones
    angle_type_mapping = {old: new for new, old in enumerate(angle_types_coeffs_cleaned['angle_type'], start=1)}
    #update angle_type values in angle_data_cleaned using the mapping
    angle_data_cleaned['angle_type'] = angle_data_cleaned['angle_type'].map(angle_type_mapping)

    # If the "cleaned" is now empty, set the dataframes to None
    if angle_types_coeffs_cleaned.empty:
        angle_types_coeffs_cleaned = None
    if angle_data_cleaned.empty:
        angle_data_cleaned = None
    if angle_coeffs_cleaned.empty:
        angle_coeffs_cleaned = None

    # update angle_counts and angle_types in the main function after cleaning
    angle_count = len(angle_data_cleaned) if angle_data_cleaned is not None else 0
    angle_types = len(angle_types_coeffs_cleaned) if angle_types_coeffs_cleaned is not None else 0
    
    return angle_types_coeffs_cleaned,angle_coeffs_cleaned, angle_data_cleaned, angle_count, angle_types

def clean_dihedral_data(dihedral_types_coeffs, dihedral_data, dihedral_rm_str, dihedral_count, dihedral_types):
    
    if dihedral_rm_str != '':
        dihedral_rm_list = [int(x) for x in dihedral_rm_str.split()]
    else:
        # select rows where every column except 'dihedral_type' and 'element_*' columns is zero
        element_cols = [col for col in dihedral_types_coeffs.columns if 'element_' in col]
        mask = (dihedral_types_coeffs.drop(columns=['dihedral_type'] + element_cols) == 0).all(axis=1)
        dihedral_rm_list = dihedral_types_coeffs.loc[mask, 'dihedral_type'].tolist()

    # Remove rows from dihedral_types_coeffs based on dihedral_rm_list
    dihedral_types_coeffs_cleaned = dihedral_types_coeffs[~dihedral_types_coeffs['dihedral_type'].isin(dihedral_rm_list)].reset_index(drop=True)
    
    # Remove rows from dihedral_data based on dihedral_rm_list
    dihedral_data_cleaned = dihedral_data[~dihedral_data['dihedral_type'].isin(dihedral_rm_list)].reset_index(drop=True)

    # Create a dictionary to map old dihedral_type values to new ones
    dihedral_type_mapping = {old: new for new, old in enumerate(dihedral_types_coeffs_cleaned['dihedral_type'], start=1)}
    #update dihedral_type values in dihedral_data_cleaned using the mapping
    dihedral_data_cleaned['dihedral_type'] = dihedral_data_cleaned['dihedral_type'].map(dihedral_type_mapping)

    # If the "cleaned" is now empty, set the dataframes to None
    if dihedral_types_coeffs_cleaned.empty:
        dihedral_types_coeffs_cleaned = None
    if dihedral_data_cleaned.empty:
        dihedral_data_cleaned = None
    # update dihedral_counts and dihedral_types in the main function after cleaning
    dihedral_count = len(dihedral_data_cleaned) if dihedral_data_cleaned is not None else 0
    dihedral_types = len(dihedral_types_coeffs_cleaned) if dihedral_types_coeffs_cleaned is not None else 0 

    return dihedral_types_coeffs_cleaned, dihedral_data_cleaned, dihedral_count, dihedral_types

def clean_improper_data(improper_types_coeffs,  improper_data, improper_rm_str, improper_count, improper_types):
    
    if improper_rm_str != '':
        improper_rm_list = [int(x) for x in improper_rm_str.split()]
    else:
        # select rows where every column except 'improper_type' and 'element_*' columns is zero
        element_cols = [col for col in improper_types_coeffs.columns if 'element_' in col]
        mask = (improper_types_coeffs.drop(columns=['improper_type'] + element_cols) == 0).all(axis=1)
        improper_rm_list = improper_types_coeffs.loc[mask, 'improper_type'].tolist()

    # Remove rows from improper_types_coeffs based on improper_rm_list
    improper_types_coeffs_cleaned = improper_types_coeffs[~improper_types_coeffs['improper_type'].isin(improper_rm_list)].reset_index(drop=True)
    
    # Remove rows from improper_data based on improper_rm_list
    improper_data_cleaned = improper_data[~improper_data['dihedral_type'].isin(improper_rm_list)].reset_index(drop=True)

    # Create a dictionary to map old dihedral_type values to new ones
    improper_type_mapping = {old: new for new, old in enumerate(improper_types_coeffs_cleaned['improper_type'], start=1)}
    #update dihedral_type values in dihedral_data_cleaned using the mapping
    improper_data_cleaned['dihedral_type'] = improper_data_cleaned['dihedral_type'].map(improper_type_mapping)

    # If the "cleaned" is now empty, set the dataframes to None
    if improper_types_coeffs_cleaned.empty:
        improper_types_coeffs_cleaned = None
    if improper_data_cleaned.empty:
        improper_data_cleaned = None

    # update dihedral_counts and dihedral_types in the main function after cleaning
    improper_count = len(improper_data_cleaned) if improper_data_cleaned is not None else 0
    improper_types = len(improper_types_coeffs_cleaned) if improper_types_coeffs_cleaned is not None else 0 

    return improper_types_coeffs_cleaned, improper_data_cleaned, improper_count, improper_types