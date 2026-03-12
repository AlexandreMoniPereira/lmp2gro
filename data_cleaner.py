import pandas as pd
import numpy as np

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