import math
import numpy as np
import pandas as pd
import helper
import fit_equations as fe

def atom_data(input_file,headers_dict,atom_count,element_data,types_dict):

    lammps_data = np.loadtxt(input_file,
                         skiprows=headers_dict['Atoms'] + 2,
                         max_rows=atom_count)

    column_names = ['atom_id', 'mol_id', 'atom_type', 'charge', 'x', 'y', 'z','vel_x','vel_y','vel_z']

    atom_data = pd.DataFrame(lammps_data, columns=column_names[:lammps_data.shape[1]])


    ff_types = generate_ff(headers_dict, open(input_file).readlines(), element_data,types_dict)
    element_mapping = dict(zip(ff_types['type_number'], ff_types['element']))
    type_mapping = dict(zip(ff_types['type_number'], ff_types['type']))

    atom_data['element'] = atom_data['atom_type'].map(element_mapping)
    atom_data['gro_type'] = atom_data['atom_type'].map(type_mapping)
    return atom_data,ff_types

def bond_data(input_file,headers_dict,bond_count,atom_data):

    lammps_data = np.loadtxt(input_file,
                         skiprows=headers_dict['Bonds'] + 2,
                         max_rows=bond_count)

    column_names = ['bond_id', 'bond_type', 'ai', 'aj']

    bond_data = pd.DataFrame(lammps_data, columns=column_names[:lammps_data.shape[1]])

    bond_data['element_i'] = bond_data['ai'].apply(lambda x: helper.atom_id_to_type(x,atom_data))
    bond_data['element_j'] = bond_data['aj'].apply(lambda x: helper.atom_id_to_type(x,atom_data))

    return bond_data

def bond_coeffs(input_file, bond_data, headers_dict, bond_types):
    lammps_data = np.loadtxt(input_file,
                             skiprows=headers_dict['Bond Coeffs'] + 2,
                             max_rows=bond_types)
    

    if lammps_data.ndim == 1:        # single row (1D)
        lammps_data = [lammps_data.tolist()]
    else:                            # already 2D
        lammps_data = lammps_data.tolist()

    column_names = ['bond_type', 'kb_lammps', 'dist_lammps']

    bond_coeffs_df = pd.DataFrame(lammps_data, columns=column_names[:len(lammps_data[0])])

    bond_types_coeffs = bond_data[['bond_type', 'element_i', 'element_j']].drop_duplicates() 
    bond_types_coeffs = bond_types_coeffs.sort_values(by='bond_type').reset_index(drop=True)

    bond_types_coeffs['kb_gro'] = bond_coeffs_df['kb_lammps'] * 4.184 * 100 * 2   
    bond_types_coeffs['dist_gro'] = bond_coeffs_df['dist_lammps'] / 10

    return bond_types_coeffs

def complete_bond_data(bond_data, bond_types_coeffs):
    bond_type_to_kgro = dict(zip(bond_types_coeffs['bond_type'], bond_types_coeffs['kb_gro']))
    bond_type_to_distgro = dict(zip(bond_types_coeffs['bond_type'], bond_types_coeffs['dist_gro']))

    bond_data['kb_gro'] = bond_data['bond_type'].map(bond_type_to_kgro)
    bond_data['dist_gro'] = bond_data['bond_type'].map(bond_type_to_distgro)

    return bond_data

def angle_data(input_file,headers_dict,angle_count,atom_data):
    lammps_data = np.loadtxt(input_file,
                         skiprows=headers_dict['Angles'] + 2,
                         max_rows=angle_count)

    column_names = ['angle_id', 'angle_type', 'ai', 'aj', 'ak']

    angle_data = pd.DataFrame(lammps_data, columns=column_names[:lammps_data.shape[1]])

    angle_data['element_i'] = angle_data['ai'].apply(lambda x: helper.atom_id_to_type(x,atom_data))
    angle_data['element_j'] = angle_data['aj'].apply(lambda x: helper.atom_id_to_type(x,atom_data))
    angle_data['element_k'] = angle_data['ak'].apply(lambda x: helper.atom_id_to_type(x,atom_data))

    return angle_data

def angle_coeffs(input_file, angle_data, headers_dict, angle_types):
    lammps_data = np.loadtxt(input_file,
                             skiprows=headers_dict['Angle Coeffs'] + 2,
                             max_rows=angle_types)

    if len(lammps_data.T) == 4:
        angle_style='cossine-periodic'
        column_names=['angle_type','C_lammps','B_lammps','n_lammps']
        angle_coeffs=pd.DataFrame(lammps_data, columns=column_names[:lammps_data.shape[1]])

        angle_coeffs[['k_eff_lammps', 'theta0_deg']] = angle_coeffs.apply(
            lambda row: fe.get_lammps_cos_per_params(C=row['C_lammps'],
                                          B=row['B_lammps'],
                                          n=row['n_lammps']), axis=1, result_type="expand")
        

    elif len(lammps_data.T) == 5:
        angle_style='fourier'
        column_names=['angle_type','k_lammps','C0_lammps','C1_lammps','C2_lammps']
        angle_coeffs=pd.DataFrame(lammps_data, columns=column_names[:lammps_data.shape[1]])

        angle_coeffs[['k_eff_lammps', 'theta0_deg']] = angle_coeffs.apply(
            lambda row: fe.get_fourier_gromos_params(K=row['k_lammps'],
                                          C0=row['C0_lammps'],
                                          C1=row['C1_lammps'],
                                          C2=row['C2_lammps']), axis=1, result_type="expand")

    elif len(lammps_data.T) == 3:
        if lammps_data.ndim == 1:        # single row (1D)
            lammps_data = [lammps_data.tolist()]
        else:                            # already 2D
            lammps_data = lammps_data.tolist()
        angle_style = 'harmonic'
        column_names = ['angle_type', 'k_eff_lammps', 'theta0_deg']
        angle_coeffs=pd.DataFrame(lammps_data, columns=column_names)
         
    angle_coeffs['k_eff_gromacs']=angle_coeffs['k_eff_lammps']*4.184
    angle_types_coeffs = angle_data[['angle_type', 'element_i', 'element_j','element_k']].drop_duplicates()
    angle_types_coeffs = angle_types_coeffs.sort_values(by='angle_type').reset_index(drop=True)

    return angle_coeffs,angle_types_coeffs,angle_style

def dihedral_improper_data(input_file, header_position, dihedral_count, atom_data):
    lammps_data = np.loadtxt(input_file,
                         skiprows=header_position + 2,
                         max_rows= dihedral_count)

    column_names = ['dihedral_id', 'dihedral_type', 'ai', 'aj', 'ak','al']

    dihedral_data = pd.DataFrame(lammps_data, columns=column_names[:lammps_data.shape[1]])

    dihedral_data['element_i'] = dihedral_data['ai'].apply(lambda x: helper.atom_id_to_type(x,atom_data))
    dihedral_data['element_j'] = dihedral_data['aj'].apply(lambda x: helper.atom_id_to_type(x,atom_data))
    dihedral_data['element_k'] = dihedral_data['ak'].apply(lambda x: helper.atom_id_to_type(x,atom_data))
    dihedral_data['element_l'] = dihedral_data['al'].apply(lambda x: helper.atom_id_to_type(x,atom_data))

    return  dihedral_data

def dihedral_coeffs(input_file, dihedral_data, headers_dict, dihedral_types):
    lammps_data = np.loadtxt(input_file,
                             skiprows=headers_dict['Dihedral Coeffs'] + 2,
                             max_rows=dihedral_types)

    column_names=['dihedral_type','K_lammps','d_lammps','n_lammps']
    dihedral_coeffs=pd.DataFrame(lammps_data, columns=column_names[:lammps_data.shape[1]])

    dihedral_coeffs['K_gromacs']=dihedral_coeffs['K_lammps']*4.184
    dihedral_coeffs['phi_s_gromacs'] = np.where(dihedral_coeffs['d_lammps'] == 1, 0, 180)
    dihedral_coeffs['n_gromacs'] = dihedral_coeffs['n_lammps']

    
    dihedral_types_coeffs =  dihedral_data[['dihedral_type', 'element_i', 'element_j','element_k','element_l']].drop_duplicates()
    dihedral_types_coeffs =  dihedral_types_coeffs.sort_values(by='dihedral_type').reset_index(drop=True)

    dihedral_types_coeffs = dihedral_types_coeffs.merge(
    dihedral_coeffs[['dihedral_type', 'K_gromacs', 'phi_s_gromacs', 'n_gromacs']],
    on='dihedral_type',
    how='left'  # Use 'left' to preserve existing rows in dihedral_types_coeffs
    )

    return dihedral_types_coeffs

def improper_coeffs(input_file, improper_data, headers_dict, improper_types):
   
    lammps_data = np.loadtxt(input_file,
                             skiprows=headers_dict['Improper Coeffs'] + 2,
                             max_rows=improper_types)

    if lammps_data.ndim == 1:        # single row (1D)
        lammps_data = [lammps_data.tolist()]

    column_names=['improper_type','k_lammps','C0_lammps','C1_lammps','C2_lammps','all']
    try:
        improper_coeffs=pd.DataFrame(lammps_data, columns=column_names[:lammps_data.shape[1]])
    except:
        improper_coeffs=pd.DataFrame(lammps_data, columns=column_names)

    improper_coeffs[['k_eff_lammps', 'theta0_deg']] =improper_coeffs.apply(
    lambda row: fe.get_fourier_harmonic_params(K_fourier=row['k_lammps'],
                                      C0=row['C0_lammps'],
                                      C1=row['C1_lammps'],
                                      C2=row['C2_lammps']), axis=1, result_type="expand")
    
    improper_coeffs['k_eff_gromacs']=improper_coeffs['k_eff_lammps']*4.184

    improper_types_coeffs =  improper_data[['dihedral_type', 'element_i', 'element_j','element_k','element_l']].drop_duplicates()
    improper_types_coeffs =  improper_types_coeffs.rename(columns={'dihedral_type':'improper_type'})
    improper_types_coeffs =  improper_types_coeffs.sort_values(by='improper_type').reset_index(drop=True)

    improper_types_coeffs = improper_types_coeffs.merge(
    improper_coeffs[['improper_type', 'k_eff_gromacs', 'theta0_deg']],
    on='improper_type',
    how='left'  # Use 'left' to preserve existing rows in improper_types_coeffs
    )
    
    return improper_types_coeffs

def extract_box_params(lines):
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

    gro_box_str=f'    {lx/10:.5f} {ly/10:.5f} {lz/10:.5f} {0:.5f} {0:.5f} {xy/10:.5f} {0/10:.5f} {xz/10:.5f} {yz/10:.5f}'

    return gro_box_str

#=============================== AUXILIARY FUNCTIONS ===============================

#================ FUNCTION TO GENERATE FORCE FIELD TYPES DATAFRAME =================
def generate_ff(headers_dict, lines, element_data,types_nr):

    element_map = dict(zip(element_data['element'], element_data['Z']))
    ff_types = pd.DataFrame({
        'type_number': [0] * types_nr['atom_types'],
        'masses': [None] * types_nr['atom_types'],
        'element': [None] * types_nr['atom_types'],
        'type': [None] * types_nr['atom_types'],
        'gro_sigma': [0.0] * types_nr['atom_types'],
        'gro_epsilon': [0.0] * types_nr['atom_types']
    })

    start_line = headers_dict['Masses'] + 2  # +1 for header, +1 for empty line

    for i in range(types_nr['atom_types']):
        masses_data = lines[start_line + i].split()
        ff_types.loc[i, 'type_number'] = int(masses_data[0])
        ff_types.loc[i, 'masses'] = float(masses_data[1])
        ff_types.loc[i, 'element'] = find_closest_element(float(masses_data[1]),element_data)

        try:
            ff_types.loc[i, 'type'] = masses_data[3]
        except:
            ff_types.loc[i, 'type'] = ff_types.loc[i, 'element'].replace(' ','')+'_'+str(ff_types.loc[i, 'type_number'])

    ff_types['atomic_number'] = ff_types['element'].map(element_map)
    #========= STARTING PAIR COEFFS =================================

    start_line = headers_dict['Pair Coeffs'] + 2  #+1 for header, +1 for empty line

    for i in range(types_nr['atom_types']):
        pair_data = lines[start_line + i].split()


        ff_types.loc[i, 'gro_sigma'] = float(pair_data[2])/10        #A to nm
        ff_types.loc[i, 'gro_epsilon'] = float(pair_data[1])*4.184   #kcal/mol to kj/mol
    return ff_types

def find_closest_element(mass,element_data):
    element_masses = element_data['A'].astype(float).values
    element_names = element_data['element'].values
    idx = np.argmin(np.abs(element_masses - mass))
    return element_names[idx]