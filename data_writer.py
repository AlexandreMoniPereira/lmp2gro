import helper

def write_bonded_info(folder_name, types_dict, 
                      bond_count, bond_types_coeffs,
                      angle_count, angle_types_coeffs, angle_coeffs, angle_style,
                      dihedral_count, dihedral_types_coeffs,
                      improper_count, improper_types_coeffs):
    
    """
    Writes the bonded interaction parameters to a GROMACS .itp file.    
    """

    ffbonded_name=f'{folder_name}/ffbonded.itp'
    
    with open(ffbonded_name, 'w') as file:
        file.write(f''';Created using lmp2gro

''')
        if bond_count !=0:
            file.write(f'''
[ bondtypes ]
;     i       j     func           b0                 kb
''')
        for i in range(types_dict['bond_types']):
            line = (f"{bond_types_coeffs['element_i'][i]:>8}"
                    f"{bond_types_coeffs['element_j'][i]:>8}"
                    f"     1 "
                    f"{bond_types_coeffs['dist_gro'][i]:>20.6f}"
                    f"{bond_types_coeffs['kb_gro'][i]:>20.6f}\n")
            file.write(line)


        if angle_count !=0:
            if angle_style=='fourier':
                ang_num=2
            else:
                ang_num=1

            file.write(f'''
[ angletypes ]
;      i       j       k  func                  th0                cth
''')
        for i in range(types_dict['angle_types']):
            line = (f"{angle_types_coeffs['element_i'][i]:>8}"
                    f"{angle_types_coeffs['element_j'][i]:>8}"
                    f"{angle_types_coeffs['element_k'][i]:>8}"
                    f"     {ang_num} "
                    f"{angle_coeffs['theta0_deg'][i]:>20.6f}"
                    f"{angle_coeffs['k_eff_gromacs'][i]:>20.6f}\n")
            file.write(line)

        if dihedral_count !=0:
            file.write(f'''
[ dihedraltypes ]
;      i       j       k       l  func    phi_s               K        n
''')
            for i in range(len(dihedral_types_coeffs)):
                line = (f"{dihedral_types_coeffs['element_i'][i]:>8}"
                        f"{dihedral_types_coeffs['element_j'][i]:>8}"
                        f"{dihedral_types_coeffs['element_k'][i]:>8}"
                        f"{dihedral_types_coeffs['element_l'][i]:>8}"
                        f"     9 "
                        f"{dihedral_types_coeffs['phi_s_gromacs'][i]:>8.1f}"
                        f"{dihedral_types_coeffs['K_gromacs'][i]:>20.6f}"
                        f"{dihedral_types_coeffs['n_gromacs'][i]:>5.0f}\n")
                file.write(line)
            
        if improper_count !=0:
            file.write(f'''
[ dihedraltypes ]
;      i       j       k       l  func    thetha              K    
''')
            for i in range(len(improper_types_coeffs)):
                line = (f"{improper_types_coeffs['element_i'][i]:>8}"
                        f"{improper_types_coeffs['element_j'][i]:>8}"
                        f"{improper_types_coeffs['element_k'][i]:>8}"
                        f"{improper_types_coeffs['element_l'][i]:>8}"
                        f"     2 "
                        f"{improper_types_coeffs['theta0_deg'][i]:>8.1f}"
                        f"{improper_types_coeffs['k_eff_gromacs'][i]:>15.10f}\n")
                file.write(line)


def write_molecule_itp(folder_name, resname, atom_count, atom_data,
                       bond_count, bond_data,
                       angle_count, angle_data, angle_style,
                       dihedral_count, dihedral_data,
                       improper_count, improper_data):
    
    """Writes the molecule information to a GROMACS .itp file."""
       
    molecule_itp_name=f'{folder_name}/conf.itp'


    with open(molecule_itp_name, 'w') as file:
        file.write(f''';Created using lmp2gro

[ moleculetype ]
; name       nrexcl
{resname}               3
''')
        file.write(f'''[ atoms ]
;  nr   type  resnr  residu  atom   cgnr    charge
''')
        
        for i in range(atom_count):
                line = (f"{int(atom_data['atom_id'][i]):>5}"
                        f"{atom_data['gro_type'][i]:>5}"
                        f"     1 "
                        f"{resname:<5}"
                        f"{atom_data['element'][i].replace(' ',''):>5}"
                        f"{int(atom_data['atom_id'][i]):<5}"
                        f"{int(atom_data['atom_id'][i]):>5}"
                        f"{atom_data['charge'][i]:>15.6f}\n")
                file.write(line)
        

        if bond_count>0:
            if helper.check_duplicate_bond_types(bond_data):
                file.write(f'''
[ bonds ]
;     ai      aj   func           b0                 kb 
''')
        
                for i in range(bond_count):
                    line = (f"{int(bond_data['ai'][i]):>8}"
                            f"{int(bond_data['aj'][i]):>8}"
                            f"     1 "
                            f"{bond_data['dist_gro'][i]:>20.6f}"
                            f"{bond_data['kb_gro'][i]:>20.6f}\n")
                    file.write(line)
            else:
                file.write(f'''
[ bonds ]
;     ai      aj   func 
''')
        
                for i in range(bond_count):
                    line = (f"{int(bond_data['ai'][i]):>8}"
                            f"{int(bond_data['aj'][i]):>8}"
                            f"     1 \n")
                    file.write(line)

        if angle_count>0:
            if angle_style=='fourier':
                ang_num=2
            else:
                ang_num=1
            file.write(f'''
[ angles ]
;     ai      aj      ak   func 
''')
        
            for i in range(angle_count):
                line = (f"{int(angle_data['ai'][i]):>8}"
                        f"{int(angle_data['aj'][i]):>8}"
                        f"{int(angle_data['ak'][i]):>8}"
                        f"     {ang_num} \n")
                file.write(line)

        if dihedral_count !=0:
            file.write(f'''
[ dihedrals ]
;     ai      aj      ak      al   func  
''')
            
            for i in range(dihedral_count):
                    line = (f"{int(dihedral_data['ai'][i]):>8}"
                            f"{int(dihedral_data['aj'][i]):>8}"
                            f"{int(dihedral_data['ak'][i]):>8}"
                            f"{int(dihedral_data['al'][i]):>8}"
                            f"     9 \n")
                    file.write(line)

        if improper_count !=0:        
            file.write(f'''
[ dihedrals ]
;     ai      aj      ak      al   func  
''')
            
            for i in range(improper_count):
                    line = (f"{int(improper_data['ai'][i]):>8}"
                            f"{int(improper_data['aj'][i]):>8}"
                            f"{int(improper_data['ak'][i]):>8}"
                            f"{int(improper_data['al'][i]):>8}"
                            f"     2 \n")
                    file.write(line)


def write_atomtypes(ff_types, atom_types, folder_name):

    '''Writes the atom types .itp file, contains LJ parameters.'''

    atomtypes_name=f'{folder_name}/atomtypes.itp'
    with open(atomtypes_name, 'w') as file:
        file.write(f''';Created using lmp2gro

[ atomtypes ]
;name     at.nr      mass  charge   ptype     sigma   epsilon
''')

        for i in range(atom_types):
            line = (f"{ff_types['type'][i]:<5}"  # Left-aligned in 5 spaces
                    f"{ff_types['atomic_number'][i]:>10}"  # Right-aligned in 10 spaces
                    f"{ff_types['masses'][i]:>10.4f}"  # Right-aligned, 10 spaces, 4 decimals
                    f"  0.0000     A  "
                    f"{ff_types['gro_sigma'][i]:>10.4f}"  # Right-aligned, 10 spaces, 4 decimals
                    f"{ff_types['gro_epsilon'][i]:>10.4f}\n")  # Right-aligned, 10 spaces, 4 decimals
            file.write(line)

def write_gro_file(folder_name, atom_count, atom_data, resname,gro_box_str):

    gro_name=f'{folder_name}/conf.gro'

    with open(gro_name, 'w') as file:
        file.write(f'''Created using lmp2gro\n''')
        file.write(f"{atom_count:>5}\n")

        for i in range(atom_count):
            line = (f"{int(atom_data['mol_id'][i]):>5}"
                    f"{resname:<5}"
                    f"{atom_data['element'][i]:>5}"
                    f"{int(atom_data['atom_id'][i]):>5}"
                    f"{atom_data['x'][i]/10:>8.3f}"
                    f"{atom_data['y'][i]/10:>8.3f}"
                    f"{atom_data['z'][i]/10:>8.3f}\n")
            file.write(line)

        file.write(gro_box_str)

def write_topology(folder_name,resname):


    topology_name=f'{folder_name}/topol.top'
    with open(topology_name, 'w') as file:
        file.write(f''';Created using lmp2gro
;
;	Example topology file
;
[ defaults ]
; nbfunc        comb-rule       gen-pairs       fudgeLJ fudgeQQ
  1             3               yes              1.0     0.0

; The force field files to be included
#include "atomtypes.itp"
#include "ffbonded.itp"
#include "conf.itp"

[ system ]
Example system

[ molecules ]
{resname}	1
''')