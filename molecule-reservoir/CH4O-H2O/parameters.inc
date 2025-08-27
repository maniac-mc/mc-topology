# LAMMPS Lennard-Jones pair coefficients
# epsilon (kcal/mol), sigma (Angstrom)

set type 1 charge 0.433743478261
set type 2 charge -0.192556521739
set type 3 charge -0.602556521739
set type 4 charge 0.158343478261
set type 5 charge 0.157043478261
set type 6 charge -0.300956521739
set type 7 charge 0.736043478261
set type 8 charge +0.4998
set type 9 charge +0.1546
set type 10 charge +0.0000
set type 11 charge -0.6544
set type 12 charge 0.000
set type 13 charge 0.527
set type 14 charge -1.054

pair_coeff 1 1 0.0860 3.400 # C1-C1
pair_coeff 2 2 0.0860 3.400 # C2-C2
pair_coeff 3 3 0.1094 3.400 # C3-C3
pair_coeff 4 4 0.0150 2.511 # H2-H2
pair_coeff 5 5 0.0157 2.650 # H3-H3
pair_coeff 6 6 0.1700 3.250 # N-N
pair_coeff 7 7 0.0125 1.960 # Zn-Zn
pair_coeff 8 8 0.2195 3.6499 # CH CH
pair_coeff 9 9 0.1791 3.1655 # OH OH
pair_coeff 10 10 0.000 1.000 # HO HO
pair_coeff 11 11 0.2108 3.6370 # Me Me
pair_coeff 12 12 0.1848 3.165 # OW-OW
pair_coeff 13 13 0.0000 1.000 # HW HW
pair_coeff 14 14 0.0000 1.000 # MW-MW

bond_coeff 1 346.543 1.490 # C1-C3
bond_coeff 2 488.000 1.335 # C1-N
bond_coeff 3 540.249 1.350 # C2-C2
bond_coeff 4 367.000 1.080 # C2-H2
bond_coeff 5 440.2210 1.370 # C2-N
bond_coeff 6 340.000 1.090 # C3-H3
bond_coeff 7 78.500 2.011 # N-Zn
bond_coeff 8 0.000 0.945 # OH-HO
bond_coeff 9 0.000 1.430 # OH-CH
bond_coeff 10 0.000 0.037 # OH-Me
bond_coeff 11 0.0000 0.1050 # OW-MW
bond_coeff 12 0.0000 0.9572 # OW-HW

angle_coeff 1 48.088 109.32 # C1-C3-H3
angle_coeff 2 71.254 105.27 # C1-N-C2
angle_coeff 3 48.680 128.33 # C1-N-Zn
angle_coeff 4 49.451 125.67 # C2-C2-H2
angle_coeff 5 76.750 108.65 # C2-C2-N
angle_coeff 6 32.477 126.40 # C2-N-Zn
angle_coeff 7 66.015 123.92 # C3-C1-N
angle_coeff 8 49.954 125.68 # H2-C2-N
angle_coeff 9 35.000 109.50 # H3-C3-H3
angle_coeff 10 75.484 112.16 # N-C1-N
angle_coeff 11 35.240 109.48 # N-Zn-N
angle_coeff 12 0.000 108.5 # CH-OH-HO
angle_coeff 13 0.000 127.53 # Me-OH-CH
angle_coeff 14 0.000 123.97 # Me-OH-HO
angle_coeff 15 0.0000 52.2600 # MW-OW-HW
angle_coeff 16 0.0000 104.5200 # HW-OW-HW

dihedral_coeff * 0.0 0 0 0.0
dihedral_coeff 1 2.325 2 180 0.0 # C2-C2-N-C1
dihedral_coeff 2 2.325 2 180 0.0 # C2-C2-N-Zn
dihedral_coeff 3 2.325 2 180 0.0 # H2-C2-N-Zn
dihedral_coeff 4 2.325 2 180 0.0 # H2-C2-N-C1
dihedral_coeff 5 5.000 2 180 0.0 # C3-C1-N-C2
dihedral_coeff 6 5.000 2 180 0.0 # C3-C1-N-Zn
dihedral_coeff 7 5.000 2 180 0.0 # N-C1-N-C2
dihedral_coeff 8 5.000 2 180 0.0 # N-C1-N-Zn
dihedral_coeff 9 5.150 2 180 0.0 # H2-C2-C2-H2
dihedral_coeff 10 5.150 2 180 0.0 # H2-C2-C2-N
dihedral_coeff 11 5.150 2 180 0.0 # N-C2-C2-N

improper_coeff * 0 0
improper_coeff 1 2.0 180.0 # C1-N-C2-Zn
improper_coeff 2 2.0 180.0 # C2-C2-H2-N
improper_coeff 3 2.0 180.0 # C3-C1-N-N
