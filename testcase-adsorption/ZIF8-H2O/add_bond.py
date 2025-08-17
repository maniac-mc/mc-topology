#!/usr/bin/env python
# coding: utf-8

import numpy as np
import MDAnalysis as mda


def write_lammps_data(filename, data, symmetrize = False):
    coords = data['coordinates']
    species_names = data['species_names']
    charges = data['charges']
    box = data['box']
    atom_names = data['atom_names']

    # Map unique species names to atom type IDs (1-based)
    unique_species = list(dict.fromkeys(species_names))  # preserve order
    species_to_typeid = {name: i + 1 for i, name in enumerate(unique_species)}

u = mda.Universe("zif8-water.data")

filename = "zif8-water.data"
symmetrize = True

with open(filename, 'w') as f:
    f.write("LAMMPS data file generated add_bond.py\n\n")

    num_atoms = u.atoms.n_atoms
    num_atom_types = len(np.unique((u.atoms.types)))
    f.write(f"{num_atoms} atoms\n")
    f.write(f"{num_atom_types} atom types\n")
    f.write(f"9 bonds\n")
    f.write(f"2 bond types\n")
    f.write(f"9 angles\n")
    f.write(f"2 angle types\n\n")


    # Box bounds
    if symmetrize:
        xlo, xhi = -u.dimensions[0]/2, u.dimensions[0]/2
        ylo, yhi = -u.dimensions[1]/2, u.dimensions[1]/2
        zlo, zhi = -u.dimensions[2]/2, u.dimensions[2]/2
    else:
        xlo, xhi = 0.0, u.dimensions[0]
        ylo, yhi = 0.0, u.dimensions[1]
        zlo, zhi = 0.0, u.dimensions[2]

    f.write(f"{xlo:.6f} {xhi:.6f} xlo xhi\n")
    f.write(f"{ylo:.6f} {yhi:.6f} ylo yhi\n")
    f.write(f"{zlo:.6f} {zhi:.6f} zlo zhi\n\n")

    # Masses (use dummy mass = 1.0 unless you have real masses)
    f.write("Masses\n\n")
    f.write("1 15.999 # Ow\n")
    f.write("2 1e-05 # Mw\n")
    f.write("3 1.008 # Hw\n")
    f.write("4 65.38 # Zn\n")
    f.write("5 12.011 # C1\n")
    f.write("6 12.011 # C2\n")
    f.write("7 1.008 # H1\n")
    f.write("8 12.011 # C3\n")
    f.write("9 1.008 # H2\n")
    f.write("10 1.008 # H3\n")
    f.write("11 14.007 # N\n\n")

    # Atoms section (atom ID, molecule-ID, atom-type, charge, x, y, z)
    f.write("Atoms # full\n\n")
    for i, (pos, name, q) in enumerate(zip(u.atoms.positions, u.atoms.types, u.atoms.charges), 1):

        f.write(f"{i} 1 {name} {q:.6f} {pos[0]:.8f} {pos[1]:.8f} {pos[2]:.8f}\n")

    bonds = []
    cpt = 0

    xyz = u.atoms.positions.T
    Lx, Ly, Lz = u.dimensions[:3]
    for x1, y1, z1, t1, i1 in zip(xyz[0], xyz[1], xyz[2], u.atoms.types, u.atoms.ids[:30]):
        for x2, y2, z2, t2, i2 in zip(xyz[0], xyz[1], xyz[2], u.atoms.types, u.atoms.ids[:30]):

            # d = np.sqrt((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)
            dx = x1 - x2 - Lx * np.round((x1 - x2) / Lx)
            dy = y1 - y2 - Ly * np.round((y1 - y2) / Ly)
            dz = z1 - z2 - Lz * np.round((z1 - z2) / Lz)
            d = np.sqrt(dx**2 + dy**2 + dz**2)

            if ((t1 == "1") & (t2 == "2")) | ((t1 == "2") & (t2 == "1")):
                if (i1 < i2) & (d < 0.4):
                    bonds.append([cpt + 1, 1, np.int32(i1), np.int32(i2)]) # O-M
                    cpt += 1

            if ((t1 == "1") & (t2 == "3")) | ((t1 == "3") & (t2 == "1")):
                if (i1 < i2) & (d < 1.2):
                    bonds.append([cpt + 1, 2, np.int32(i1), np.int32(i2)]) # O-H
                    cpt += 1

    bonds = np.array(bonds)

    f.write("\n")
    f.write("Bonds\n\n")
    for a, b, c, d in bonds:

        f.write(str(a) + " " + str(b) + " " + str(c) + " " + str(d) + "\n")

    angles = []
    cpt = 0
    for _, _, a1, b1 in bonds:
        for _, _, a2, b2 in bonds:
            if (a1 == a2):
                if b1 < b2:
                    if (u.atoms.types[u.atoms.ids == b1] == "2") | (u.atoms.types[u.atoms.ids == b2] == "2"):
                        angles.append([cpt+1, 1, b1, a1, b2]) # M-O-M
                        cpt += 1
                    else:
                        angles.append([cpt+1, 2, b1, a1, b2]) # H-O-H
                        cpt += 1

    angles = np.array(angles)

    f.write("\n")
    f.write("Angles\n\n")

    for a, b, c, d, e in angles:

        f.write(str(a) + " " + str(b) + " " + str(c) + " " + str(d) + " " + str(e) + "\n")

    f.write("\n")

