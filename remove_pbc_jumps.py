#!/usr/bin/env python3

# This Python code assembles a molecule which was virtually
# split by imposing PBC. That process of "re-solving" the
# molecule topology is supported by a list of the chemical
# bonds defined between the atoms of the molecule. Note that
# the process of resolving the molecule topology is not
# possible to perform without having the bond list and the
# sizes of the simulation box. A bond is considered split
# if its size is equal or great then the half of the box
# size in the respective direction (x, y, or z).
#
# Author  : Veselin Kolev <vesso.kolev@gmail.com>
# Version : 2018013100
# License : GPLv2
# Requires: Python 3,NumPy

import helper
import sys

# The file 'gro_in' contains the split molecule.

gro_in ='split.gro'
gro_out='assembled.gro'

comment,num_atoms,atoms,box=helper.read_gro_file(gro_in)

# Compute the half of the box in all directions:

half_box=[i/2 for i in box]

anums =[i[4] for i in atoms]
anames=[i[3].split()[0] for i in atoms]

bonds_=[]

for bond in helper.bonds:
   bonds_.append([anums[anames.index(bond.split()[0])],\
                  anums[anames.index(bond.split()[1])]])

first_check=True

indexes=helper.gen_rand_sequence(0,len(bonds_)-1)

bonds=[bonds_[i] for i in indexes]

moved=[False for i in range(len(atoms))]

flag=True

while flag:

   flag=False

   for i in bonds:

      indx_1=anums.index(i[0])
      indx_2=anums.index(i[1])

      indx=None

      if first_check:
         first_check=False
         moved[indx_1]=True

      if not (moved[indx_1] and moved[indx_2]):

         if moved[indx_1]:

            indx=[indx_1,indx_2,indx_2,indx_2,indx_2,indx_2,\
                  indx_1,indx_2,indx_1,indx_1,indx_2]

         if moved[indx_2]:

            indx=[indx_2,indx_1,indx_1,indx_1,indx_1,indx_1,\
                  indx_2,indx_1,indx_2,indx_2,indx_1]

         if indx is not None:

            # Get the distance between the atoms and check if it is
            # greater than the half of the simulation box size in
            # each direction:

            diff=helper.coord_diff(atoms,indx[0],indx[1])

            check=[abs(diff[j])>half_box[j] for j in [0,1,2]]

            if sum(check)>0:
               for j in range(3):
                  if check[j]:
                     if diff[j]>0.0:
                        atoms[indx[2]][5][j]=atoms[indx[3]][5][j]+box[j]
                        flag=True
                     else:
                        atoms[indx[4]][5][j]=atoms[indx[5]][5][j]-box[j]
                        flag=True
            else:
               if moved[indx[6]]:
                  moved[indx[7]]=moved[indx[8]]
                  flag=True
               else:
                  moved[indx[9]]=moved[indx[10]]
                  flag=True


# Write out the GRO file containing the assembled molecule.

helper.write_gro(gro_out,comment,num_atoms,atoms,box)

sys.exit(0)
