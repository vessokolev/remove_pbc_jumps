
bonds=[\
'        H2      O2',\
'        O2      C20',\
'        C20     H6',\
'        C20     H7',\
'        C20     C19',\
'        C19     H5',\
'        C19     O3',\
'        O3      H1',\
'        C19     C18',\
'        C18     H3',\
'        C18     H4',\
'        C18     O1',\
'        O1      C17',\
'        C17     O4',\
'        C17     C16',\
'        C16     H9',\
'        C16     H8',\
'        C16     C15',\
'        C15     H10',\
'        C15     H11',\
'        C15     C14',\
'        C14     H13',\
'        C14     H12',\
'        C14     C13',\
'        C13     H14',\
'        C13     H15',\
'        C13     C12',\
'        C12     H17',\
'        C12     H16',\
'        C12     C11',\
'        C11     H18',\
'        C11     H19',\
'        C11     C10',\
'        C10     H20',\
'        C10     H21',\
'        C10     C9',\
'        C9      H22',\
'        C9      C8',\
'        C8      H23',\
'        C8      C7',\
'        C7      H24',\
'        C7      H25',\
'        C7      C6',\
'        C6      H27',\
'        C6      H26',\
'        C6      C5',\
'        C5      H28',\
'        C5      H29',\
'        C5      C4',\
'        C4      H31',\
'        C4      H30',\
'        C4      C3',\
'        C3      H32',\
'        C3      H33',\
'        C3      C2',\
'        C2      H35',\
'        C2      H34',\
'        C2      C1',\
'        C1      H36',\
'        C1      H37',\
'        C1      C21',\
'        C21     H39',\
'        C21     H40',\
'        C21     H38'\
]


def get_atom_description(current_e_resnum,current_resnum,line):

   # Supply the 'line' without the "chomp" symbol at its end!!!
   # That means: line.rsplit()!

   flag=True

   try:
      resnum=int(line[0:5])
   except:
      flag=False

   if flag:
      if current_e_resnum==-1:
         current_resnum=resnum
         current_e_resnum=0
      else:
         if current_resnum!=resnum:
            current_resnum=resnum
            current_e_resnum+=1

      resnum=int(line[0:5])
      resname='{:<5}'.format(line[5:10].strip(' '))
      aname='{:<5}'.format(line[10:15].strip(' '))

      anum=int(line[15:20])

      coord=[float(line[20:28]),float(line[28:36]),float(line[36:44])]

      try:
         anum=int(line[15:20])
      except:
         flag=False

      if flag:
         return_list=[current_e_resnum,current_resnum,resname,aname,anum,coord]

   if not flag:
      return None
   else:
      return return_list


def read_gro_file(gro_in_file):

   counter=0

   to_read_comment=True
   to_read_numatom=True

   flag_read_num_line=False
   flag_read_comment=False

   comment=''

   # Initialize:
   current_e_resnum=-1
   current_resnum=-1
   res=[-1,'','']

   atoms=[]

   # Open the file without loading all its content into the memory:

   with open(gro_in_file,'r') as file_obj:

       for line in file_obj:
          if to_read_comment and counter==0:
             to_read_comment=False
             comment=line.rstrip()
          else:
             if to_read_numatom and counter==1:
                if to_read_comment:
                   print('FATAL ERROR! Badly formatted GRO file!')
                   sys.exit(1)
                else:
                   try:
                      num_atoms=int(line.rstrip())
                   except:
                      print('FATAL ERROR! Badly formatted GRO file!')
                      sys.exit(1)
             else:
                res=get_atom_description(res[0],res[1],line.rstrip())
                if counter==num_atoms+2:
                   tmp=line.split()
                   box=[float(tmp[0]),float(tmp[1]),float(tmp[2])]
                else:
                   if res is not None:
                      atoms.append(res)
          counter+=1

   return comment,num_atoms,atoms,box


def split_molecule(atoms,box):

   coords=[]

   for i in atoms:

      coord=[0.0,0.0,0.0]

      for j,k,l in zip([0,1,2],box,[0,1,2]):
         if i[5][j]>k:
            coord[l]=i[5][j]-k
         else:
            coord[l]=i[5][j]

      coords.append(coord)

   return coords


def write_gro(gro_out,comment,num_atoms,atoms,box):

   # Very simple Gromos87 GRO writer.

   file_obj=open(gro_out,'w')

   file_obj.write(comment+'\n')
   file_obj.write(str(num_atoms)+'\n')

   for i in atoms:
      line='%5d' % i[1]
      line+='%5s' % i[2]
      line+='%-5s' % i[3]
      line+='%5d' % i[4]
      line+='%8.3f' % i[5][0]
      line+='%8.3f' % i[5][1]
      line+='%8.3f' % i[5][2]

      file_obj.write(line+'\n')

   box_record= '%10.5f' % box[0]
   box_record+='%10.5f' % box[1]
   box_record+='%10.5f' % box[2]

   file_obj.write(box_record+'\n')

   file_obj.close()

   return None


def split_molecule_pbc(gro_in,gro_out):

   comment,num_atoms,atoms,box=read_gro_file(gro_in)

   half_box=[i/2 for i in box]

   coords=split_molecule(atoms,box)

   # Write the GRO file:

   write_gro(gro_out,comment,num_atoms,atoms,coords,box)

   return None


def coord_diff(atoms,i,j):

   return [atoms[i][5][k]-atoms[j][5][k] for k in [0,1,2]]


def gen_rand_sequence(lower_limit,upper_limit):

   import numpy

   num_members=upper_limit-lower_limit+1

   # Shift.
   upper_limit+=1

   sequence=[]

   flag=True

   while flag:
      # Generate a random integer in [lower_limit,upper_limit]
      rnd=lower_limit+numpy.random.randint(upper_limit)
      if not rnd in sequence:
         sequence.append(rnd)

      if len(sequence)==num_members:
         flag=False
      else:
         flag=True


   return sequence

