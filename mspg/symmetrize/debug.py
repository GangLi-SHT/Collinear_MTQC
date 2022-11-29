#!/usr/bin/env python3
# coding:utf-8
from fireworks import LaunchPad
from pymatgen import Structure
import numpy as np
import mspglib as mspg
import spglib as spg
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from symmetrize import symmetrize
# create the atomate db from your db.json
#PATH_TO_MY_DB_JSON = "/nfs/home/jinsh/apps/atomate/config/db.json"
#atomate_db = VaspCalcDb.from_db_file(PATH_TO_MY_DB_JSON)
sym_prec = 0.01
filename = "./structs/Mn2As98d.mcif"
#filename = "temp.mcif"
struct = Structure.from_file(filename)
num = mspg.magnetic_spacegroup_number(struct,symprec=sym_prec)
struct_s = symmetrize(struct,symprec=sym_prec)
num_s = mspg.magnetic_spacegroup_number(struct_s,symprec=1e-2)
num_sg = SpacegroupAnalyzer(struct_s,1e-8).get_space_group_number()
num_sg1 = SpacegroupAnalyzer(struct,1e-2).get_space_group_number()
print(num,num_s)
print(num_sg1,num_sg)
struct1 = mspg.standard_primitive(struct,0.01)
#print(struct1)
#print(struct_s)
symmetries = mspg.magnetic_symmetries(struct1,0.01)
for ele in symmetries:
    print(ele)
struct1.to('mcif','temp_p.mcif')
print(struct1)
struct_s.to('mcif','temp_ps.mcif')
print(struct_s)
#ind = np.invert(np.isclose(struct_s.frac_coords,struct1.frac_coords,atol=1e-3))
