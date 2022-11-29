#!/usr/bin/env python3
# coding:utf-8
from fireworks import LaunchPad
from pymatgen import Structure
import numpy as np
import mspglib as mspg
import spglib as spg
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
#from symmetrize import symmetrize
# create the atomate db from your db.json
#PATH_TO_MY_DB_JSON = "/nfs/home/jinsh/apps/atomate/config/db.json"
#atomate_db = VaspCalcDb.from_db_file(PATH_TO_MY_DB_JSON)
query = {"state":"COMPLETED"}
sym_prec = 0.01
lpad = LaunchPad.auto_load()
wf_ids = lpad.get_wf_ids(query=query)
for wf_id in wf_ids:
    uuid = None
    wf = lpad.get_wf_by_fw_id(wf_id)
    for fw in wf.fws:
        if "Magnetic Orderings Analysis" in fw.name:
            uuid = fw.tasks[0]["wf_uuid"]
    if uuid:
        analysis_entry = lpad.db.magnetic_orderings.find_one({'wf_meta.wf_uuid':{"$regex":uuid}})
        id_stable = analysis_entry["decomposes_to"]
        # id_stable may point to null means no structure is stable
        if id_stable:
            analysis_entry = lpad.db.magnetic_orderings.find_one({'task_id':id_stable})
        magnetism = analysis_entry["ordering"]
        print(f'{analysis_entry["formula"]}:{magnetism}')
        if magnetism != 'NM':
            uuid = analysis_entry["wf_meta"]["wf_uuid"]
            print(uuid) 
            struct = Structure.from_dict(analysis_entry["structure"])
            struct.to('mcif','temp.mcif')
            num = mspg.magnetic_spacegroup_number(struct,symprec=sym_prec)
            print(num)
            struct_s = mspg.symmetrize(struct,symprec=sym_prec)
            num_s = mspg.magnetic_spacegroup_number(struct_s,symprec=1e-4)
            print(num_s)
            if num_s != num:
                struct.to('mcif',f'./structs/{analysis_entry["formula_pretty"]}{uuid[-3:]}.mcif') 
