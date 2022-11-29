from pymatgen.core import Structure
import mspglib as mspg
struct = Structure.from_file("temp.mcif")
num = mspg.magnetic_spacegroup_number(struct)
print(num)
struct1 = Structure.from_file("temp1.mcif")
num = mspg.magnetic_spacegroup_number(struct1)
print(num)
