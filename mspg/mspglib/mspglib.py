from pymatgen.core import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.electronic_structure.core import Magmom
from pymatgen.core.operations import MagSymmOp
import numpy as np
import numpy.linalg as LA
import collections
from  . import _mspglib as mspg
import os
from ._permutation  import *

datafilepath = os.path.join(os.path.dirname(__file__), "Magnetic_Sym_El/")
data_filename = os.listdir(datafilepath)
data_filename = [ ele  for ele in data_filename if '._Magnetic_Sym_El' not in  ele ]

class SpglibError(object):
    message = "no error"

spglib_error = SpglibError()


def standard_primitive(struct,symprec=1e-5):
    Analyzer = SpacegroupAnalyzer(struct)
    dataset = get_symmetry_dataset(struct,symprec)
    
    lattice = dataset["primitive_lattice"]
    positions = dataset["positions"]
    types = dataset["primitive_types"]
    magmoms = dataset["magmoms"]
    _, T_mat,trans,_ = _compare_with_magnetic_database(struct,symprec)
    
    species = [Analyzer._unique_species[i - 1] for i in types]
    s = Structure(LA.inv(T_mat).T.dot(lattice), species, np.dot(positions,np.array(T_mat).T)+np.array(trans),\
                  site_properties = {"magmom":magmoms},to_unit_cell = True)
    return s

def magnetic_spacegroup_number(struct,symprec = 1e-5):
    num,_,_,_ = _compare_with_magnetic_database(struct,symprec)
    return num

def _compare_with_magnetic_database(struct,symprec=1e-5):    

    dataset = get_symmetry_dataset(struct,symprec)
    
    rotations = dataset["rotations"]
    translations = dataset["translations"]
    time_reversal = dataset["time_reversal"]
    
    Mtype = dataset['msgtype']
    SG_num = dataset['number']
    T_mat = np.eye(3)
    tran = [0.0,0.0,0.0]
    translations = np.array([ [ rational(t) for t in ele]  for ele in translations])
    magnetic_symmetries = [MagSymmOp.from_rotation_and_translation_and_time_reversal(rot,np.mod(trans,1),t) 
                           for rot,trans,t in zip(rotations,translations,time_reversal)]
    MSG_database = MSG(SG_num,Mtype)
    if Mtype >2 and SG_num in SG_ambi:
        magnetic_symmetries_c = magnetic_symmetries
        perm = permutation_orthorhombic[SG_ambi.index(SG_num)]
        shift = origin_shift_orthorhombic[SG_ambi.index(SG_num)]
        for i,p in enumerate([5]+perm):
            if p!=5:
                T_mat = np.array(permu_mat_orthorhombic[p]).T
                centering = dataset["international"][0]
                M = Mat2P[centering]
                T_mat = LA.inv(M).dot(T_mat).dot(M)
                tran = LA.inv(M).dot(shift[i-1])
                magnetic_symmetries_c = [ trans_magsymmop(ele,T_mat,tran) for ele in magnetic_symmetries]
            compare = [compare_magnetic_space_groups(magnetic_symmetries_c,ele,symprec) for ele in MSG_database.set_MSG]
            if np.any(compare) :
                ind = np.nonzero(compare)[0][0]
                break

    elif Mtype == 4 and SG_num<3:
        # determine lattice vector c
        # determine a and b
        c = 2*magnetic_symmetries[SG_num].translation_vector
        c = np.around(c,10).astype(int)
        T_mat = np.zeros((3,3)).astype(int)
        if np.allclose(c,[1,1,1]):
            T_mat = np.eye(3).astype(np.int)
            T_mat[2] = c
        else:
            ind = np.nonzero(c==0)[0]
            if len(ind)==2:
                T_mat[[0,1],ind] = 1
            if len(ind)==1:
                T_mat[0,ind[0]] = 1
                T_mat[1,np.mod(ind[0]+1,3)] = 1
            T_mat[2] = c
        T_mat = LA.inv(T_mat.T)
        magnetic_symmetries = [ trans_magsymmop(ele,T_mat) for ele in magnetic_symmetries]
        ind = np.nonzero([compare_magnetic_space_groups(magnetic_symmetries,ele,symprec) for ele in MSG_database.set_MSG])[0][0]
    elif Mtype == 4 and SG_num>2 and SG_num<16:
        magnetic_symmetries_c = magnetic_symmetries
        perm = permutation_monoclinic[SG_num-3]
        shift = origin_shift_monoclinic[SG_num-3]

        for i,p in enumerate([5]+perm):
            if p!=5:
                T_mat = LA.inv(np.array(permu_mat_monoclinic[p]).T)
                centering = dataset["international"][0]
                M = Mat2P[centering]
                T_mat = LA.inv(M).dot(T_mat).dot(M)
                tran = LA.inv(M).dot(shift[i-1])

                magnetic_symmetries_c = [ trans_magsymmop(ele,T_mat,tran) for ele in magnetic_symmetries]
            compare = [compare_magnetic_space_groups(magnetic_symmetries_c,ele,symprec) for ele in MSG_database.set_MSG]
            if np.any(compare) :
                ind = np.nonzero(compare)[0][0]
                break
        
    else:
        ind = np.nonzero([compare_magnetic_space_groups(magnetic_symmetries,ele,symprec) for ele in MSG_database.set_MSG])[0][0]
    MSG_num = MSG_database.set_MSG_num[ind].split('.')[0][16:].replace('_','.')
    
    return MSG_num,T_mat,tran,MSG_database.set_MSG[ind]
    
def magnetic_symmetries(struct,symprec=1e-5):

    # dataset = get_symmetry_dataset(struct,symprec)
    # 
    # rotations = dataset["rotations"]
    # translations = dataset["translations"]
    # time_reversal = dataset["time_reversal"]
    # 
    # symmetries = [MagSymmOp.from_rotation_and_translation_and_time_reversal(rot,trans,t) 
    #                        for rot,trans,t in zip(rotations,translations,time_reversal)]
    _, _, _, symmetries = _compare_with_magnetic_database(struct,symprec=symprec)
    return symmetries
    

def get_symmetry_dataset(struct,symprec=1e-5):
    _set_no_error()
    
    struct = standardize_magmoms(struct)
    Analyzer = SpacegroupAnalyzer(struct)
    cell = Analyzer._cell
    
    lattice, positions, numbers, magmoms = _expand_cell(cell)
    if lattice is None:
        return None

    spg_ds = mspg.dataset(lattice, positions, magmoms, numbers, symprec)
    if spg_ds is None:
        _set_error_message()
        return None

    keys = ('number',
            'hall_number',
            'international',
            'hall',
            'choice',
            'transformation_matrix',
            'origin_shift',
            'rotations',
            'translations',
            'time_reversal',
            'primitive_lattice',
            'mapping_to_primitive',
            'primitive_types',
            'positions',
            'magmoms',
            'pointgroup',
            'msgtype')
    dataset = {}
    for key, data in zip(keys, spg_ds):
        dataset[key] = data

    dataset['international'] = dataset['international'].strip()
    dataset['hall'] = dataset['hall'].strip()
    dataset['choice'] = dataset['choice'].strip()
    dataset['transformation_matrix'] = np.array(
        dataset['transformation_matrix'], dtype='double', order='C')
    dataset['origin_shift'] = np.array(dataset['origin_shift'], dtype='double')
    dataset['rotations'] = np.array(dataset['rotations'],
                                    dtype='intc', order='C')
    dataset['translations'] = np.array(dataset['translations'],
                                       dtype='double', order='C')
    dataset['time_reversal'] = np.array(dataset['time_reversal'],
                                       dtype='intc', order='C')
    dataset['primitive_lattice'] = np.array(np.transpose(dataset['primitive_lattice']),
                                            dtype='double', order='C')
    dataset['mapping_to_primitive'] = np.array(dataset['mapping_to_primitive'],
                                               dtype='intc')
    dataset['primitive_types'] = np.array(dataset['primitive_types'], dtype='intc')
    dataset['positions'] = np.array(dataset['positions'],
                                        dtype='double', order='C')
    dataset['magmoms'] = np.array(dataset['magmoms'],
                                        dtype='double', order='C')
    dataset['pointgroup'] = dataset['pointgroup'].strip()

    _set_error_message()
    return dataset
    
class MSG():
    def __init__(self, num_former, Mtype):
        self.num_former = num_former
        self.Mtype = Mtype
        self.set_MSG,self.set_MSG_num= self._get_set_MSG()
    def _get_set_MSG(self):
        ind = [f'Magnetic_Sym_El_{self.num_former}_' in ele for ele in data_filename ]
        filenames = np.array(data_filename)[ind]
        set_sp_ops = []
        MSG_nums = []
        for ele in filenames:
            with open(f'{datafilepath}{ele}') as f:
                data = [line.split() for line in f]
                data = np.array(data[1:]).astype(np.float)
                rotation = data[:,:9].reshape((len(data),3,3)).astype(np.int)
                translation = data[:,9:12]
                time_reversal = data[:,-1].astype(np.int)
                
            sp_ops = [ MagSymmOp.from_rotation_and_translation_and_time_reversal(rot,tran,tim) 
                      for rot,tran,tim in zip(rotation,np.mod(translation,1),time_reversal)]
            if MSG_type(sp_ops) == self.Mtype:
                set_sp_ops.append(sp_ops)
                MSG_nums.append(ele)
        return set_sp_ops,MSG_nums
   
    
def MSG_type(sp_ops):
    D = [ ele for ele in sp_ops if ele.time_reversal==1 ]
    DT = [ ele for ele in sp_ops if ele.time_reversal==-1 ]
    if DT == []:
        Mtype = 1
    elif np.allclose(DT[0].affine_matrix,np.eye(4)):
        Mtype = 2
    elif np.any([ np.allclose(ele.rotation_matrix,np.eye(3)) and not np.allclose(ele.translation_vector,[0,0,0]) for ele in DT]):
        Mtype = 4
    else :
        Mtype = 3
        
    return Mtype

def compare_magnetic_space_groups(MSGA,MSGB,symperc=1e-10):
    def compare(OperationA,OperationB):
        rot_trans_eq = np.allclose(OperationA.affine_matrix,OperationB.affine_matrix,atol = symperc*1e-1) 
        time_eq = OperationA.time_reversal==OperationB.time_reversal
        result = rot_trans_eq and time_eq
        return result
    result = np.all([ np.any([ compare(A,B) for A in MSGA ]) for B in MSGB])
    return result
        
    
def standardize_magmoms(struct):
    if "magmom" in struct.site_properties.keys():
        mags = struct.site_properties["magmom"]
    else:
        mags = [0.0,]*len(struct)
    struct_n = struct.copy()
    if isinstance(mags[0],Magmom):
        mags = list(map(lambda x: x.global_moment,mags))
        struct_n.add_site_property("magmom",mags)
    if np.array(mags[0]).size == 1:
        mags = [[0.0,0.0,ele] for ele in mags]
        struct_n.add_site_property("magmom",mags)
    return struct_n


def trans_magsymmop(symmop,Tmat,tau = [0.0,0.0,0.0]):
    rot = Tmat.dot(symmop.rotation_matrix).dot(LA.inv(Tmat))
    trans = Tmat.dot(symmop.translation_vector)
    trans = trans+np.array(tau)-rot.dot(tau)
    trans = np.mod(np.around(trans,5),1)
    return MagSymmOp.from_rotation_and_translation_and_time_reversal(rot,trans,symmop.time_reversal)

def _expand_cell(cell):
    latt = np.array(np.transpose(cell[0]), dtype='double', order='C')
    pos = np.array(cell[1], dtype='double', order='C')
    types = np.array(cell[2], dtype='intc')
    magmom = np.array(cell[3], order='C', dtype='double')
    return (latt,pos,types,magmom)

def get_error_message():
    return spglib_error.message


def _set_error_message():
    spglib_error.message = mspg.error_message()


def _set_no_error():
    spglib_error.message = "no error"

def lattice_mod(coords,lattice,symprec=0.01):
    lattice_abs = LA.norm(lattice,axis = 1)
    diff = np.abs(coords-1)
    diff_cart = diff*lattice_abs
    coords_n = coords.copy()
    coords_n[diff_cart<symprec] = 0.0
    return coords_n

def rational(t):
    if np.isclose(t,0,atol=2e-3):
        r = 0
    else:
        for i in range(1,1000):
            if np.isclose(i*t,np.round(i*t),atol = 2e-3):
                r = np.round(i*t)/i
                break
        else:
            r = t
    return r

def magsym_op_frac2cart(op_in,latt):
    L = latt.T.copy()
    mat_af = op_in.affine_matrix.copy()
    mat_af[0:3,0:3] = np.dot(np.dot(L,mat_af[0:3,0:3]),LA.inv(L))
    mat_af[0:3,3] = np.dot(L,mat_af[0:3,3])
    return MagSymmOp(mat_af,op_in.time_reversal)

def symmetrize_structure(struct_p,symprec = 0.01):       
    data = SpacegroupAnalyzer(struct_p,symprec = symprec).get_symmetry_dataset()
    P = data['transformation_matrix']
    rot_mat = data['std_rotation_matrix']
    latt = data['std_lattice'].T.dot(P).T
    latt_rot = LA.inv(rot_mat).dot(latt.T).T 
    ind = [ np.allclose(ele,np.eye(3),atol = symprec) for ele in data['rotations']]
    
    pure_trans = np.array([ [ rational(i) for i in t] for t in data['translations'][ind]])
    #print(data['std_positions'].round(3))    
    
    coords = data['std_positions']
    coords = LA.inv(P).dot(coords.T).T
    origin_shift = np.array([ rational(i) for i in LA.inv(P).dot(data['origin_shift'])])
    coords = coords-origin_shift
    coords = np.vstack([ np.mod(coords+t,1) for t in pure_trans ])
    coords_n = coords.copy()
    coords = lattice_mod(coords,latt,symprec=symprec)
    coords_cart = coords.dot(latt_rot)
    _,ind =  np.unique(coords_cart.round(2),axis = 0, return_index = True)
    coords = coords[ind]
    coords_n = coords_n[ind]

    coords0 = np.mod(struct_p.frac_coords,1)
    coords0 = lattice_mod(coords0,struct_p.lattice.matrix,symprec=symprec)

    mags = np.array([ rot_mat.dot(m) for m in struct_p.site_properties['magmom']])
    coords0_cart = coords0.dot(struct_p.lattice.matrix)
    coords_cart = coords.dot(latt_rot)
    atol = np.max(LA.norm(struct_p.lattice.matrix,axis=1))*symprec
    ind = [ np.isclose(coords0_cart,ele,atol=atol).all(axis=1) for ele in coords_cart]
    check = np.all([ np.sum(i.astype(int))==1 for i in ind])
    if not check:
        raise IndexError('Can not match all coordinates of standard and original structure')
    
    coords = np.array([ coords_n[i][0] for i in np.array(ind).T])
    
    return Structure(latt,struct_p.species,coords,site_properties = {"magmom":mags})



def site_symmetry(struct_p,symmetries,symprec=0.01):
    ind = np.arange(len(struct_p))
    ind_sym = np.arange(len(symmetries))
    equal_sites_ind = []
    rep_element_ind = []
    coords0 = np.around([ele.frac_coords for ele in struct_p],4)
    coords0 = lattice_mod(coords0,struct_p.lattice.matrix,symprec=symprec)
    coords_input = coords0
    site_sym_ind = []
    while len(ind)>0:
        coords0 = coords_input[ind]
        coords0_cart = coords0.dot(struct_p.lattice.matrix)
        coords = np.array([ np.mod(op.operate(coords0[0]),1) for op in symmetries])
        coords = lattice_mod(coords,struct_p.lattice.matrix,symprec=symprec)
        coords_cart = coords.dot(struct_p.lattice.matrix)
        mapping = np.array([ np.isclose(coords_cart,ele,atol=9*symprec).all(axis = 1) for ele in coords0_cart])
        equal_sites_ind.append(ind[mapping.any(axis = 1)])
        site_sym_ind.append(ind_sym[mapping[0]])
        rep_element_ind.append([ np.nonzero(ele)[0] for ele in mapping[1:][mapping[1:].any(axis=1)] ] )
        # remove indices of grouped elements
        ind = np.array([ i for i in ind if i not in equal_sites_ind[-1]])
    
    g = len(symmetries)
    check = []
    for  equal_ind,site_ind,rep_ind in zip(equal_sites_ind,site_sym_ind,rep_element_ind):
        check.append(len(site_ind)== g/len(equal_ind) and len(rep_ind)==len(equal_ind)-1 and np.all([ len(ele)==g/len(equal_ind) for ele in rep_ind]))
    if not np.all(check):
        raise ValueError('Wrong decomposition into site symmetries')
    mags = np.array(struct_p.site_properties['magmom'])
    symmetries_cart = [ magsym_op_frac2cart(op,struct_p.lattice.matrix) for op in symmetries]
    return equal_sites_ind,site_sym_ind,rep_element_ind
    
    

def symmetrize_magmom(struct_p,symmetries,symprec=0.01,nround = 10):
    # struct_p must be a standard primitive structure 
    # and magnetic moment is set by a global axis in form of array
    
    # round magnetic moments first
    magmoms = np.around(struct_p.site_properties['magmom'],nround)
    struct_p.add_site_property('magmom',magmoms)

    symmetries_cart_mag = [ magsym_op_frac2cart(op,struct_p.lattice.matrix) for op in symmetries]
    time_reversals = [ op.time_reversal for op in symmetries]
    symmetries_cart_mag = [ t*LA.det(ele.rotation_matrix)*ele.rotation_matrix 
                           for ele,t in zip(symmetries_cart_mag,time_reversals)]
    equal_sites_ind,site_sym_ind,rep_element_ind = site_symmetry(struct_p,symmetries,symprec=symprec)    
    
    magmoms = []
    for site_ind,rep_ind in zip(equal_sites_ind,rep_element_ind):
        magmom = [struct_p[site_ind[0]].magmom] + [LA.inv(symmetries_cart_mag[rep[0]]).dot(struct_p[i].magmom)
                   for rep,i in zip(rep_ind,site_ind[1:])]
        #magmoms = [ m.moment for m in magmoms]
        m = [struct_p[site_ind[0]].magmom] + [struct_p[i].magmom
                   for rep,i in zip(rep_ind,site_ind[1:])]
        magmoms.append(np.average(magmom,axis=0))
    magmom_ideal = []
    # constrain magnetic moments by site symmetries
    
    for ind,magmom in zip(site_sym_ind,magmoms):
        if not np.allclose(magmom,[0,0,0],atol=1e-8):
            mag_axis0 = magmom/LA.norm(magmom)
            #print(magmom)
            for i in ind:
                if not np.allclose(LA.det(symmetries_cart_mag[i])*symmetries_cart_mag[i],np.eye(3)):
                    #print('here')
                    v,axis = LA.eig(np.around(symmetries_cart_mag[i],6))
                    mag_axis = axis.T[np.isclose(v,1)].real
                    if len(mag_axis)>1:
                        magmom_ideal.append(np.sum([c*m for c,m in zip(mag_axis.dot(magmom),mag_axis)],axis=0))
                    else:
                        magmom_ideal.append(np.sign(magmom.dot(mag_axis[0]))*mag_axis[0]*LA.norm(magmom))
                    break
            else:
                magmom_ideal.append(magmom)
        else:
            magmom_ideal.append([0.0,0.0,0.0])

    magmoms = np.zeros((len(struct_p),3))
    for ind,reps,mag in zip(equal_sites_ind,rep_element_ind,magmom_ideal):
        magmoms[ind[0]] = mag
        if len(ind)>1:
            magmoms[ind[1:]] = np.array([symmetries_cart_mag[r[0]].dot(mag) for r in reps])
    magmoms = magmoms.round(nround+2)            
    return magmoms


def symmetrize(struct,symprec=0.01,nround = 10):
    # magnetic symmetry may be improved after only symmetrized nonmagnetic lattice
    # symmetries is obtained with respect to standard primitive cell of symmetrized structure
    # Primitive cell of non symmetrized structure does not match symmetries of standard primitive cell of symmetrized structure
    struct_p = standard_primitive(struct,symprec = symprec)
    #print(mspg.magnetic_spacegroup_number(struct,symprec = symprec))
    symmetries = magnetic_symmetries(struct,symprec = symprec)
    struct_p = symmetrize_structure(struct_p,symprec = symprec)
    magmoms = symmetrize_magmom(struct_p,symmetries,symprec = symprec,nround = nround)
    struct_p.add_site_property('magmom',magmoms)
    return struct_p
