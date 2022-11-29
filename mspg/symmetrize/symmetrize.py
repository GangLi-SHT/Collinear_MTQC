# get standard primitive cell 
# get standard conventional lattice
# map std_lattice and std_positions to original cell
# transform magmoms by std_rotations
# get standard primitive cell by mspglib
# get standard primitive magnetic symmetry operations
# symmetrize magmoms by symmetries
import numpy as np
from pymatgen.core.operations import MagSymmOp
from pymatgen.core.structure import Structure
import mspglib as mspg
import numpy.linalg as LA
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

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
    ind = [ np.allclose(ele,np.eye(3)) for ele in data['rotations']]
    
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
    ind = [ np.isclose(coords0_cart,ele,atol=9*symprec).all(axis=1) for ele in coords_cart]
    #print(np.array(ind).astype(int).sum(axis = 1))    
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
    
    

def symmetrize_magmom(struct_p,symmetries,symprec=0.01):
    # struct_p must be a standard primitive structure 
    # and magnetic moment is set by a global axis in form of array
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
            
    return magmoms


def symmetrize(struct,symprec=0.01):
    # magnetic symmetry may be improved after only symmetrized nonmagnetic lattice
    # symmetries is obtained with respect to standard primitive cell of symmetrized structure
    # Primitive cell of non symmetrized structure does not match symmetries of standard primitive cell of symmetrized structure
    struct_p = mspg.standard_primitive(struct,symprec = symprec)
    #print(mspg.magnetic_spacegroup_number(struct,symprec = symprec))
    symmetries = mspg.magnetic_symmetries(struct,symprec = symprec)
    struct_p = symmetrize_structure(struct_p,symprec = symprec)
    magmoms = symmetrize_magmom(struct_p,symmetries,symprec = symprec)
    struct_p.add_site_property('magmom',magmoms)
    return struct_p
