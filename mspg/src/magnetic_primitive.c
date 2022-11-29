#include <stdio.h>
#include <stdlib.h>
#include "cell.h"
#include "magnetic_cell.h"
#include "delaunay.h"
#include "mathfunc.h"
#include "primitive.h"
#include "magnetic_primitive.h"
#include "magnetic_symmetry.h"
#include "symmetry.h"
#include "debug.h"

#ifdef DMALLOC
#include "dmalloc.h"
#endif

#define REDUCE_RATE 0.95
#define NUM_ATTEMPT 20

static Magnetic_Symmetry * collect_magnetic_primitive_symmetry(const Magnetic_Symmetry *symmetry,
                                             const int primsym_size);

static Cell * get_cell_with_smallest_lattice(const Cell * cell,
                                             const double symprec);

Cell * relabeling_types_according_magmom(const Magnetic_Cell *cell, const double symprec);

static VecDBL * collect_pure_translations(const Symmetry *symmetry);

static int vector_check_d3(const double a[3], const double b[3], const double symprec);
static int maximum(int arr[], int n);

static int get_primitive_in_translation_space(double t_mat_inv[3][3],
                                              const VecDBL *pure_trans,
                                              const int symmetry_size,
                                              const double symprec);

static Primitive * get_primitive(const Cell * cell,
                                 const double symprec,
                                 const double angle_tolerance);

static Cell * get_primitive_cell(int * mapping_table,
                                 const Cell * cell,
                                 const VecDBL * pure_trans,
                                 const double symprec,
                                 const double angle_tolerance);

static int get_primitive_lattice_vectors(double prim_lattice[3][3],
                                         const Cell * cell,
                                         const VecDBL * pure_trans,
                                         const double symprec,
                                         const double angle_tolerance);

static int find_primitive_lattice_vectors(double prim_lattice[3][3],
                                          const VecDBL * vectors,
                                          const Cell * cell,
                                          const double symprec);

/* Return 0 if failed */
static int get_primitive_lattice_vectors(double prim_lattice[3][3],
                                         const Cell * cell,
                                         const VecDBL * pure_trans,
                                         const double symprec,
                                         const double angle_tolerance);

static int restore_types(Magnetic_Primitive *primitive,const Magnetic_Cell *cell_old,const Cell *cell_new);

static VecDBL * get_translation_candidates(const VecDBL * pure_trans);


Magnetic_Primitive * prm_get_magnetic_primitive(const Magnetic_Cell * cell,
                                       const double symprec,
                                       const double angle_tolerance);
Primitive * Magnetic_Primitive_to_Primitive(const Magnetic_Primitive * primitive);

Magnetic_Primitive * set_magnetic_primitive(const Primitive * primitive, const double magmom[][3]);

Magnetic_Primitive * magn_prm_alloc_primitive(const int size);

void magn_prm_free_primitive(Magnetic_Primitive * primitive);

Magnetic_Primitive * prm_get_magnetic_primitive(const Magnetic_Cell * cell,
                                       const double symprec,
                                       const double angle_tolerance)
{
  Cell *cell_new;
  Magnetic_Primitive *primitive_new;
  Primitive *primitive;
  
  cell_new = NULL;
  primitive_new = NULL;
  primitive = NULL;

  cell_new = relabeling_types_according_magmom(cell,symprec);

  primitive = prm_get_primitive(cell_new,symprec,angle_tolerance);
  

  if ((primitive_new = set_magnetic_primitive(primitive,cell->magmom) ) == NULL){
    prm_free_primitive(primitive);
    return NULL;
  }
  

  if (!restore_types(primitive_new,cell,cell_new))
  {
    prm_free_primitive(primitive);
    return NULL;
  }

  prm_free_primitive(primitive);
  cel_free_cell(cell_new);
  

  return primitive_new;
}

static int restore_types(Magnetic_Primitive *primitive,const Magnetic_Cell *cell_old,const Cell *cell_new)
{
  int *types_matching;
  int i,ind;
  if ((types_matching = (int *) malloc(sizeof(int)*maximum(cell_new->types,cell_new->size))) == NULL) {
    warning_print("spglib: Memory could not be allocated ");
    return 0;
  }

  for (i=0; i<maximum(cell_new->types,cell_new->size); i++)
  {
    types_matching[i] = -1;
  }
  for (i=0; i<cell_old->size; i++)
  { ind = cell_new->types[i]-1;
    if (types_matching[ind]==-1)
    {
      types_matching[ind] = cell_old->types[i];
    }
  }


  for (i=0; i<primitive->cell->size; i++)
  { ind = primitive->cell->types[i]-1;
    primitive->cell->types[i] = types_matching[ind];
  }

  free(types_matching);
  return 1;
}

Magnetic_Primitive * set_magnetic_primitive(const Primitive * primitive, const double magmom[][3])
{
  Magnetic_Primitive *primitive_new;
  int i,ind;
  
  if ((primitive_new = magn_prm_alloc_primitive(primitive->size)) == NULL) {
    return NULL;
  }   
  primitive_new->cell = initialize_magnetic_cell_from_cell(primitive->cell);

  primitive_new->tolerance = primitive->tolerance;
  primitive_new->angle_tolerance = primitive->angle_tolerance;
  primitive_new->size = primitive->size;

  if ((primitive_new->orig_lattice =
       (double (*)[3]) malloc(sizeof(double[3]) * 3)) == NULL) {
    warning_print("spglib: Memory could not be allocated.");
    return NULL;
  }


  mat_copy_matrix_d3(primitive_new->orig_lattice,primitive->orig_lattice);


  for (i=0;i<primitive->size;i++)
  {
    ind = primitive->mapping_table[i];
    mat_copy_vector_d3(primitive_new->cell->magmom[ind],magmom[i]);
  }

  for (i=0;i<primitive->size;i++)
  {
    primitive_new->mapping_table[i] = primitive->mapping_table[i];
  }
  
  return primitive_new;

}


Primitive * Magnetic_Primitive_to_Primitive(const Magnetic_Primitive * primitive)
{ 
  Primitive *primitive_new;
  
 
  if ((primitive_new = prm_alloc_primitive(primitive->size)) == NULL) {
    return NULL;
  }

 
  primitive_new->cell = Magnetic_Cell_to_Cell(primitive->cell);
  
  if ((primitive_new->orig_lattice =
       (double (*)[3]) malloc(sizeof(double[3]) * 3)) == NULL) {
    warning_print("spglib: Memory could not be allocated.");
    return NULL;
  }

  mat_copy_matrix_d3(primitive_new->orig_lattice,primitive->orig_lattice);

  primitive_new->size = primitive->size;
  primitive_new->tolerance = primitive->tolerance;
  primitive_new->angle_tolerance = primitive->angle_tolerance;


  return primitive_new;
}

Magnetic_Symmetry * prm_get_magnetic_primitive_symmetry(const Magnetic_Symmetry *symmetry,
                                      const double symprec)
{
  int i, primsym_size;
  VecDBL *pure_trans;
  Magnetic_Symmetry *prim_symmetry;
  Symmetry *symmetry_new;
  double t_mat[3][3], t_mat_inv[3][3], tmp_mat[3][3];

  pure_trans = NULL;
  prim_symmetry = NULL;
  symmetry_new = magnetic_symmetry_to_symmetry(symmetry);

  if ((pure_trans = collect_pure_translations(symmetry_new)) == NULL) {
    return NULL;
  }
  primsym_size = symmetry->size / pure_trans->size;

  /* t_mat: T=(Lp^-1.L) where L is identity matrix. */
  if (get_primitive_in_translation_space
      (t_mat_inv, pure_trans, symmetry->size, symprec) == 0) {
    mat_free_VecDBL(pure_trans);
    pure_trans = NULL;
    return NULL;
  }

  mat_free_VecDBL(pure_trans);
  pure_trans = NULL;

  if (!mat_inverse_matrix_d3(t_mat, t_mat_inv, symprec)) {
    return NULL;
  }

  /* Collect operations for primitive cell from operations in 'symmetry' */
  if ((prim_symmetry = collect_magnetic_primitive_symmetry(symmetry, primsym_size)) ==
      NULL) {
    return NULL;
  }

  /* Overwrite prim_symmetry by R_p = TRT^-1, t_p = T.t */
  for (i = 0; i < prim_symmetry->size; i++) {
    mat_multiply_matrix_di3(tmp_mat, t_mat, prim_symmetry->rot[i]);
    mat_multiply_matrix_d3(tmp_mat, tmp_mat, t_mat_inv);
    mat_cast_matrix_3d_to_3i(prim_symmetry->rot[i], tmp_mat);
    mat_multiply_matrix_vector_d3(prim_symmetry->trans[i],
                                  t_mat, prim_symmetry->trans[i]);
  }

#ifdef SPGDEBUG
  int j;
  for (i = 0; i < prim_symmetry->size; i++) {
    fprintf(stderr, "--- %d ---\n", i + 1);
    for (j = 0; j < 3; j++) {
      fprintf(stderr, "%d %d %d\n",
              prim_symmetry->rot[i][j][0],
              prim_symmetry->rot[i][j][1],
              prim_symmetry->rot[i][j][2]);
    }
    fprintf(stderr, "%f %f %f\n",
            prim_symmetry->trans[i][0],
            prim_symmetry->trans[i][1],
            prim_symmetry->trans[i][2]);
  }
#endif

  return prim_symmetry;
}

static Magnetic_Symmetry * collect_magnetic_primitive_symmetry(const Magnetic_Symmetry *symmetry,
                                             const int primsym_size)
{
  int i, j, num_psym, is_found;
  Magnetic_Symmetry *prim_symmetry;

  prim_symmetry = NULL;

  prim_symmetry = sym_alloc_magnetic_symmetry(primsym_size);
  num_psym = 1;
  mat_copy_matrix_i3(prim_symmetry->rot[0], symmetry->rot[0]);
  mat_copy_vector_d3(prim_symmetry->trans[0], symmetry->trans[0]);
  for (i = 1; i < symmetry->size; i++) {
    is_found = 1;
    for (j = 0; j < num_psym; j++) {
      if (mat_check_identity_matrix_i3(prim_symmetry->rot[j], symmetry->rot[i]))
      {
        is_found = 0;
        break;
      }
    }
    if (is_found) {
      if (num_psym == primsym_size) {
        sym_free_magnetic_symmetry(prim_symmetry);
        prim_symmetry = NULL;
        break;
      }
      mat_copy_matrix_i3(prim_symmetry->rot[num_psym], symmetry->rot[i]);
      mat_copy_vector_d3(prim_symmetry->trans[num_psym], symmetry->trans[i]);
      prim_symmetry->time_reversal[num_psym]=symmetry->time_reversal[i];
      num_psym++;
    }
  }

  if (prim_symmetry == NULL) {
    return NULL;
  }

  if (num_psym != primsym_size) {
    sym_free_magnetic_symmetry(prim_symmetry);
    prim_symmetry = NULL;
    return NULL;
  }

  return prim_symmetry;
}

Cell * relabeling_types_according_magmom(const Magnetic_Cell *cell, const double symprec)
{
  Cell *cell_new;
  int *types;
  int i,j,num_independent;
  
  cell_new = NULL;

  if ((types = (int *) malloc(sizeof(int)*cell->size)) == NULL) {
    warning_print("spglib: Memory could not be allocated ");
    return NULL;
  }

  cell_new = Magnetic_Cell_to_Cell(cell);


  for (i=0;i<cell->size;i++)
  {
    types[i] = -1;
  }
  

  num_independent = 1;
  for (i=0; i<cell->size; i++)
  { 
    if (types[i]==-1){
      for (j=i;j<cell->size;j++)
      { 
      
        if ((cell->types[i] == cell->types[j]) && vector_check_d3(cell->magmom[i],cell->magmom[j],symprec))
          {
            types[j] = num_independent;
        
          }

      }
      
      num_independent++;
    }else
    {   

        continue;
    }
  }

  cell_new->size = cell->size ;
  mat_copy_matrix_d3(cell_new->lattice,cell->lattice);
  for (i=0; i<cell->size; i++)
  {
    mat_copy_vector_d3(cell_new->position[i],cell->position[i]);
    cell_new->types[i] = types[i];
  }
  
  free(types);
  return cell_new;
}

static int vector_check_d3(const double a[3], const double b[3], const double symprec)
{
  
  if ( mat_Dabs( a[0] - b[0] ) > symprec ||
       mat_Dabs( a[1] - b[1] ) > symprec ||
       mat_Dabs( a[2] - b[2] ) > symprec  ) {
    return 0;
  }
  else{
    return 1;
  }

  
}


static VecDBL * collect_pure_translations(const Symmetry *symmetry)
{
  int i, num_pure_trans;
  VecDBL *pure_trans;
  VecDBL *ret_pure_trans;
  static int identity[3][3] = {{ 1, 0, 0 },
                               { 0, 1, 0 },
                               { 0, 0, 1 }};
  num_pure_trans = 0;
  pure_trans = NULL;
  ret_pure_trans = NULL;

  if ((pure_trans = mat_alloc_VecDBL(symmetry->size)) == NULL) {
    return NULL;
  }

  for (i = 0; i < symmetry->size; i++) {
    if (mat_check_identity_matrix_i3(symmetry->rot[i], identity)) {
      mat_copy_vector_d3(pure_trans->vec[num_pure_trans], symmetry->trans[i]);
      num_pure_trans++;
    }
  }

  if ((ret_pure_trans = mat_alloc_VecDBL(num_pure_trans)) == NULL) {
    mat_free_VecDBL(pure_trans);
    pure_trans = NULL;
    return NULL;
  }

  for (i = 0; i < num_pure_trans; i++) {
    mat_copy_vector_d3(ret_pure_trans->vec[i], pure_trans->vec[i]);
  }

  mat_free_VecDBL(pure_trans);
  pure_trans = NULL;

  return ret_pure_trans;
}

static int get_primitive_in_translation_space(double t_mat_inv[3][3],
                                              const VecDBL *pure_trans,
                                              const int symmetry_size,
                                              const double symprec)
{
  int i, j, primsym_size;
  Primitive *primitive;
  Cell *cell;

  cell = NULL;
  primitive = NULL;

  if ((cell = cel_alloc_cell(pure_trans->size)) == NULL) {
    return 0;
  }

  primsym_size = symmetry_size / pure_trans->size;
  if (symmetry_size != primsym_size * pure_trans->size) {
    cel_free_cell(cell);
    cell = NULL;
    return 0;
  }

  for (i = 0; i < pure_trans->size; i++) {
    cell->types[i] = 1;
    for (j = 0; j < 3; j++) {
      cell->position[i][j] = pure_trans->vec[i][j];
    }
  }

  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      if (i == j) {
        cell->lattice[i][j] = 1;
      } else {
        cell->lattice[i][j] = 0;
      }
    }
  }

  primitive = get_primitive(cell, symprec, -1.0);
  cel_free_cell(cell);
  cell = NULL;

  if (primitive->cell->size != 1) {
    prm_free_primitive(primitive);
    primitive = NULL;
    return 0;
  }

  mat_copy_matrix_d3(t_mat_inv, primitive->cell->lattice);
  prm_free_primitive(primitive);
  primitive = NULL;

  return 1;
}

/* Return NULL if failed */
static Primitive * get_primitive(const Cell * cell,
                                 const double symprec,
                                 const double angle_tolerance)
{
  int i, attempt;
  double tolerance;
  Primitive *primitive;
  VecDBL * pure_trans;

  debug_print("get_primitive (tolerance = %f):\n", symprec);

  primitive = NULL;
  pure_trans = NULL;

  if ((primitive = prm_alloc_primitive(cell->size)) == NULL) {
    goto notfound;
  }

  tolerance = symprec;
  for (attempt = 0; attempt < NUM_ATTEMPT; attempt++) {
    debug_print("get_primitive (attempt = %d):\n", attempt);
    if ((pure_trans = sym_get_pure_translation(cell, tolerance)) != NULL) {
      if (pure_trans->size == 1) {
        if ((primitive->cell = get_cell_with_smallest_lattice(cell, tolerance))
            != NULL) {
          for (i = 0; i < cell->size; i++) {
            primitive->mapping_table[i] = i;
          }
          goto found;
        }
      } else {
        if ((primitive->cell = get_primitive_cell(primitive->mapping_table,
                                                  cell,
                                                  pure_trans,
                                                  tolerance,
                                                  angle_tolerance)) != NULL) {
          goto found;
        }
      }
    }

    mat_free_VecDBL(pure_trans);
    pure_trans = NULL;

    tolerance *= REDUCE_RATE;
    warning_print("spglib: Reduce tolerance to %f ", tolerance);
    warning_print("(line %d, %s).\n", __LINE__, __FILE__);
  }

  prm_free_primitive(primitive);
  primitive = NULL;

 notfound:
  return NULL;

 found:
  primitive->tolerance = tolerance;
  primitive->angle_tolerance = angle_tolerance;
  if ((primitive->orig_lattice =
       (double (*)[3]) malloc(sizeof(double[3]) * 3)) == NULL) {
    warning_print("spglib: Memory could not be allocated.");
    return NULL;
  }
  mat_copy_matrix_d3(primitive->orig_lattice, cell->lattice);
  mat_free_VecDBL(pure_trans);
  pure_trans = NULL;
  return primitive;
}

/* Return NULL if failed */
static Cell * get_cell_with_smallest_lattice(const Cell * cell,
                                             const double symprec)
{
  int i, j;
  double min_lat[3][3], trans_mat[3][3], inv_lat[3][3];
  Cell * smallest_cell;

  debug_print("get_cell_with_smallest_lattice:\n");

  smallest_cell = NULL;

  if (!del_delaunay_reduce(min_lat, cell->lattice, symprec)) {
    return NULL;
  }

  mat_inverse_matrix_d3(inv_lat, min_lat, 0);
  mat_multiply_matrix_d3(trans_mat, inv_lat, cell->lattice);

  if ((smallest_cell = cel_alloc_cell(cell->size)) == NULL) {
    return NULL;
  }

  mat_copy_matrix_d3(smallest_cell->lattice, min_lat);
  for (i = 0; i < cell->size; i++) {
    smallest_cell->types[i] = cell->types[i];
    mat_multiply_matrix_vector_d3(smallest_cell->position[i],
                                  trans_mat, cell->position[i]);
    for (j = 0; j < 3; j++) {
      smallest_cell->position[i][j] = mat_Dmod1(smallest_cell->position[i][j]);
    }
  }

  return smallest_cell;
}

/* Return NULL if failed */
static Cell * get_primitive_cell(int * mapping_table,
                                 const Cell * cell,
                                 const VecDBL * pure_trans,
                                 const double symprec,
                                 const double angle_tolerance)
{
  int multi;
  double prim_lattice[3][3];
  Cell * primitive_cell;

  debug_print("get_primitive_cell:\n");

  primitive_cell = NULL;

  /* Primitive lattice vectors are searched. */
  /* To be consistent, sometimes tolerance is decreased iteratively. */
  /* The descreased tolerance is stored in 'static double tolerance'. */
  multi = get_primitive_lattice_vectors(prim_lattice,
                                        cell,
                                        pure_trans,
                                        symprec,
                                        angle_tolerance);
  if (! multi) {
    goto not_found;
  }

  /* Fit atoms into new primitive cell */
  if ((primitive_cell = cel_trim_cell(mapping_table,
                                      prim_lattice,
                                      cell,
                                      symprec)) == NULL) {
    goto not_found;
  }

  /* found */
  return primitive_cell;

 not_found:
  warning_print("spglib: Primitive cell could not be found ");
  warning_print("(line %d, %s).\n", __LINE__, __FILE__);
  return NULL;
}

/* Return 0 if failed */
static int get_primitive_lattice_vectors(double prim_lattice[3][3],
                                         const Cell * cell,
                                         const VecDBL * pure_trans,
                                         const double symprec,
                                         const double angle_tolerance)
{
  int i, multi, attempt;
  double tolerance;
  VecDBL * vectors, * pure_trans_reduced, *tmp_vec;

  vectors = NULL;
  pure_trans_reduced = NULL;
  tmp_vec = NULL;

  tolerance = symprec;

  if ((pure_trans_reduced = mat_alloc_VecDBL(pure_trans->size)) == NULL) {
    goto fail;
  }

  for (i = 0; i < pure_trans->size; i++) {
    mat_copy_vector_d3(pure_trans_reduced->vec[i], pure_trans->vec[i]);
  }

  for (attempt = 0; attempt < NUM_ATTEMPT; attempt++) {
    multi = pure_trans_reduced->size;

    if ((vectors = get_translation_candidates(pure_trans_reduced)) == NULL) {
      mat_free_VecDBL(pure_trans_reduced);
      pure_trans_reduced = NULL;
      goto fail;
    }

    /* Lattice of primitive cell is found among pure translation vectors */
    if (find_primitive_lattice_vectors(prim_lattice,
                                       vectors,
                                       cell,
                                       tolerance)) {
      mat_free_VecDBL(vectors);
      vectors = NULL;
      mat_free_VecDBL(pure_trans_reduced);
      pure_trans_reduced = NULL;

      if (! del_delaunay_reduce(prim_lattice, prim_lattice, symprec)) {
        goto fail;
      }

      goto found;

    } else {

      if ((tmp_vec = mat_alloc_VecDBL(multi)) == NULL) {
        mat_free_VecDBL(vectors);
        vectors = NULL;
        mat_free_VecDBL(pure_trans_reduced);
        pure_trans_reduced = NULL;
        goto fail;
      }

      for (i = 0; i < multi; i++) {
        mat_copy_vector_d3(tmp_vec->vec[i], pure_trans_reduced->vec[i]);
      }
      mat_free_VecDBL(pure_trans_reduced);
      pure_trans_reduced = NULL;

      pure_trans_reduced = sym_reduce_pure_translation(cell,
                                                       tmp_vec,
                                                       tolerance,
                                                       angle_tolerance);

      mat_free_VecDBL(tmp_vec);
      tmp_vec = NULL;
      mat_free_VecDBL(vectors);
      vectors = NULL;

      if (pure_trans_reduced == NULL) {
        goto fail;
      }

      warning_print("spglib: Tolerance is reduced to %f (%d), ",
                    tolerance, attempt);
      warning_print("num_pure_trans = %d\n", pure_trans_reduced->size);

      tolerance *= REDUCE_RATE;
    }
  }

  mat_free_VecDBL(pure_trans_reduced);
  pure_trans_reduced = NULL;

 fail:
  return 0;

 found:
  return multi;
}

/* Return 0 if failed */
static int find_primitive_lattice_vectors(double prim_lattice[3][3],
                                          const VecDBL * vectors,
                                          const Cell * cell,
                                          const double symprec)
{
  int i, j, k, size;
  double initial_volume, volume;
  double relative_lattice[3][3], min_vectors[3][3], tmp_lattice[3][3];
  double inv_mat_dbl[3][3];
  int inv_mat_int[3][3];

  debug_print("find_primitive_lattice_vectors:\n");

  size = vectors->size;
  initial_volume = mat_Dabs(mat_get_determinant_d3(cell->lattice));

  /* check volumes of all possible lattices, find smallest volume */
  for (i = 0; i < size; i++) {
    for (j = i + 1; j < size; j++) {
      for (k = j + 1; k < size; k++) {
        mat_multiply_matrix_vector_d3(tmp_lattice[0],
                                      cell->lattice,
                                      vectors->vec[i]);
        mat_multiply_matrix_vector_d3(tmp_lattice[1],
                                      cell->lattice,
                                      vectors->vec[j]);
        mat_multiply_matrix_vector_d3(tmp_lattice[2],
                                      cell->lattice,
                                      vectors->vec[k]);
        volume = mat_Dabs(mat_get_determinant_d3(tmp_lattice));
        if (volume > symprec) {
          if (mat_Nint(initial_volume / volume) == size - 2) {
            mat_copy_vector_d3(min_vectors[0], vectors->vec[i]);
            mat_copy_vector_d3(min_vectors[1], vectors->vec[j]);
            mat_copy_vector_d3(min_vectors[2], vectors->vec[k]);
            goto ret;
          }
        }
      }
    }
  }

  /* Not found */
  warning_print("spglib: Primitive lattice vectors cound not be found ");
  warning_print("(line %d, %s).\n", __LINE__, __FILE__);
  return 0;

  /* Found */
 ret:
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      relative_lattice[j][i] = min_vectors[i][j];
    }
  }

  mat_inverse_matrix_d3(inv_mat_dbl, relative_lattice, 0);
  mat_cast_matrix_3d_to_3i(inv_mat_int, inv_mat_dbl);
  if (abs(mat_get_determinant_i3(inv_mat_int)) == size-2) {
    mat_cast_matrix_3i_to_3d(inv_mat_dbl, inv_mat_int);
    mat_inverse_matrix_d3(relative_lattice, inv_mat_dbl, 0);
  } else {
    warning_print("spglib: Primitive lattice cleaning is incomplete ");
    warning_print("(line %d, %s).\n", __LINE__, __FILE__);
  }
  mat_multiply_matrix_d3(prim_lattice, cell->lattice, relative_lattice);

  return 1;
}

static VecDBL * get_translation_candidates(const VecDBL * pure_trans)
{
  int i, j, multi;
  VecDBL * vectors;

  vectors = NULL;
  multi = pure_trans->size;

  if ((vectors = mat_alloc_VecDBL(multi + 2)) == NULL) {
    return NULL;
  }

  /* store pure translations in original cell */
  /* as trial primitive lattice vectors */
  for (i = 0; i < multi - 1; i++) {
    mat_copy_vector_d3(vectors->vec[i], pure_trans->vec[i + 1]);
  }

  /* store lattice translations of original cell */
  /* as trial primitive lattice vectors */
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      if (i == j) {
        vectors->vec[i+multi-1][j] = 1;
      } else {
        vectors->vec[i+multi-1][j] = 0;
      }
    }
  }

  return vectors;
}

/* return NULL if failed */
Magnetic_Primitive * magn_prm_alloc_primitive(const int size)
{
  Magnetic_Primitive *primitive;
  int i;

  primitive = NULL;

  if ((primitive = (Magnetic_Primitive*) malloc(sizeof(Magnetic_Primitive))) == NULL) {
    warning_print("spglib: Memory could not be allocated ");
    return NULL;
  }

  primitive->cell = NULL;
  primitive->mapping_table = NULL;
  primitive->size = size;
  primitive->tolerance = 0;
  primitive->angle_tolerance = -1.0;
  primitive->orig_lattice = NULL;

  if (size > 0) {
    if ((primitive->mapping_table = (int*) malloc(sizeof(int) * size)) == NULL) {
      warning_print("spglib: Memory could not be allocated ");
      warning_print("(Primitive, line %d, %s).\n", __LINE__, __FILE__);
      free(primitive);
      primitive = NULL;
      return NULL;
    }
  }

  for (i = 0; i < size; i++) {
    primitive->mapping_table[i] = -1;
  }

  return primitive;
}

void magn_prm_free_primitive(Magnetic_Primitive * primitive)
{
  if (primitive != NULL) {
    if (primitive->mapping_table != NULL) {
      free(primitive->mapping_table);
      primitive->mapping_table = NULL;
    }

    if (primitive->cell != NULL) {
      magnetic_cel_free_cell(primitive->cell);
      primitive->cell = NULL;
    }

    if (primitive->orig_lattice != NULL) {
      free(primitive->orig_lattice);
      primitive->orig_lattice = NULL;
    }

    free(primitive);
  }
}

static int maximum(int arr[], int n)
{
    int i;
     
    // Initialize maximum element
    int max = arr[0];
 
    // Traverse array elements 
    // from second and compare
    // every element with current max 
    for (i = 1; i < n; i++)
        if (arr[i] > max)
            max = arr[i];
 
    return max;
}
