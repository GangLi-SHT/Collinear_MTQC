#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stddef.h>
#include "arithmetic.h"
#include "mathfunc.h"
#include "magnetic_primitive.h"
#include "magnetic_cell.h"
#include "magnetic_determination.h"
#include "spacegroup.h"
#include "mspacegroup.h"
#include "mspglib.h"
#include "debug.h"
#include "pointgroup.h"

#ifdef DMALLOC
#include "dmalloc.h"
#endif

static MspglibError spglib_error_code = MSPGLIB_SUCCESS;

static MspglibErrorMessage mspglib_error_message[] = {
  {MSPGLIB_SUCCESS, "no error"},
  {MSPGERR_SPACEGROUP_SEARCH_FAILED, "spacegroup search failed"},
  {MSPGERR_CELL_STANDARDIZATION_FAILED, "cell standardization failed"},
  {MSPGERR_SYMMETRY_OPERATION_SEARCH_FAILED, "symmetry operation search failed"},
  {MSPGERR_ATOMS_TOO_CLOSE, "too close distance between atoms"},
  {MSPGERR_POINTGROUP_NOT_FOUND, "pointgroup not found"},
  {MSPGERR_NIGGLI_FAILED, "Niggli reduction failed"},
  {MSPGERR_DELAUNAY_FAILED, "Delaunay reduction failed"},
  {MSPGERR_ARRAY_SIZE_SHORTAGE, "array size shortage"},
  {MSPGERR_NONE, ""},
};

static MspglibDataset * init_magnetic_dataset(void);
static int set_magnetic_dataset(SPGCONST MspglibDataset * dataset,
                                SPGCONST Magnetic_Cell * cell,
                                SPGCONST MagneticDataContainer * container);

static MspglibDataset * get_magnetic_dataset(SPGCONST double lattice[3][3],
                                      SPGCONST double position[][3],
                                      SPGCONST double magmom[][3],
                                      const int types[],
                                      const int num_atom,
                                      const int hall_number,
                                      const double symprec,
                                      const double angle_tolerance);



MspglibDataset * mspg_get_magnetic_dataset(SPGCONST double lattice[3][3],
                                          SPGCONST double position[][3],
                                          SPGCONST double magmom[][3],
                                          const int types[],
                                          const int num_atom,
                                          const double symprec)
{
  return get_magnetic_dataset(lattice,
                              position,
                              magmom,
                              types,
                              num_atom,
                              0,
                              symprec,
                              -1.0);
}


/* Return NULL if failed */
static MspglibDataset * get_magnetic_dataset(SPGCONST double lattice[3][3],
                                             SPGCONST double position[][3],
                                             SPGCONST double magmom[][3],
                                             const int types[],
                                             const int num_atom,
                                             const int hall_number,
                                             const double symprec,
                                             const double angle_tolerance)
{
  MspglibDataset *dataset;
  Magnetic_Cell *cell;
  Cell *cell_nonmagnetic;
  MagneticDataContainer *container;

  dataset = NULL;
  cell = NULL;
  cell_nonmagnetic = NULL;
  container = NULL;

  int i;

  if ((dataset = init_magnetic_dataset()) == NULL) {
    goto not_found;
  }


  if ((cell_nonmagnetic = cel_alloc_cell(num_atom)) == NULL) {
    free(dataset);
    dataset = NULL;
    goto not_found;
  }

  cel_set_cell(cell_nonmagnetic, lattice, position, types);
  if (cel_any_overlap_with_same_type(cell_nonmagnetic, symprec)) {
      cel_free_cell(cell_nonmagnetic);
      magnetic_cel_free_cell(cell);
      cell = NULL;
      free(dataset);
      dataset = NULL;
      goto atoms_too_close;
  }
  
  cell = initialize_magnetic_cell_from_cell(cell_nonmagnetic);
  cel_free_cell(cell_nonmagnetic);
  cell_nonmagnetic = NULL;


  for (i = 0;i<num_atom;i++)
  {
    mat_copy_vector_d3(cell->magmom[i],magmom[i]);
  }

  if ((container = det_magnetic_determine_all(cell,
                                              hall_number,
                                              symprec,
                                              angle_tolerance))
      != NULL) {
    if (set_magnetic_dataset(dataset,
                             cell,
                             container)) {
      det_free_magneticcontainer(container);
      container = NULL;
      magnetic_cel_free_cell(cell);
      cell = NULL;
      goto found;
    }
    det_free_magneticcontainer(container);
    container = NULL;
  }

  magnetic_cel_free_cell(cell);
  cell = NULL;
  free(dataset);
  dataset = NULL;

 not_found:
  spglib_error_code = MSPGERR_SPACEGROUP_SEARCH_FAILED;
  return NULL;

 atoms_too_close:
  spglib_error_code = MSPGERR_ATOMS_TOO_CLOSE;
  return NULL;

 found:

  spglib_error_code = MSPGLIB_SUCCESS;
  return dataset;
}

static MspglibDataset * init_magnetic_dataset(void)
{
  int i, j;
  MspglibDataset *dataset;

  dataset = NULL;

  if ((dataset = (MspglibDataset*) malloc(sizeof(MspglibDataset))) == NULL) {
    warning_print("spglib: Memory could not be allocated.");
    return NULL;
  }

  dataset->spacegroup_number = 0;
  dataset->hall_number = 0;
  strcpy(dataset->international_symbol, "");
  strcpy(dataset->hall_symbol, "");
  strcpy(dataset->choice, "");
  dataset->origin_shift[0] = 0;
  dataset->origin_shift[1] = 0;
  dataset->origin_shift[2] = 0;
  dataset->n_atoms = 0;
  dataset->mapping_to_primitive = NULL;
  dataset->n_operations = 0;
  dataset->rotations = NULL;
  dataset->translations = NULL;
  dataset->magmom = NULL;
  dataset->primitive_types = NULL;
  dataset->time_reversal = NULL;
  dataset->msgtype = 0;
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      dataset->transformation_matrix[i][j] = 0;
      dataset->primitive_lattice[i][j] = 0.0;
    }
  }
  /* dataset->pointgroup_number = 0; */
  strcpy(dataset->pointgroup_symbol, "");

  return dataset;

}

/* Return 0 if failed */
static int set_magnetic_dataset(SPGCONST MspglibDataset * dataset,
                                SPGCONST Magnetic_Cell * cell,
                                SPGCONST MagneticDataContainer * container)
{  
  int i;
  double inv_lat[3][3];
  Pointgroup pointgroup;

  /* Spacegroup type, transformation matrix, origin shift */
  dataset->msgtype = container->msg_type;
  dataset->n_atoms = container->primitive->size;
  dataset->n_std_atoms = container->primitive->cell->size;
  dataset->spacegroup_number = container->spacegroup->number;
  dataset->hall_number = container->spacegroup->hall_number;
  memcpy(dataset->international_symbol, container->spacegroup->international_short, 11);
  memcpy(dataset->hall_symbol, container->spacegroup->hall_symbol, 17);
  memcpy(dataset->choice, container->spacegroup->choice, 6);
  mat_inverse_matrix_d3(inv_lat, container->spacegroup->bravais_lattice, 0);
  mat_multiply_matrix_d3(dataset->transformation_matrix,
                         inv_lat, cell->lattice);
  mat_copy_vector_d3(dataset->origin_shift, container->spacegroup->origin_shift);
  dataset->n_operations = container->symmetry->size;

  if ((dataset->rotations =
       (int (*)[3][3]) malloc(sizeof(int[3][3]) * dataset->n_operations))
      == NULL) {
    warning_print("spglib: Memory could not be allocated.");
    goto err;
  }
  

  if ((dataset->translations =
       (double (*)[3]) malloc(sizeof(double[3]) * dataset->n_operations))
      == NULL) {
    warning_print("spglib: Memory could not be allocated.");
    goto err;
  }

  if ((dataset->time_reversal =
       (int *) malloc(sizeof(int) * dataset->n_operations))
      == NULL) {
    warning_print("spglib: Memory could not be allocated.");
    goto err;
  }

  for (i = 0; i < container->symmetry->size; i++) {
    mat_copy_matrix_i3(dataset->rotations[i], container->symmetry->rot[i]);
    mat_copy_vector_d3(dataset->translations[i], container->symmetry->trans[i]);
    dataset->time_reversal[i] = container->symmetry->time_reversal[i];
  }

  if ((dataset->mapping_to_primitive =
       (int*) malloc(sizeof(int) * dataset->n_atoms)) == NULL) {
    warning_print("spglib: Memory could not be allocated.");
    goto err;
  }

  debug_print("[line %d, %s]\n", __LINE__, __FILE__);
  debug_print("refined cell after ref_get_Wyckoff_positions\n");
  debug_print_matrix_d3(primitive->cell->lattice);
  #ifdef SPGDEBUG
  for (i = 0; i < container->primitive->cell->size; i++) {
    printf("%d: %f %f %f\n",
           container->primitive->cell->types[i],
           container->primitive->cell->position[i][0],
           container->primitive->cell->position[i][1],
           container->primitive->cell->position[i][2]);
  }
#endif

  mat_copy_matrix_d3(dataset->primitive_lattice, container->primitive->cell->lattice);
  for (i = 0; i < dataset->n_atoms; i++) {
    dataset->mapping_to_primitive[i] = container->primitive->mapping_table[i];
  }


  if ((dataset->positions =
       (double (*)[3]) malloc(sizeof(double[3]) * container->primitive->cell->size))
      == NULL) {
    warning_print("spglib: Memory could not be allocated.");
    goto err;
  }

  if ((dataset->magmom =
       (double (*)[3]) malloc(sizeof(double[3]) * container->primitive->cell->size))
      == NULL) {
    warning_print("spglib: Memory could not be allocated.");
    goto err;
  }

  if ((dataset->primitive_types =
       (int*) malloc(sizeof(int) * container->primitive->cell->size)) == NULL) {
    warning_print("spglib: Memory could not be allocated.");
    goto err;
  }

  for (i = 0; i < container->primitive->cell->size; i++) {
    mat_copy_vector_d3(dataset->positions[i], container->primitive->cell->position[i]);
    mat_copy_vector_d3(dataset->magmom[i], container->primitive->cell->magmom[i]);
    dataset->primitive_types[i] = container->primitive->cell->types[i];
  }
  /* dataset->pointgroup_number = spacegroup->pointgroup_number; */
  pointgroup = ptg_get_pointgroup(container->spacegroup->pointgroup_number);
  memcpy(dataset->pointgroup_symbol, pointgroup.symbol, 6);

  return 1;

 err:
    if (dataset->positions != NULL) {
      free(dataset->positions);
      dataset->positions = NULL;
    }
    if (dataset->primitive_types != NULL) {
      free(dataset->primitive_types);
      dataset->primitive_types = NULL;
    }
    if (dataset->magmom != NULL) {
      free(dataset->magmom);
      dataset->magmom = NULL;
    }
    if (dataset->mapping_to_primitive != NULL) {
      free(dataset->mapping_to_primitive);
      dataset->mapping_to_primitive = NULL;
    }
    if (dataset->translations != NULL) {
      free(dataset->translations);
      dataset->translations = NULL;
    }
    if (dataset->rotations != NULL) {
      free(dataset->rotations);
      dataset->rotations = NULL;
    }
    if (dataset->time_reversal != NULL) {
      free(dataset->time_reversal);
      dataset->time_reversal = NULL;
    }



  return 0;

}

void spg_free_magneticdataset(MspglibDataset *dataset)
{
  if (dataset->n_operations > 0) {
    free(dataset->rotations);
    dataset->rotations = NULL;
    free(dataset->translations);
    dataset->translations = NULL;
    free(dataset->time_reversal);
    dataset->n_operations = 0;
  }

  if (dataset->n_atoms > 0) {
    free(dataset->mapping_to_primitive);
    dataset->mapping_to_primitive = NULL;
    dataset->n_atoms = 0;
  }

  if (dataset->n_std_atoms > 0) {
    free(dataset->positions);
    dataset->positions = NULL;
    free(dataset->primitive_types);
    dataset->primitive_types = NULL;
    free(dataset->magmom);
    dataset->magmom = NULL;
    dataset->n_std_atoms = 0;
  }

  dataset->spacegroup_number = 0;
  dataset->hall_number = 0;
  strcpy(dataset->international_symbol, "");
  strcpy(dataset->hall_symbol, "");
  strcpy(dataset->choice, "");
  strcpy(dataset->pointgroup_symbol, "");

  free(dataset);
}

/*-------*/
/* error */
/*-------*/
MspglibError mspg_get_error_code(void)
{
  return spglib_error_code;
}

char * mspg_get_error_message(MspglibError error)
{
  int i;

  for (i = 0; i < 100; i++) {
    if (MSPGERR_NONE == mspglib_error_message[i].error) {
      break;
    }

    if (error == mspglib_error_message[i].error) {
      return mspglib_error_message[i].message;
    }
  }

  return NULL;
}




