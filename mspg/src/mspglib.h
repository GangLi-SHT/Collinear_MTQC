#ifndef __mspglib_H__
#define __mspglib_H__

#ifdef __cplusplus
extern "C" {
#endif
/* SPGCONST is used instead of 'const' so to avoid gcc warning. */
/* However there should be better way than this way.... */
#ifndef SPGCONST
#define SPGCONST
#endif

#include <stddef.h>

/*
  ------------------------------------------------------------------

  lattice: Lattice vectors (in Cartesian)

  [ [ a_x, b_x, c_x ],
  [ a_y, b_y, c_y ],
  [ a_z, b_z, c_z ] ]

  position: Atomic positions (in fractional coordinates)

  [ [ x1_a, x1_b, x1_c ],
  [ x2_a, x2_b, x2_c ],
  [ x3_a, x3_b, x3_c ],
  ...                   ]

  types: Atom types, i.e., species identified by number

  [ type_1, type_2, type_3, ... ]

  rotation: Rotation matricies of symmetry operations

  each rotation is:
  [ [ r_aa, r_ab, r_ac ],
  [ r_ba, r_bb, r_bc ],
  [ r_ca, r_cb, r_cc ] ]

  translation: Translation vectors of symmetry operations

  each translation is:
  [ t_a, t_b, t_c ]

  symprec: Distance tolerance in Cartesian coordinates to find crystal
           symmetry.

  ------------------------------------------------------------------

  Definitio of the operation:
  r : rotation     3x3 matrix
  t : translation  vector

  x_new = r * x + t:
  [ x_new_a ]   [ r_aa, r_ab, r_ac ]   [ x_a ]   [ t_a ]
  [ x_new_b ] = [ r_ba, r_bb, r_bc ] * [ x_b ] + [ t_b ]
  [ x_new_c ]   [ r_ca, r_cb, r_cc ]   [ x_c ]   [ t_c ]

  ------------------------------------------------------------------
*/

  typedef enum {
    MSPGLIB_SUCCESS = 0,
    MSPGERR_SPACEGROUP_SEARCH_FAILED,
    MSPGERR_CELL_STANDARDIZATION_FAILED,
    MSPGERR_SYMMETRY_OPERATION_SEARCH_FAILED,
    MSPGERR_ATOMS_TOO_CLOSE,
    MSPGERR_POINTGROUP_NOT_FOUND,
    MSPGERR_NIGGLI_FAILED,
    MSPGERR_DELAUNAY_FAILED,
    MSPGERR_ARRAY_SIZE_SHORTAGE,
    MSPGERR_NONE,
  } MspglibError;

  typedef struct {
  MspglibError error;
  char *message;
} MspglibErrorMessage;


  typedef struct {
    int spacegroup_number;
    int hall_number;
    char international_symbol[11];
    char hall_symbol[17];
    char choice[6];
    double transformation_matrix[3][3];
    double origin_shift[3];
    int n_operations;
    int (*rotations)[3][3];
    double (*translations)[3];
    int *time_reversal;
    int n_atoms;
    double primitive_lattice[3][3];
    int n_std_atoms;
    int *primitive_types;
    double (*positions)[3];
    double (*magmom)[3];
    int *mapping_to_primitive;
    /* int pointgroup_number; */
    char pointgroup_symbol[6];
    int msgtype;
  } MspglibDataset;
  MspglibDataset * mspg_get_magnetic_dataset(SPGCONST double lattice[3][3],
                                             SPGCONST double position[][3],
                                             SPGCONST double magmom[][3],
                                             const int types[],
                                             const int num_atom,
                                             const double symprec);
  void spg_free_magneticdataset(MspglibDataset *dataset);

  MspglibError mspg_get_error_code(void);
  char * mspg_get_error_message(MspglibError error);

#endif
