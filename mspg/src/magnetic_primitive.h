#ifndef __magnetic_primitive_H__
#define __magnetic_primitive_H__

#include "cell.h"
#include "mathfunc.h"
#include "symmetry.h"
#include "primitive.h"
#include "magnetic_symmetry.h"
#include "magnetic_cell.h"



typedef struct {
  Magnetic_Cell *cell;
  int *mapping_table;
  int size;
  double tolerance;
  double angle_tolerance;
  double (*orig_lattice)[3]; /* 3x3 matrix */
} Magnetic_Primitive;

Magnetic_Primitive * prm_get_magnetic_primitive(const Magnetic_Cell * cell,
                                       const double symprec,
                                       const double angle_tolerance);
Magnetic_Symmetry * prm_get_magnetic_primitive_symmetry(const Magnetic_Symmetry *symmetry,
                                       const double symprec);
Cell * relabeling_types_according_magmom(const Magnetic_Cell *cell, const double symprec);

Primitive * Magnetic_Primitive_to_Primitive(const Magnetic_Primitive * primitive);
Magnetic_Primitive * magn_prm_alloc_primitive(const int size);
void magn_prm_free_primitive(Magnetic_Primitive * primitive);

#endif