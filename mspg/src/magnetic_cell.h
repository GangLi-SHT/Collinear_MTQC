#ifndef __magnetic_cell_H__
#define __magnetic_cell_H__

#include "mathfunc.h"
#include "cell.h"

typedef struct {
  double (*magmom)[3];
  int size;
  double (*lattice)[3]; /* 3x3 matrix */
  int *types;
  double (*position)[3];
} Magnetic_Cell;

Magnetic_Cell * magnetic_cel_alloc_cell(const int size);
void magnetic_cel_free_cell(Magnetic_Cell * cell);
Cell * Magnetic_Cell_to_Cell(const Magnetic_Cell * cell);
Magnetic_Cell * initialize_magnetic_cell_from_cell(const Cell * cell);


#endif
