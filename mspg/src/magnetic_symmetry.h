/* Extension to spglib to include magnetic space group */
/* 2021/2/7 */
/* Author: Su Yunlong*/

#ifndef __magnetic_symmetry_H__
#define __magnetic_symmetry_H__

#include "cell.h"
#include "mathfunc.h"
#include "symmetry.h"

typedef struct {
  int size;
  int *time_reversal;
  int (*rot)[3][3];
  double (*trans)[3];
} Magnetic_Symmetry;

Magnetic_Symmetry * sym_alloc_magnetic_symmetry(const int size);
void sym_free_magnetic_symmetry(Magnetic_Symmetry * symmetry);

Symmetry * magnetic_symmetry_to_symmetry(const Magnetic_Symmetry * symmetry);

#endif
