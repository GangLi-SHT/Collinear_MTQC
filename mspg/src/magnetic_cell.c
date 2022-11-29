#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "cell.h"
#include "magnetic_cell.h"
#include "mathfunc.h"

#include "debug.h"

#ifdef DMALLOC
#include "dmalloc.h"
#endif

#define INCREASE_RATE 2.0
#define REDUCE_RATE 0.95
#define NUM_ATTEMPT 100

/* NULL is returned if faied */
Magnetic_Cell * magnetic_cel_alloc_cell(const int size)
{
  Magnetic_Cell *cell;

  cell = NULL;

  if (size < 1) {
    return NULL;
  }

  cell = NULL;

  if ((cell = (Magnetic_Cell*) malloc(sizeof(Magnetic_Cell))) == NULL) {
    warning_print("spglib: Memory could not be allocated.");
    return NULL;
  }

  if ((cell->lattice = (double (*)[3]) malloc(sizeof(double[3]) * 3)) == NULL) {
    warning_print("spglib: Memory could not be allocated.");
    free(cell);
    cell = NULL;
    return NULL;
  }

  cell->size = size;

  if ((cell->types = (int *) malloc(sizeof(int) * size)) == NULL) {
    warning_print("spglib: Memory could not be allocated.");
    free(cell->lattice);
    cell->lattice = NULL;
    free(cell);
    cell = NULL;
    return NULL;
  }
  if ((cell->position =
       (double (*)[3]) malloc(sizeof(double[3]) * size)) == NULL) {
    warning_print("spglib: Memory could not be allocated.");
    free(cell->types);
    cell->types = NULL;
    free(cell->lattice);
    cell->lattice = NULL;
    free(cell);
    cell = NULL;
    return NULL;
  }

  if ((cell->magmom =
       (double (*)[3]) malloc(sizeof(double[3]) * size)) == NULL) {
    warning_print("spglib: Memory could not be allocated.");
    free(cell->types);
    cell->types = NULL;
    free(cell->lattice);
    cell->lattice = NULL;
    free(cell->position);
    cell->position = NULL;
    free(cell);
    cell = NULL;
    return NULL;
  }

  return cell;
}

void magnetic_cel_free_cell(Magnetic_Cell * cell)
{
  if (cell != NULL) {
    if (cell->lattice != NULL) {
      free(cell->lattice);
      cell->lattice = NULL;
    }
    if (cell->position != NULL) {
      free(cell->position);
      cell->position = NULL;
    }
    if (cell->types != NULL) {
      free(cell->types);
      cell->types = NULL;
    }
    if (cell->magmom != NULL) {
      free(cell->magmom);
      cell->magmom = NULL;
    }
    free(cell);
  }
}

Cell * Magnetic_Cell_to_Cell(const Magnetic_Cell * cell)
{
	Cell *cell_new;
	int i;

	if ((cell_new= cel_alloc_cell(cell->size)) == NULL) {
    return NULL;
    }
    
    mat_copy_matrix_d3(cell_new->lattice,cell->lattice);

    for (i=0;i<cell->size;i++)
    {
    	mat_copy_vector_d3(cell_new->position[i],cell->position[i]);
    	cell_new->types[i] = cell->types[i];
    }
    return cell_new;
}

Magnetic_Cell * initialize_magnetic_cell_from_cell(const Cell * cell)
{
	Magnetic_Cell *cell_new;
	double zero_magmom[3] = {0.0,0.0,0.0};
	int i;

	if ((cell_new = magnetic_cel_alloc_cell(cell->size)) == NULL) {
    return NULL;
    }
    
    mat_copy_matrix_d3(cell_new->lattice,cell->lattice);

    for (i=0;i<cell->size;i++)
    {
    	mat_copy_vector_d3(cell_new->position[i],cell->position[i]);
    	cell_new->types[i] = cell->types[i];
    }
    for (i=0;i<cell->size;i++)
    {
    	mat_copy_vector_d3(cell_new->magmom[i],zero_magmom);
    }
    
    return cell_new;

}
