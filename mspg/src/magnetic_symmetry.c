#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "cell.h"
#include "mathfunc.h"
#include "symmetry.h"
#include "magnetic_symmetry.h"
#include "overlap.h"
#include "debug.h"

#ifdef DMALLOC
#include "dmalloc.h"
#endif

#define NUM_ATOMS_CRITERION_FOR_OPENMP 1000
#define ANGLE_REDUCE_RATE 0.95
#define NUM_ATTEMPT 100
#define PI 3.14159265358979323846

/* Return NULL if failed */
Magnetic_Symmetry * sym_alloc_magnetic_symmetry(const int size)
{
  Magnetic_Symmetry *symmetry;

  symmetry = NULL;

  if (size < 1) {
    return NULL;
  }

  if ((symmetry = (Magnetic_Symmetry*) malloc(sizeof(Magnetic_Symmetry))) == NULL) {
    warning_print("spglib: Memory could not be allocated ");
    return NULL;
  }

  symmetry->size = size;
  symmetry->time_reversal = NULL;
  symmetry->rot = NULL;
  symmetry->trans = NULL;

  if ((symmetry->time_reversal =
       (int *) malloc(sizeof(int)*size)) == NULL) {
    warning_print("spglib: Memory could not be allocated ");
    warning_print("(line %d, %s).\n", __LINE__, __FILE__);
    free(symmetry);
    symmetry = NULL;
    return NULL;
  }

  if ((symmetry->rot =
       (int (*)[3][3]) malloc(sizeof(int[3][3]) * size)) == NULL) {
    warning_print("spglib: Memory could not be allocated ");
    warning_print("(line %d, %s).\n", __LINE__, __FILE__);
    free(symmetry->time_reversal);
    symmetry->time_reversal = NULL;
    free(symmetry);
    symmetry = NULL;
    return NULL;
  }
  if ((symmetry->trans =
       (double (*)[3]) malloc(sizeof(double[3]) * size)) == NULL) {
    warning_print("spglib: Memory could not be allocated ");
    warning_print("(line %d, %s).\n", __LINE__, __FILE__);
    free(symmetry->time_reversal);
    symmetry->time_reversal = NULL;
    free(symmetry->rot);
    symmetry->rot = NULL;
    free(symmetry);
    symmetry = NULL;
    return NULL;
  }

  return symmetry;
}

void sym_free_magnetic_symmetry(Magnetic_Symmetry *symmetry)
{
    if (symmetry->size > 0) {

    free(symmetry->rot);
    symmetry->rot = NULL;
    free(symmetry->trans);
    symmetry->trans = NULL;
    free(symmetry->time_reversal);
    
    symmetry->time_reversal=NULL;
  }
  free(symmetry);
}

Symmetry * magnetic_symmetry_to_symmetry(const Magnetic_Symmetry * symmetry)
{  int i;
   Symmetry *symmetry_new;

   if (symmetry!=NULL){
    symmetry_new = sym_alloc_symmetry(symmetry->size);
    for (i=0; i<symmetry->size ; i++)
    {
      mat_copy_matrix_i3(symmetry_new->rot[i],symmetry->rot[i]);
      mat_copy_vector_d3(symmetry_new->trans[i],symmetry->trans[i]);
    }
    return symmetry_new;
   }

   return NULL;
}

Magnetic_Symmetry * get_operation(void);
