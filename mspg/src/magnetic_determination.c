#include <stdlib.h>
#include <stdio.h>
#include "cell.h"
#include "mspacegroup.h"
#include "magnetic_determination.h"
#include "primitive.h"
#include "magnetic_cell.h"
#include "magnetic_primitive.h"
#include "spacegroup.h"
#include "debug.h"

#ifdef DMALLOC
#include "dmalloc.h"
#endif

#define REDUCE_RATE_OUTER 0.9
#define NUM_ATTEMPT_OUTER 10
#define REDUCE_RATE 0.95
#define ANGLE_REDUCE_RATE 0.95
#define NUM_ATTEMPT 20


void det_free_magneticcontainer(MagneticDataContainer * container);


static MagneticDataContainer * get_magnetic_spacegroup_and_primitive(const Magnetic_Cell * cell,
                                                                    const int hall_number,
                                                                    const double symprec,
                                                                    const double angle_symprec)
{
  int attempt;
  double tolerance, angle_tolerance;
  Magnetic_Primitive *primitive;
  MagneticDataContainer *container;

  debug_print("get_spacegroup_and_primitive (tolerance = %f):\n", symprec);

  container = NULL;

 if ((container = (MagneticDataContainer*) malloc(sizeof(MagneticDataContainer))) == NULL) {
    warning_print("spglib: Memory could not be allocated.");
    return NULL;
  }
  
  primitive = NULL;
  container->primitive = NULL;
  container->spacegroup = NULL;
  container->symmetry = NULL;
  container->msg_type = 0; 

  tolerance = symprec;
  angle_tolerance = angle_symprec;
  
  for (attempt = 0; attempt < NUM_ATTEMPT; attempt++) {
    if ((primitive = prm_get_magnetic_primitive(cell,
                                                tolerance,
                                                angle_tolerance)) != NULL) {


      debug_print("[line %d, %s]\n", __LINE__, __FILE__);
      debug_print("primitive lattice\n");
      debug_print_matrix_d3(container->primitive->cell->lattice);

      container->primitive = primitive;


      if ((container->spacegroup = spa_search_mspacegroup(
             primitive,
             hall_number,
             container->primitive->tolerance,
             container->primitive->angle_tolerance)) != NULL) {
        goto found;
      }
      magn_prm_free_primitive(container->primitive);
      container->primitive = NULL;
    }

    warning_print("spglib: Attempt %d tolerance = %f failed.",
                  attempt, tolerance);
    warning_print(" (line %d, %s).\n", __LINE__, __FILE__);

    tolerance *= REDUCE_RATE;
    if (angle_tolerance > 0) {
      angle_tolerance *= ANGLE_REDUCE_RATE;
    }
  }
  magn_prm_free_primitive(primitive);
  primitive = NULL;

  det_free_magneticcontainer(container);
  container = NULL;

  return NULL;

found:
  container->primitive = NULL;
  container->primitive = standard_primitive_cell(primitive,container->spacegroup);
  container->symmetry = sym_get_magnetic_operation(container->primitive->cell, tolerance, angle_tolerance);
  container->msg_type = magnetic_type_check(container->symmetry,tolerance);
  magn_prm_free_primitive(primitive);
  primitive = NULL;

  return container;
}

void det_free_magneticcontainer(MagneticDataContainer * container)
{
  if (container != NULL) {
    if (container->spacegroup != NULL) {
      free(container->spacegroup);
      container->spacegroup = NULL;
    }
    if (container->primitive != NULL) {
      magn_prm_free_primitive(container->primitive);
      container->primitive = NULL;
    }
    if (container->symmetry != NULL) {
      sym_free_magnetic_symmetry(container->symmetry);
      container->primitive = NULL;
    }
    free(container);
  }
}

MagneticDataContainer * det_magnetic_determine_all(const Magnetic_Cell * cell,
                                                   const int hall_number,
                                                   const double symprec,
                                                   const double angle_symprec)
{
  int attempt;
  double tolerance;
  MagneticDataContainer *container;

  container = NULL;


  if (hall_number < 0 || hall_number > 530) {
    return NULL;
  }

  tolerance = symprec;
  for (attempt = 0; attempt < NUM_ATTEMPT_OUTER; attempt++) {
    if ((container = get_magnetic_spacegroup_and_primitive(cell,
                                                  hall_number,
                                                  tolerance,
                                                  angle_symprec)) != NULL) {
     
      goto found;
    }
    tolerance *= REDUCE_RATE_OUTER;
  }

found:
  return container;
}
