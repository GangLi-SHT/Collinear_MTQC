#ifndef __mspacegroup_H__
#define __mspacegroup_H__

#include "cell.h"
#include "magnetic_primitive.h"
#include "mathfunc.h"
#include "primitive.h"
#include "magnetic_primitive.h"
#include "spacegroup.h"
#include "magnetic_cell.h"
#include "symmetry.h"

Magnetic_Symmetry * sym_get_magnetic_operation(const Magnetic_Cell * primitive,
                                               const double symprec,
                                               const double angle_tolerance);
Spacegroup * spa_search_mspacegroup(const Magnetic_Primitive * primitive,
                                    const int hall_number,
                                    const double symprec,
                                    const double angle_tolerance);
int magnetic_type_check( const Magnetic_Symmetry * symmetry, const double symprec);
Magnetic_Primitive * standard_primitive_cell(const Magnetic_Primitive * primitive, SPGCONST Spacegroup * spacegroup);

#endif
