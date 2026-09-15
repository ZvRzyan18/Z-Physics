#ifndef ZP_SOLVER2D_H
#define ZP_SOLVER2D_H

#include "z_physics/hint.h"
#include "z_physics/core2d/manifold2d.h"

typedef struct {
 float   _dt;
 float   _inv_dt;
  /* soft constraint values */
 float   _bias_coeff;
 float   _mass_coeff;
 float   _impulse_coeff;

 uint8_t _solver_substeps;
} zp_solver_input2d;

ZP_CPP_BEGIN

void zp_manifold2d_soft_prepare_contact(zp_manifold2d *const zp_restrict m, void *const zp_restrict w, const zp_solver_input2d *const input);
void zp_manifold2d_soft_presolve_contact(zp_manifold2d *const zp_restrict m, void *const zp_restrict w, const zp_solver_input2d *const input);
void zp_manifold2d_soft_solve_contact(zp_manifold2d *const zp_restrict m, void *const zp_restrict w, const zp_solver_input2d *const input);

ZP_CPP_END

#endif

