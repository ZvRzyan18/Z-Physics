#ifndef ZP_WORLD2D_H
#define ZP_WORLD2D_H

#include "z_physics/hint.h"
#include "z_physics/type.h"
#include "z_physics/container.h"
#include "z_physics/core2d/contacthash2d.h"

typedef struct {
 zp_vec2  _gravity;
 float    _growth_base_rate;
 uint16_t _bodies_initial_reserve;
 uint16_t _joints_initial_reserve;
 uint16_t _contacts_initial_reserve;
 uint8_t  _hash_max_buckets;
 uint8_t  _solver_substeps;
 uint8_t  _time_substeps;
} zp_create_world2d;

typedef struct {
 /* NOTE : this can change anytime, breaking compatibility, so do not access */
 zp_container     _body_container;
 zp_contacthash2d _contacts;
 zp_vec2          _gravity;
 float            _inv_timestep_substep;
 float            _growth_base;
 uint16_t         _solver_substeps;
 uint16_t         _time_substeps;
} zp_world2d;

typedef struct {
 zp_vec2    _position;
 zp_complex _rotation;
 zp_vec2    _velocity;
 float      _omega;
 uint8_t    _is_sleeping;
} zp_bodydata2d;

ZP_CPP_BEGIN

zp_cold int zp_world2d_init(zp_world2d *const zp_restrict world, const zp_create_world2d *const zp_restrict data);
zp_cold void zp_world2d_destroy(zp_world2d *const zp_restrict world);
zp_hot void zp_world2d_update(zp_world2d *const zp_restrict world, const float dt);
zp_handle zp_world2d_create_body(zp_world2d *const zp_restrict world, const void *const zp_restrict value);
void zp_world2d_remove_body(zp_world2d *const zp_restrict world, const zp_handle id);
zp_hot void zp_world2d_get_bodydata(zp_world2d *const zp_restrict world, const zp_handle id, zp_bodydata2d *const zp_restrict dat);


ZP_CPP_END

#endif

