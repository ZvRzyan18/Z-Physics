#ifndef ZP_BODY2D_H
#define ZP_BODY2D_H

#include "z_physics/hint.h"
#include "z_physics/type.h"
#include "z_physics/complex.h"


/*
 create infos
*/
typedef struct {
 uint16_t _flags; /* TODO : place it at the bottom and then use offset of?*/
 zp_vec2  _position;
 zp_vec2  _velocity;
 zp_vec2  _force;
 zp_vec2  _half_size;
 float    _linear_damping;
 float    _rotation;
 float    _omega;
 float    _torque;
 float    _density;
 float    _angular_damping;
 float    _restitution;
 float    _friction;
} zp_create_body2d;


/*
 bodies
*/

typedef struct {
 uint16_t         _flags;
 zp_vec2          _aabb[2];
 zp_vec2          _position;
 zp_vec2          _velocity;
 zp_vec2          _force;
 zp_complex       _rotation;
 float            _linear_damping;
 float            _omega;
 float            _torque;
 float            _angular_damping;
 float            _inv_inertia;
 float            _inv_mass;
 float            _restitution;
 float            _friction;
 float            _idle_time;
 
 zp_container_id  _id; /* minmize the byte size */
} zp_head2d;

typedef struct {
 zp_head2d _head;
 zp_vec2   _half_size;
} zp_box2d;

typedef union {
 zp_head2d _head;
 zp_box2d  _box;
} zp_body2d;

ZP_CPP_BEGIN 

zp_cold zp_noinline void zp_body2d_init(zp_body2d *const zp_restrict body, const void *const zp_restrict data);
void zp_body2d_updatev(zp_body2d *const zp_restrict body, const void *const zp_restrict world, const float dt);
void zp_body2d_updatep(zp_body2d *const zp_restrict body, const void *const zp_restrict world, const float dt);

ZP_CPP_END

#endif


