#ifndef ZP_TYPES_H
#define ZP_TYPES_H

#include "z_physics/vector.h"
#include <stdint.h>
#include <stddef.h>

typedef size_t   zp_bump_ptr;
typedef uint32_t zp_handle;
typedef uint32_t zp_pool_id;
#define ZP_POOL_NULL_ID 0xFFFFFFFF

typedef uint16_t zp_container_id;
#define ZP_CONTAINER_NULL_ID 0xFFFF

typedef zp_vec2 zp_complex;
typedef zp_vec4 zp_quaternion;

#endif

