#include "zp_physics/core2d/collision2d.h"
#include "zp_physics/core2d.h"
#include <assert.h>

typedef uint8_t (*zp_collider_function)(zp_manifold2d *const zp_restrict out, const zp_body2d *const zp_restrict a, const zp_body2d *const zp_restrict b);

#define BOX_BOX (zp_collider_function)zp_manifold2d_box_vs_box
/*
 lookup table of function pointer is much more easier to scale
 than manual switch.
*/
static zp_collider_function function_luts[3][3] = {
 {NULL, NULL, NULL},
 {NULL, BOX_BOX, NULL},
 {NULL, NULL, NULL},
};


/*
 collision dynamic dispatch
*/
uint8_t zp_collision2d_collide(zp_manifold2d *const zp_restrict out, const zp_body2d *zp_restrict a, const zp_body2d *zp_restrict b) {
 if(((a->_head._flags & ZP_BODY_MOVEMENT_MASK_2D) == ZP_BODY_MOVEMENT_STATIC_2D) && ((b->_head._flags & ZP_BODY_MOVEMENT_MASK_2D) == ZP_BODY_MOVEMENT_STATIC_2D)) 
  return 0;
 
 const zp_vec2 *const ba = a->_head._aabb;
 const zp_vec2 *const bb = b->_head._aabb;
 /* aabb checks */
 if((ba[1].x < bb[0].x || ba[0].x > bb[1].x || ba[1].y < bb[0].y || ba[0].y > bb[1].y))
  return 0;
 return function_luts[a->_head._flags & ZP_BODY_MASK_2D][b->_head._flags & ZP_BODY_MASK_2D](out, a, b);

}

