#ifndef ZP_CONTAINER_H
#define ZP_CONTAINER_H

#include "z_physics/hint.h"
#include "z_physics/type.h"
#include <assert.h>

typedef struct {
 uint8_t  *_bytes;
 zp_container_id *_to_index_lut;
 zp_container_id *_to_id_lut;
 zp_container_id *_free_list;
 float           _growth_base;
 zp_container_id _stride;
 zp_container_id _reserve;
 zp_container_id _max_size;
 zp_container_id _size;
 uint16_t        _increase_count;
} zp_container;


ZP_CPP_BEGIN

zp_cold int zp_container_init(zp_container *const zp_restrict c, const zp_container_id stride, const zp_container_id reserve, const float growth_base);
zp_cold void zp_container_destroy(zp_container *const zp_restrict c);
zp_container_id zp_container_acquire(zp_container *const zp_restrict c);
void zp_container_release(zp_container *const zp_restrict c, zp_container_id id);

zp_hot zp_inline void* zp_container_get(zp_container *const zp_restrict c, const zp_container_id id) {
 assert(id != ZP_CONTAINER_NULL_ID);
 assert(c->_to_index_lut[id] < c->_size);
 return (void*)(c->_bytes + ((size_t)c->_to_index_lut[id] * (size_t)c->_stride));
}

ZP_CPP_END

#endif

