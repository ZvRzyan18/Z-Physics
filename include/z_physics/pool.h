#ifndef ZP_POOL_H
#define ZP_POOL_H

#include "z_physics/hint.h"
#include "z_physics/type.h"
#include <assert.h>

typedef struct {
 uint8_t    *_bytes;
 zp_pool_id *_free_list;
 zp_pool_id  _size;
 zp_pool_id  _reserve;
 zp_pool_id  _max_size;
 float       _growth_base;
 uint16_t    _stride;
 uint16_t    _increase_count;
} zp_pool;

ZP_CPP_BEGIN

zp_cold int zp_pool_init(zp_pool *const zp_restrict pool, const uint16_t stride, const uint16_t reserve, const float growth_base);
zp_cold void zp_pool_destroy(zp_pool *const zp_restrict pool);
zp_pool_id zp_pool_acquire(zp_pool *const zp_restrict pool);
void zp_pool_release(zp_pool *const zp_restrict pool, const zp_pool_id id);
zp_inline void* zp_pool_get(zp_pool *const zp_restrict pool, zp_pool_id id) {
 assert(id != ZP_POOL_NULL_ID);
 return pool->_bytes + ((size_t)id) * ((size_t)pool->_stride);
}
ZP_CPP_END

#endif

