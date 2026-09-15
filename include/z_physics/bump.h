#ifndef ZP_BUMP_H
#define ZP_BUMP_H

#include "z_physics/hint.h"

typedef struct {
 uint8_t *_bytes;
 size_t _size;
 size_t _max_size;
} zp_bump;

ZP_CPP_BEGIN

zp_cold int zp_bump_init(zp_bump *zp_restrict b, size_t max_size);
zp_cold void zp_bump_destroy(zp_bump *zp_restrict b);
size_t zp_bump_acquire(zp_bump *zp_restrict b, uint16_t requested_bytes, uint8_t alignment);
void zp_bump_reset(zp_bump *zp_restrict b, size_t offset);
zp_inline void* zp_bump_get(zp_bump *const zp_restrict b, const size_t offset) {
 return b->_bytes + offset;
}

ZP_CPP_END

#endif

