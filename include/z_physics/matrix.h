#ifndef ZP_MATRIX_H
#define ZP_MATRIX_H

#include "z_physics/vector.h"


typedef union {
 zp_vec2 basis[2];
 float v[2][2];
 float arr[4];
} zp_mat2x2;

typedef union {
 zp_vec3 basis[3];
 float v[3][3];
 float arr[9];
} zp_mat3x3;

typedef union {
 zp_vec4 basis[4];
 float v[4][4];
 float arr[16];
} zp_mat4x4;


ZP_CPP_BEGIN

extern zp_mat2x2 zp_mat2x2_identity;
extern zp_mat3x3 zp_mat3x3_identity;
extern zp_mat4x4 zp_mat4x4_identity;


zp_mat3x3 zp_fromquaternion3(const zp_vec4 a);
zp_mat4x4 zp_fromquaternion4(const zp_vec4 a);

zp_mat2x2 zp_transpose2(const zp_mat2x2 a);
zp_mat3x3 zp_transpose3(const zp_mat3x3 a);
zp_mat4x4 zp_transpose4(const zp_mat4x4 a);

zp_mat2x2 zp_mulmm2(const zp_mat2x2 a, const zp_mat2x2 b);
zp_mat3x3 zp_mulmm3(const zp_mat3x3 a, const zp_mat3x3 b);
zp_mat4x4 zp_mulmm4(const zp_mat4x4 a, const zp_mat4x4 b);

zp_vec2 zp_mulm2v2(const zp_mat2x2 a, const zp_vec2 b);
zp_vec3 zp_mulm3v3(const zp_mat3x3 a, const zp_vec3 b);
zp_vec2 zp_mulm3v2(const zp_mat3x3 a, const zp_vec2 b);
zp_vec4 zp_mulm4v4(const zp_mat4x4 a, const zp_vec4 b);
zp_vec3 zp_mulm4v3(const zp_mat4x4 a, const zp_vec3 b);



ZP_CPP_END


#endif



