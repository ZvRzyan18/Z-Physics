#ifndef ZP_COMPLEX_H
#define ZP_COMPLEX_H

#include "zp_physics/types.h"
#include "zp_physics/math.h"

ZP_CPP_BEGIN

extern zp_complex zp_complex_identity;
extern zp_quaternion zp_quaternion_identity;


/*
 complex for 2d rotation
*/
zp_const zp_inline float zp_creal(const zp_complex c) {
 return c.x;
}

zp_const zp_inline float zp_cimag(const zp_complex c) {
 return c.y;
}

zp_const zp_inline zp_complex zp_cmul(const zp_complex a, const zp_complex b) {
 zp_complex out;
 out.x = zp_fma(zp_creal(a), zp_creal(b), -(zp_cimag(a) * zp_cimag(b)));
 out.y = zp_fma(zp_creal(a), zp_cimag(b), (zp_cimag(a) * zp_creal(b)));
 return out;
} 

zp_const zp_inline zp_complex zp_cconj(const zp_complex a) {
 zp_complex out;
 out.x = zp_creal(a);
 out.y = -zp_cimag(a);
 return out;
} 

zp_complex zp_crotate(const float radians);
float zp_cangle(const zp_complex c);




/*
 quaternion for 3d rotation
*/
zp_const zp_inline zp_quaternion zp_qmul(const zp_quaternion a, const zp_quaternion b) {
 zp_quaternion out;
/*
 out.x = a.x * b.w + a.y * b.z - a.z * b.y + a.w * b.x;
 out.y = -a.x * b.z + a.y * b.w + a.z * b.x + a.w * b.y;
 out.z = a.x * b.y - a.y * b.x + a.z * b.w + a.w * b.z;
 out.w = -a.x * b.x - a.y * b.y - a.z * b.z + a.w * b.w;
*/
 out.x = zp_fma(a.x, b.w, zp_fma(a.y, b.z, -zp_fma(a.z, b.y, a.w * b.x)));
 out.y = zp_fma(-a.x, b.z, zp_fma(a.y, b.w, zp_fma(a.z, b.x, a.w * b.y)));
 out.z = zp_fma(a.x, b.y, -zp_fma(a.y, b.x, zp_fma(a.z, b.w, a.w * b.z)));
 out.w = zp_fma(-a.x, b.x, -zp_fma(a.y, b.y, -zp_fma(a.z, b.z, a.w * b.w)));

 return out;
}


zp_const zp_inline zp_quaternion zp_qconj(const zp_quaternion a) {
 zp_quaternion out;
 out.x = -a.x;
 out.y = -a.y;
 out.z = -a.z;
 out.w = a.w;
 return out;
}

zp_quaternion zp_qfromaxisangle(const float a, const float x, const float y, const float z);
zp_vec4 zp_qtoaxisangle(const zp_quaternion a);


ZP_CPP_END

#endif

