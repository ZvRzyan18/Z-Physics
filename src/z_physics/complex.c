#include "z_physics/complex.h"

zp_complex zp_complex_identity = {.x = 1.0f, .y = 0.0f};
zp_quaternion zp_quaternion_identity = {.x = 0.0f, .y = 0.0f, .z = 0.0f, .w = 1.0f};


/*
 this function requires a little bit of accuracy
 that is why exact values must be handled.
 PS. zp_sincos is not accurate enough to handle initial rotation
 because it is designed for speed.
*/
zp_complex zp_crotate(const float radians) {
 zp_complex out;

 const float compare_epsilon = 1e-2f;
 float pos_radians = zp_abs(radians);
 
 float sine_sign = zp_copysign(1.0f, radians); /* needs to perform eor sign bit */
 /* handle exact values */
 if(pos_radians < compare_epsilon) {
  out.x = 1.0f;
  out.y = 0.0f * sine_sign;
 } else if(zp_abs(pos_radians - 1.570796f) < compare_epsilon) {
  out.x = 0.0f;
  out.y = 1.0f * sine_sign;
 } else if(zp_abs(pos_radians - 3.14159f) < compare_epsilon) {
  out.x = -1.0f;
  out.y = 0.0f * sine_sign;
 } else if(zp_abs(pos_radians - 4.71238898f) < compare_epsilon) {
  out.x = 0.0f;
  out.y = -1.0f * sine_sign;
 } else if(zp_abs(pos_radians - 6.28318531f) < compare_epsilon) {
  out.x = 1.0f;
  out.y = 0.0f * sine_sign;
 }
 else /* fallback */
  zp_sincos(radians, &out.y, &out.x);
 return out;
}


float zp_cangle(const zp_complex c) {
 /* do not require much accuracy */
 return zp_atan2(c.y, c.x);
}


zp_quaternion zp_qfromeulerangle(float roll, const float pitch, const float yaw) {
 zp_quaternion out;
 float cy, sy;
 float cp, sp;
 float cr, sr;

 zp_sincos(yaw * 0.5f, &sy, &cy);
 zp_sincos(pitch * 0.5f, &sp, &cp);
 zp_sincos(roll * 0.5f, &sr, &cr);
 
 const float crcp = cr * cp;
 const float crsp = cr * sp;
 const float srcp = sr * cp;
 const float srsp = sr * sp;
 
 out.x = zp_fma(srcp, cy, -(crsp * sy));
 out.y = zp_fma(crsp, cy, srcp * sy);
 out.z = zp_fma(crcp, sy, -(srsp * cy));
 out.w = zp_fma(crcp, cy, srsp * sy);
 return out;
}


zp_vec3 zp_qtoeulerangle(const zp_quaternion a) {
 zp_vec3 out;
 const float y2 = a.y * a.y;
 float mpitch = 2.0f * zp_fma(a.w, a.y, -(a.z * a.x));
 /* asin is already clamped and will never return inf/nan */
 out.y = zp_asin(mpitch);
 /* TODO : use vector instruction in atan2 */
 out.x = zp_atan2(2.0f * zp_fma(a.w, a.x, a.y * a.z), 1.0f - 2.0f * (a.x * a.x + y2));
 out.z = zp_atan2(2.0f * zp_fma(a.w, a.z, a.x * a.y), 1.0f - 2.0f * (y2 + a.z * a.z));
 return out;
}


zp_quaternion zp_qfromaxisangle(const float a, const float x, const float y, const float z) {
 zp_quaternion out;
 float s;
 zp_sincos(a * 0.5f, &s, &out.w);
	out.x = x * s;
 out.y = y * s;
	out.z = z * s;
	return out;
}

zp_vec4 zp_qtoaxisangle(const zp_quaternion a) {
 zp_vec4 out;
 out.w = zp_asin(a.w);
 out.w += out.w;
 const float b = zp_rsqrt(1.0f - a.w * a.w);
 out.x = a.x * b;
 out.y = a.y * b;
 out.z = a.z * b;
 return out;
}




