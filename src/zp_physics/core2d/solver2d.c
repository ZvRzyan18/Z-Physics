#include "zp_physics/core2d/solver2d.h"
#include "zp_physics/core2d/world2d.h"
#include "zp_physics/core2d.h"

/*
 src : https://github.com/erincatto/solver2d
*/


void zp_manifold2d_soft_prepare_contact(zp_manifold2d *const zp_restrict m, void *const zp_restrict w, const zp_solver_input2d *const input) {
 zp_world2d *const world = (zp_world2d*)w;

 /* 
  checking for missing body, because if you remove a body, that would 
  results in a stale object id, but it is fine, because by adding these checks
  can remove a body as fast as O(1)
 */
 if(zp_unlikely(world->_body_container._to_index_lut[m->_body_a] >= world->_body_container._size))
  return;
 zp_body2d *const body_a = (zp_body2d*)zp_container_get(&world->_body_container, m->_body_a);
 if(zp_unlikely(world->_body_container._to_index_lut[m->_body_b] >= world->_body_container._size))
  return;
 zp_body2d *const body_b = (zp_body2d*)zp_container_get(&world->_body_container, m->_body_b);

 zp_compiler_memory_barrier();

 zp_complex a_rot = body_a->_head._rotation;
 zp_vec2 a_velocity = body_a->_head._velocity;
 float a_inv_mass = body_a->_head._inv_mass;
 float a_inv_inertia = body_a->_head._inv_inertia;
 float a_restitution = body_a->_head._restitution;

 zp_compiler_memory_barrier();

 zp_complex b_rot = body_b->_head._rotation;
 zp_vec2 b_velocity = body_b->_head._velocity;
 float b_inv_mass = body_b->_head._inv_mass;
 float b_inv_inertia = body_b->_head._inv_inertia;
 float b_restitution = body_b->_head._restitution;


 float restitution_threshold = zp_sqrt(zp_dot2(world->_gravity, world->_gravity)) * 0.4f * input->_dt;


 float inv_mass = a_inv_mass + b_inv_mass;
 float e = zp_min(a_restitution, b_restitution);

 for(uint8_t i = 0; i < m->_contact_count; i++) {
 	zp_contact2d *const contact = m->_contacts + i;

  zp_vec2 r1 = zp_cmul(a_rot, contact->_r1);
  zp_vec2 r2 = zp_cmul(b_rot, contact->_r2);

	 float rn1 = zp_dot2(r1, contact->_normal);
	 float rn2 = zp_dot2(r2, contact->_normal);
		float k_normal = inv_mass;
		k_normal += a_inv_inertia * (zp_dot2(r1, r1) - rn1 * rn1) + b_inv_inertia * (zp_dot2(r2, r2) - rn2 * rn2);
 	contact->_mass_normal = 1.0f / k_normal;

	 zp_vec2 tangent = zp_perp2(contact->_normal);
	 float rt1 = zp_dot2(r1, tangent);
	 float rt2 = zp_dot2(r2, tangent);
	 float k_tangent = inv_mass;
	 k_tangent += a_inv_inertia * (zp_dot2(r1, r1) - rt1 * rt1) + b_inv_inertia * (zp_dot2(r2, r2) - rt2 * rt2);
	 contact->_mass_tangent = 1.0f / k_tangent;

  /*
   Not correct. but by allowing only the linear relative velocity
   makes it more predictable.
  */
  zp_vec2 relative_vel = zp_sub2(b_velocity, a_velocity);
  float impact_speed = zp_dot2(relative_vel, contact->_normal);

  /* 
   hyperbolic tangent gives a smooth transition betweew low up to max resitution coeffs


   tanh(x)
   1.0 ---------------------------------
                      ____
                  ___
              ___
           __
         __
       _
     -
    -
  -------------------------------------

  conditional by threshold
  -------------------------------------
                  |-------------------
                  |
                  |
                  |
                  |
                  |
   _______________
  -------------------------------------
  
  */
  if(impact_speed < -0.1f) {
   float x = (-impact_speed - restitution_threshold);
   contact->_bias = (-impact_speed * e) * zp_tanh(zp_copysign(x, restitution_threshold));
  }
 }
}



void zp_manifold2d_soft_presolve_contact(zp_manifold2d *const zp_restrict m, void *const zp_restrict w, const zp_solver_input2d *const input) {
 (void)input;
 
 zp_world2d *const world = (zp_world2d*)w;
 
 if(zp_unlikely(world->_body_container._to_index_lut[m->_body_a] >= world->_body_container._size))
  return;
 zp_body2d *const body_a = (zp_body2d*)zp_container_get(&world->_body_container, m->_body_a);
  
 if(zp_unlikely(world->_body_container._to_index_lut[m->_body_b] >= world->_body_container._size))
  return;
 zp_body2d *const body_b = (zp_body2d*)zp_container_get(&world->_body_container, m->_body_b);

 zp_compiler_memory_barrier();

 zp_complex a_rot = body_a->_head._rotation;
 zp_vec2 a_velocity = body_a->_head._velocity;
 float a_omega = body_a->_head._omega;
 float a_inv_mass = body_a->_head._inv_mass;
 float a_inv_inertia = body_a->_head._inv_inertia;

 zp_compiler_memory_barrier();

 zp_complex b_rot = body_b->_head._rotation;
 zp_vec2 b_velocity = body_b->_head._velocity;
 float b_omega = body_b->_head._omega;
 float b_inv_mass = body_b->_head._inv_mass;
 float b_inv_inertia = body_b->_head._inv_inertia;

 for(uint8_t i = 0; i < m->_contact_count; i++) {
 	zp_contact2d *const contact = m->_contacts + i;

  zp_vec2 r1 = zp_cmul(a_rot, contact->_r1);
  zp_vec2 r2 = zp_cmul(b_rot, contact->_r2);
  
  float j = contact->_accumulated_normal;
	 zp_vec2 tangent = zp_perp2(contact->_normal);
	 
	 zp_vec2 impulse = zp_add2(zp_mul2(zp_stv2(j), contact->_normal), zp_mul2(zp_stv2(contact->_accumulated_tangent), tangent));
 	a_velocity = zp_sub2(a_velocity, zp_mul2(zp_stv2(a_inv_mass), impulse));
	 a_omega -= a_inv_inertia * zp_cross2(r1, impulse);
		 
	 b_velocity = zp_add2(b_velocity, zp_mul2(zp_stv2(b_inv_mass), impulse));
	 b_omega += b_inv_inertia * zp_cross2(r2, impulse);
 } 
 body_a->_head._velocity = a_velocity;
 body_a->_head._omega = a_omega;
 
 zp_compiler_memory_barrier();

 body_b->_head._velocity = b_velocity;
 body_b->_head._omega = b_omega;
}


void zp_manifold2d_soft_solve_contact(zp_manifold2d *const zp_restrict m, void *const zp_restrict w, const zp_solver_input2d *const input) {
 zp_world2d *const world = (zp_world2d*)w;
 
 if(zp_unlikely(world->_body_container._to_index_lut[m->_body_a] >= world->_body_container._size))
  return;
 zp_body2d *const body_a = (zp_body2d*)zp_container_get(&world->_body_container, m->_body_a);
  
 if(zp_unlikely(world->_body_container._to_index_lut[m->_body_b] >= world->_body_container._size))
  return;
 zp_body2d *const body_b = (zp_body2d*)zp_container_get(&world->_body_container, m->_body_b);

 zp_compiler_memory_barrier();

 zp_vec2 a_pos = body_a->_head._position;
 zp_complex a_rot = body_a->_head._rotation;
 zp_vec2 a_velocity = body_a->_head._velocity;
 float a_omega = body_a->_head._omega;
 float a_inv_mass = body_a->_head._inv_mass;
 float a_inv_inertia = body_a->_head._inv_inertia;
 float a_friction = body_a->_head._friction;

 zp_compiler_memory_barrier();

 zp_vec2 b_pos = body_b->_head._position;
 zp_complex b_rot = body_b->_head._rotation;
 zp_vec2 b_velocity = body_b->_head._velocity;
 float b_omega = body_b->_head._omega;
 float b_inv_mass = body_b->_head._inv_mass;
 float b_inv_inertia = body_b->_head._inv_inertia;
 float b_friction = body_b->_head._friction;


 float friction = zp_sqrt(a_friction * b_friction); 

 zp_vec2 r1, r2;
 zp_vec2 va, vb;
 zp_vec2 relative_vel;
 zp_vec2 impulse;
 float j, depth;


 for(uint8_t i = 0; i < m->_contact_count; i++) {
 	zp_contact2d *const contact = m->_contacts + i;

  r1 = contact->_r1;
  r2 = contact->_r2;


  r1 = zp_cmul(a_rot, r1);
  r2 = zp_cmul(b_rot, r2);

  float adjusted_depth = contact->_depth + zp_dot2(zp_sub2(r2, r1), contact->_normal);

  zp_vec2 p1 = zp_add2(a_pos, contact->_r1);
  zp_vec2 p2 = zp_add2(b_pos, contact->_r2);
  depth = -(zp_dot2(zp_sub2(p2, p1), contact->_normal) - adjusted_depth);
  

	 va = zp_add2(a_velocity, zp_cross_sv2(a_omega, r1));
  vb = zp_add2(b_velocity, zp_cross_sv2(b_omega, r2));
	 relative_vel = zp_sub2(vb, va);

		float vn = zp_dot2(relative_vel, contact->_normal);


  float slop = 0.004f;
  float penetration_error = 0.0f;
  
  if(depth < 0.0f) {
   float error = zp_min((contact->_depth + slop), 0.0f);
   penetration_error = -error * input->_bias_ratio;
  } else {
   penetration_error = contact->_depth * input->_inv_dt;
  }

  /* logarithmic with damping-like behaviour */
  penetration_error = zp_log2(penetration_error + 1.0f) * 0.69314718f;
  

  if(vn < contact->_bias) {
   penetration_error = -(vn - contact->_bias) + penetration_error;
  } else {
   penetration_error = -vn + penetration_error;
  }

	
  float mass = contact->_mass_normal * input->_mass_coeff;
	 j = mass * penetration_error;
  j -= input->_impulse_coeff * contact->_accumulated_normal;

  
		float accumulated_normal = contact->_accumulated_normal;
		contact->_accumulated_normal = zp_max(accumulated_normal + j, 0.0f);
		j = contact->_accumulated_normal - accumulated_normal;
  
	 impulse = zp_mul2(zp_stv2(j), contact->_normal);
	 
 	a_velocity = zp_sub2(a_velocity, zp_mul2(zp_stv2(a_inv_mass), impulse));
	 a_omega -= a_inv_inertia * zp_cross2(r1, impulse);
		 
	 b_velocity = zp_add2(b_velocity, zp_mul2(zp_stv2(b_inv_mass), impulse));
	 b_omega += b_inv_inertia * zp_cross2(r2, impulse);


 /*
  tangent / friction impulse
  
   ^
  /|\
   |
   | --- normal direction
   | 
   |       ______ tangent direction
   |       |
   •------------- >
  
  • resist sliding behaviour 
  • always perpendicular or 90°
 */
  
  va = zp_add2(a_velocity, zp_cross_sv2(a_omega, r1));
  vb = zp_add2(b_velocity, zp_cross_sv2(b_omega, r2));
	 
	 relative_vel = zp_sub2(vb, va);


		zp_vec2 tangent = zp_perp2(contact->_normal);
		float vt = zp_dot2(relative_vel, tangent);
 	j = contact->_mass_tangent * -vt;

		float friction_range = friction * contact->_accumulated_normal;

		float prev_tangent = contact->_accumulated_tangent;
		contact->_accumulated_tangent = zp_min(zp_max(prev_tangent + j, -friction_range), friction_range);
		j = contact->_accumulated_tangent - prev_tangent;

  impulse = zp_mul2(zp_stv2(j), tangent);
 	a_velocity = zp_sub2(a_velocity, zp_mul2(zp_stv2(a_inv_mass), impulse));
	 a_omega -= a_inv_inertia * zp_cross2(r1, impulse);
		 
	 b_velocity = zp_add2(b_velocity, zp_mul2(zp_stv2(b_inv_mass), impulse));
	 b_omega += b_inv_inertia * zp_cross2(r2, impulse);
 }

 body_a->_head._velocity = a_velocity;
 body_a->_head._omega = a_omega;
 
 zp_compiler_memory_barrier();

 body_b->_head._velocity = b_velocity;
 body_b->_head._omega = b_omega;
}





void zp_manifold2d_soft_relaxation(zp_manifold2d *const zp_restrict m, void *const zp_restrict w, const zp_solver_input2d *const input) {
 zp_world2d *const world = (zp_world2d*)w;
 
 if(zp_unlikely(world->_body_container._to_index_lut[m->_body_a] >= world->_body_container._size))
  return;
 zp_body2d *const body_a = (zp_body2d*)zp_container_get(&world->_body_container, m->_body_a);
  
 if(zp_unlikely(world->_body_container._to_index_lut[m->_body_b] >= world->_body_container._size))
  return;
 zp_body2d *const body_b = (zp_body2d*)zp_container_get(&world->_body_container, m->_body_b);

 zp_compiler_memory_barrier();

 zp_vec2 a_pos = body_a->_head._position;
 zp_complex a_rot = body_a->_head._rotation;
 zp_vec2 a_velocity = body_a->_head._velocity;
 float a_omega = body_a->_head._omega;
 float a_inv_mass = body_a->_head._inv_mass;
 float a_inv_inertia = body_a->_head._inv_inertia;
 float a_friction = body_a->_head._friction;

 zp_compiler_memory_barrier();

 zp_vec2 b_pos = body_b->_head._position;
 zp_complex b_rot = body_b->_head._rotation;
 zp_vec2 b_velocity = body_b->_head._velocity;
 float b_omega = body_b->_head._omega;
 float b_inv_mass = body_b->_head._inv_mass;
 float b_inv_inertia = body_b->_head._inv_inertia;
 float b_friction = body_b->_head._friction;


 float friction = zp_sqrt(a_friction * b_friction); 

 zp_vec2 r1, r2;
 zp_vec2 va, vb;
 zp_vec2 relative_vel;
 zp_vec2 impulse;
 float j, depth;


 for(uint8_t i = 0; i < m->_contact_count; i++) {
 	zp_contact2d *const contact = m->_contacts + i;

  r1 = contact->_r1;
  r2 = contact->_r2;


  r1 = zp_cmul(a_rot, r1);
  r2 = zp_cmul(b_rot, r2);

  float adjusted_depth = contact->_depth + zp_dot2(zp_sub2(r2, r1), contact->_normal);

  zp_vec2 p1 = zp_add2(a_pos, contact->_r1);
  zp_vec2 p2 = zp_add2(b_pos, contact->_r2);
  depth = -(zp_dot2(zp_sub2(p2, p1), contact->_normal) - adjusted_depth);
  

	 va = zp_add2(a_velocity, zp_cross_sv2(a_omega, r1));
  vb = zp_add2(b_velocity, zp_cross_sv2(b_omega, r2));
	 relative_vel = zp_sub2(vb, va);

		float vn = zp_dot2(relative_vel, contact->_normal);




  float slop = 0.004f;
  float penetration_error = 0.0f;

  if(depth < 0.0f) {
   float error = zp_min((contact->_depth + slop), 0.0f);
   penetration_error = -error;

  } else {
   penetration_error = contact->_depth * input->_inv_dt;
  }

  /* logarithmic with damping-like behaviour */
  penetration_error = zp_log2(penetration_error + 1.0f) * 0.69314718f;
  
  
  if(vn < contact->_bias) {
   penetration_error = -(vn - contact->_bias) + penetration_error;
  } else {
   penetration_error = -vn + penetration_error;
  }
  
  
	 j = contact->_mass_normal * penetration_error;

  
		float accumulated_normal = contact->_accumulated_normal;
		contact->_accumulated_normal = zp_max(accumulated_normal + j, 0.0f);
		j = contact->_accumulated_normal - accumulated_normal;
  
	 impulse = zp_mul2(zp_stv2(j), contact->_normal);
 	a_velocity = zp_sub2(a_velocity, zp_mul2(zp_stv2(a_inv_mass), impulse));
	 a_omega -= a_inv_inertia * zp_cross2(r1, impulse);
		 
	 b_velocity = zp_add2(b_velocity, zp_mul2(zp_stv2(b_inv_mass), impulse));
	 b_omega += b_inv_inertia * zp_cross2(r2, impulse);

  
  va = zp_add2(a_velocity, zp_cross_sv2(a_omega, r1));
  vb = zp_add2(b_velocity, zp_cross_sv2(b_omega, r2));
	 relative_vel = zp_sub2(vb, va);


		zp_vec2 tangent = zp_perp2(contact->_normal);
		float vt = zp_dot2(relative_vel, tangent);
 	j = contact->_mass_tangent * -vt;

		float friction_range = friction * contact->_accumulated_normal;

		float prev_tangent = contact->_accumulated_tangent;
		contact->_accumulated_tangent = zp_min(zp_max(prev_tangent + j, -friction_range), friction_range);
		j = contact->_accumulated_tangent - prev_tangent;

  impulse = zp_mul2(zp_stv2(j), tangent);
 	a_velocity = zp_sub2(a_velocity, zp_mul2(zp_stv2(a_inv_mass), impulse));
	 a_omega -= a_inv_inertia * zp_cross2(r1, impulse);
		 
	 b_velocity = zp_add2(b_velocity, zp_mul2(zp_stv2(b_inv_mass), impulse));
	 b_omega += b_inv_inertia * zp_cross2(r2, impulse);
 }

 body_a->_head._velocity = a_velocity;
 body_a->_head._omega = a_omega;
 
 zp_compiler_memory_barrier();

 body_b->_head._velocity = b_velocity;
 body_b->_head._omega = b_omega;
}



