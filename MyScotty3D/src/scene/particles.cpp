
#include "particles.h"

bool Particles::Particle::update(const PT::Aggregate &scene, Vec3 const &gravity, const float radius, const float dt) {

	//A4T4: particle update
	// Compute the trajectory of this particle for the next dt seconds.
	// dt is small
	
	// (1) Build a ray representing the particle's path as if it travelled at constant velocity.
	// Ray ray_path = Ray();
	// ray_path.point = position;
	// ray_path.dir = velocity;

	// float v_norm = velocity.norm();
	// float max_distance = dt*v_norm;//radius;
	// ray_path.dist_bounds = Vec2{0.0f, max_distance};
	// // max distance!
	// // (2) Intersect the ray with the scene and account for collisions. Be careful when placing
	// // collision points using the particle radius. Move the particle to its next position.
	// auto collision = scene.hit(ray_path);
	// // float t_cost;
	// // check if hit
	// if( !collision.hit){
	// 	// 	// printf("1:%f,%f,%f\n",position.x,position.y,position.z);
	// 	// 	// printf("2:%f,%f\n",t_cost,velocity.y);
	// 	// first do position, then do velocity!!!!
	// 	position+=velocity*dt;
	// 	// (3) Account for acceleration due to gravity after updating position.
	// 	velocity += gravity *dt;
	// 	age -= dt;

	// 	printf("not hit %f,%f,%f\n",position.x,position.y,position.z);
	// 	if(age<=0)
	// 		return false;
	// 	return true;
	// }else{
	// 	auto col_norm = collision.normal;
	// 	auto d_tmp = collision.position - (dt*velocity+position);
	// 	auto check = dot(col_norm,-d_tmp);
	// 	if(check>radius){
	// 		position+=velocity*dt;
	// 		// (3) Account for acceleration due to gravity after updating position.
	// 		velocity += gravity *dt;
	// 		age -= dt;

	// 		printf("not hit %f,%f,%f\n",position.x,position.y,position.z);
	// 		if(age<=0)
	// 			return false;
	// 		return true;
	// 	}
	// }

	float t_left = dt;
	while(t_left>0){
		Ray ray_path = Ray();
		// ray_path = Ray();
		// printf("old v %f\n",velocity.y);
		ray_path.point = position;
		ray_path.dir = velocity.unit();
		// max distance!
		auto v_norm = velocity.norm();
		auto max_distance = t_left*v_norm;//radius;
		// ray_path.dist_bounds = Vec2{0.0f, max_distance};
		// (2) Intersect the ray with the scene and account for collisions. Be careful when placing
		// collision points using the particle radius. Move the particle to its next position.
		auto collision = scene.hit(ray_path);
		// float t_cost;
		// check if hit
		float t_consume;
		// printf("hit or not %d\n",collision.hit);
		// printf("old v2 %f\n",velocity.y);
		if(!collision.hit){
			position+=velocity*t_left;
			// (3) Account for acceleration due to gravity after updating position.
			velocity += gravity *t_left;
			// t_consume = t_left;
			break;
		}else{
			auto col_norm2 = collision.normal;

			auto cos_t = std::abs(dot(col_norm2,velocity))/v_norm;
			auto new_d = radius/cos_t;

			// auto d_tmp2 = collision.position - (t_left*velocity+position);

			// auto check2 = dot(col_norm2,-d_tmp2);
			// printf("check2 %f\n",check2);
			// if(collision.distance<new_d){
			// 	velocity = velocity - 2 * dot(velocity, col_norm2) * col_norm2;
			// 	t_consume = 0.f;
			// }
			if((collision.distance-new_d)>max_distance){
				position+=velocity*t_left;
				// (3) Account for acceleration due to gravity after updating position.
				velocity += gravity *t_left;
				// t_consume = t_left;
				break;
			}
			else{
				// auto cos_t = std::abs(dot(col_norm2,velocity))/(velocity.norm());
				if(v_norm==0){
					position+=velocity*t_left;
					// (3) Account for acceleration due to gravity after updating position.
					velocity += gravity *t_left;
					// t_consume = t_left;
					break;
				}
				t_consume = (collision.distance-new_d)/v_norm;
				// printf("t %f\n",t_consume);
				// position = collision.position -new_d * velocity.normalize();
				if(t_consume>0){
					position += velocity * t_consume;
				}else{
					// printf("here v %f\n",velocity.y);
					t_consume = 0.f;
				}
				// printf("new v %f\n",col_norm2);
				velocity = velocity - 2 * dot(velocity, col_norm2) * col_norm2;
				// printf("v2 %f,%f,%f\n",position.x,position.y,position.z);
				// printf("new v %f\n",velocity.y);
				// break;
			}
		}
		velocity += gravity * t_consume;
		t_left -=t_consume;
	}
	// time cost?
	// velocity += gravity * t_cost;
	// 		t_cost = collision.distance/v_norm; 
	// 	// (3) Account for acceleration due to gravity after updating position.
	// 	// printf("1:%f,%f,%f\n",position.x,position.y,position.z);
	// 	// (4) Repeat until the entire time step has been consumed.
	// 	t_left-=t_cost;
	// }
	// (5) Decrease the particle's age and return 'false' if it should be removed.
	// printf("new v %f\n",velocity.y);
	// printf("hit %f,%f,%f\n",position.x,position.y,position.z);
	age -= dt;
	if(age<=0)
		return false;
	return true;
}

void Particles::advance(const PT::Aggregate& scene, const Mat4& to_world, float dt) {

	if(step_size < EPS_F) return;

	step_accum += dt;

	while(step_accum > step_size) {
		step(scene, to_world);
		step_accum -= step_size;
	}
}

void Particles::step(const PT::Aggregate& scene, const Mat4& to_world) {

	std::vector<Particle> next;
	next.reserve(particles.size());

	for(Particle& p : particles) {
		if(p.update(scene, gravity, radius, step_size)) {
			next.emplace_back(p);
		}
	}

	if(rate > 0.0f) {

		//helpful when emitting particles:
		float cos = std::cos(Radians(spread_angle) / 2.0f);

		//will emit particle i when i == time * rate
		//(i.e., will emit particle when time * rate hits an integer value.)
		//so need to figure out all integers in [current_step, current_step+1) * step_size * rate
		//compute the range:
		double begin_t = current_step * double(step_size) * double(rate);
		double end_t = (current_step + 1) * double(step_size) * double(rate);

		uint64_t begin_i = uint64_t(std::max(0.0, std::ceil(begin_t)));
		uint64_t end_i = uint64_t(std::max(0.0, std::ceil(end_t)));

		//iterate all integers in [begin, end):
		for (uint64_t i = begin_i; i < end_i; ++i) {
			//spawn particle 'i':

			float y = lerp(cos, 1.0f, rng.unit());
			float t = 2 * PI_F * rng.unit();
			float d = std::sqrt(1.0f - y * y);
			Vec3 dir = initial_velocity * Vec3(d * std::cos(t), y, d * std::sin(t));

			Particle p;
			p.position = to_world * Vec3(0.0f, 0.0f, 0.0f);
			p.velocity = to_world.rotate(dir);
			p.age = lifetime; //NOTE: could adjust lifetime based on index
			next.push_back(p);
		}
	}

	particles = std::move(next);
	current_step += 1;
}

void Particles::reset() {
	particles.clear();
	step_accum = 0.0f;
	current_step = 0;
	rng.seed(seed);
}

bool operator!=(const Particles& a, const Particles& b) {
	return a.gravity != b.gravity
	|| a.radius != b.radius
	|| a.initial_velocity != b.initial_velocity
	|| a.spread_angle != b.spread_angle
	|| a.lifetime != b.lifetime
	|| a.rate != b.rate
	|| a.step_size != b.step_size
	|| a.seed != b.seed;
}
