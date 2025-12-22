
#include "material.h"
#include "../util/rand.h"

namespace Materials {

Vec3 reflect(Vec3 dir) {
	//A3T5 Materials - reflect helper

    // Return direction to incoming light that would be
	// reflected out in direction dir from surface
	// with normal (0,1,0)

    return Vec3{-dir.x,dir.y,-dir.z};
}

Vec3 refract(Vec3 out_dir, float index_of_refraction, bool& was_internal) {
	//A3T5 Materials - refract helper

	// Use Snell's Law to refract out_dir through the surface.
	// Return the refracted direction. Set was_internal to true if
	// refraction does not occur due to total internal reflection,
	// and false otherwise.

	// The surface normal is (0,1,0)
	// out_dir=out_dir.normalize();
	float cos_i = out_dir.y;
	// !!!! >0 then divide
	if(cos_i>0){
		index_of_refraction = 1.f/index_of_refraction;
	}

	float sin_i_2 = 1.f-cos_i*cos_i; 
	float ratio_z = abs(out_dir.z)/(abs(out_dir.x)+abs(out_dir.z));
	float ratio_x = abs(out_dir.x)/(abs(out_dir.x)+abs(out_dir.z));
	// shouldn't happen
	// if(sin_i_2<0){
	// 	sin_i_2 = 0.0f;
	// }
	float index_2 = index_of_refraction * index_of_refraction;
	float cos_t_2 = 1 - index_2*sin_i_2;
	// should not be larger than 90
	if(cos_t_2<=0){
		was_internal = true;
		return reflect(out_dir);
	}

	float ref_x;
	float ref_z;
	was_internal = false;

	float cos_t = sqrt(cos_t_2);
	if(out_dir.y>0){
		cos_t = -cos_t;
	}

	float sin_t_2 = 1-cos_t_2;
	float sin_t = sqrt(sin_t_2);
	ref_z = sin_t * ratio_z;
	ref_x = ratio_x*sin_t;

	if(out_dir.x>0){
		ref_x = -ref_x;
	}
	if(out_dir.z>0){
		ref_z = -ref_z;
	}

	Vec3 ret = Vec3{ref_x,cos_t,ref_z}.normalize();

	return ret;
}

float schlick(Vec3 in_dir, float index_of_refraction) {
	//A3T5 Materials - Schlick's approximation helper

	// Implement Schlick's approximation of the Fresnel reflection factor.
	// Schlick’s approximation
	// r0 = n2/n1
	float cos_i = in_dir.y;
	// !!!! >0 then divide
	if(cos_i>0){
		index_of_refraction = 1.f/index_of_refraction;
	}

	float n1 = 1.f;
	float n2 = index_of_refraction;
	float r0 = std::pow(((n1-n2)/(n1+n2)),2.f);
	float r = r0 + (1-r0) * powf((1-abs(in_dir.y)),5.f);

	// bool is_internal;
    // float cos_t = std::abs(refract(in_dir,index_of_refraction,is_internal).y);
	// if(is_internal){
	// 	return 1.;
	// }
	// float cos_i = in_dir.y;
	// float ni = 1.f;
	// float nt = index_of_refraction;
	// float r1 = (nt*cos_i - ni*cos_t)/(nt*cos_i + ni*cos_t);
	// float r2 = (ni*cos_i - nt*cos_t)/(ni*cos_i + nt*cos_t);

	// rf = (r1*r1 + r2*r2)/2.f;



	// printf("%f\n",r);

	return r;
}

Spectrum Lambertian::evaluate(Vec3 out, Vec3 in, Vec2 uv) const {
	//A3T4: Materials - Lambertian BSDF evaluation

    // Compute the ratio of outgoing/incoming radiance when light from in_dir
    // is reflected through out_dir: (albedo / PI_F) * cos(theta).
    // Note that for Scotty3D, y is the 'up' direction.
	float cos_theta = in.y;
	// Call albedo.lock()->evaluate(uv) to get the albedo
	Spectrum albedo_s = albedo.lock()->evaluate(uv);
	Spectrum reflection;
	if(cos_theta>0)
		reflection = (albedo_s / PI_F) * cos_theta;
	else 
		reflection = Spectrum{};
    return reflection;
}

Scatter Lambertian::scatter(RNG &rng, Vec3 out, Vec2 uv) const {
	//A3T4: Materials - Lambertian BSDF scattering
	//Select a scattered light direction at random from the Lambertian BSDF

	[[maybe_unused]] Samplers::Hemisphere::Cosine sampler; //this will be useful

	// with direction and attenuation components
	Scatter ret;
	//TODO: sample the direction the light was scatter from from a cosine-weighted hemisphere distribution:
	ret.direction = sampler.sample(rng);

	// attenuation component via Lambertian::evaluate
	//TODO: compute the attenuation of the light using Lambertian::evaluate():
	ret.attenuation = evaluate(out, ret.direction, uv);

	return ret;
}

float Lambertian::pdf(Vec3 out, Vec3 in) const {
	//A3T4: Materials - Lambertian BSDF probability density function
    // Compute the PDF for sampling in_dir from the cosine-weighted hemisphere distribution.
	[[maybe_unused]] Samplers::Hemisphere::Cosine sampler; //this might be handy!
	float res = in.y/PI_F;
	if(res<=0){
		res = 0;
	}
    return res;
}

Spectrum Lambertian::emission(Vec2 uv) const {
	return {};
}

std::weak_ptr<Texture> Lambertian::display() const {
	return albedo;
}

void Lambertian::for_each(const std::function<void(std::weak_ptr<Texture>&)>& f) {
	f(albedo);
}

Spectrum Mirror::evaluate(Vec3 out, Vec3 in, Vec2 uv) const {
	return {};
}

Scatter Mirror::scatter(RNG &rng, Vec3 out, Vec2 uv) const {
	//A3T5: mirror

	// Use reflect to compute the new direction
	// Don't forget that this is a discrete material!
	// Similar to albedo, reflectance represents the ratio of incoming light to reflected light

    Scatter ret;
    ret.direction = reflect(out);
    ret.attenuation = reflectance.lock()->evaluate(uv);
    return ret;
}

float Mirror::pdf(Vec3 out, Vec3 in) const {
	return 0.0f;
}

Spectrum Mirror::emission(Vec2 uv) const {
	return {};
}

std::weak_ptr<Texture> Mirror::display() const {
	return reflectance;
}

void Mirror::for_each(const std::function<void(std::weak_ptr<Texture>&)>& f) {
	f(reflectance);
}

Spectrum Refract::evaluate(Vec3 out, Vec3 in, Vec2 uv) const {
	return {};
}

Scatter Refract::scatter(RNG &rng, Vec3 out, Vec2 uv) const {
	//A3T5 - refract

	// Use refract to determine the new direction - what happens in the total internal reflection case?
    // Be wary of your eta1/eta2 ratio - are you entering or leaving the surface?
	// Don't forget that this is a discrete material!
	// For attenuation, be sure to take a look at the Specular Transimission section of the PBRT textbook for a derivation
	//  You do not need to scale by the Fresnel Coefficient - you'll only need to account for the correct ratio of indices of refraction

    Scatter ret;
	float index = ior; // 
	// index = ior;
	bool is_internal;
    ret.direction = refract(out,index,is_internal);
	if(!is_internal){
    	ret.attenuation = transmittance.lock()->evaluate(uv);
		if(out.y>0){
			ret.attenuation *= ((1.f/index) * (1.f/index));
		}
		else{
			ret.attenuation *= (index*index);
		}
	}
	else {
		ret.attenuation = Spectrum{1.f};
	}

	// float cos_t = ret.direction.y;
	// float cos_i = out.y;
	// float ni = 1.f;
	// float nt = index;
	// float r1 = (nt*cos_i - ni*cos_t)/(nt*cos_i + ni*cos_t);
	// float r2 = (ni*cos_i - nt*cos_t)/(ni*cos_i + nt*cos_t);
	// float rf = (r1*r1 + r2*r2)/2.f;

	// ret.attenuation *= rf;

    return ret;
}

float Refract::pdf(Vec3 out, Vec3 in) const {
	return 0.0f;
}

Spectrum Refract::emission(Vec2 uv) const {
	return {};
}

bool Refract::is_emissive() const {
	return false;
}

bool Refract::is_specular() const {
	return true;
}

bool Refract::is_sided() const {
	return true;
}

std::weak_ptr<Texture> Refract::display() const {
	return transmittance;
}

void Refract::for_each(const std::function<void(std::weak_ptr<Texture>&)>& f) {
	f(transmittance);
}

Spectrum Glass::evaluate(Vec3 out, Vec3 in, Vec2 uv) const {
	return {};
}

Scatter Glass::scatter(RNG &rng, Vec3 out, Vec2 uv) const {
	//A3T5 - glass

    // (1) Compute Fresnel coefficient. Tip: Schlick's approximation.
	float index_of_refraction = ior; // ior?
	float fresnel_coef = schlick(out,index_of_refraction);
    // (2) Reflect or refract probabilistically based on Fresnel coefficient. Tip: RNG::coin_flip
	bool choose = rng.coin_flip(fresnel_coef);
    // (3) Compute attenuation based on reflectance or transmittance

    // Be wary of your eta1/eta2 ratio - are you entering or leaving the surface?
    // What happens upon total internal reflection?
    // When debugging Glass, it may be useful to compare to a pure-refraction BSDF
	// For attenuation, be sure to take a look at the Specular Transimission section of the PBRT textbook for a derivation
	//  You do not need to scale by the Fresnel Coefficient - you'll only need to account for the correct ratio of indices of refraction

    Scatter ret;
    ret.direction = Vec3();

	if(choose){
		ret.direction = reflect(out);
		ret.attenuation = reflectance.lock()->evaluate(uv);
	}else{
		bool is_internal;
    	ret.direction = refract(out,index_of_refraction,is_internal);
		if(!is_internal){
			ret.attenuation = transmittance.lock()->evaluate(uv);
			if(out.y>0){
				ret.attenuation *= ((1.f/index_of_refraction) * (1.f/index_of_refraction));
			}
			else{
				ret.attenuation *= (index_of_refraction*index_of_refraction);
			}
		}else{
				ret.attenuation = Spectrum{1.f};
		}
	}
    
    return ret;
}

float Glass::pdf(Vec3 out, Vec3 in) const {
	return 0.0f;
}

Spectrum Glass::emission(Vec2 uv) const {
	return {};
}

bool Glass::is_emissive() const {
	return false;
}

bool Glass::is_specular() const {
	return true;
}

bool Glass::is_sided() const {
	return true;
}

std::weak_ptr<Texture> Glass::display() const {
	return transmittance;
}

void Glass::for_each(const std::function<void(std::weak_ptr<Texture>&)>& f) {
	f(reflectance);
	f(transmittance);
}

Spectrum Emissive::evaluate(Vec3 out, Vec3 in, Vec2 uv) const {
	return {};
}

Scatter Emissive::scatter(RNG &rng, Vec3 out, Vec2 uv) const {
	Scatter ret;
	ret.direction = {};
	ret.attenuation = {};
	return ret;
}

float Emissive::pdf(Vec3 out, Vec3 in) const {
	return 0.0f;
}

Spectrum Emissive::emission(Vec2 uv) const {
	return emissive.lock()->evaluate(uv);
}

bool Emissive::is_emissive() const {
	return true;
}

bool Emissive::is_specular() const {
	return true;
}

bool Emissive::is_sided() const {
	return false;
}

std::weak_ptr<Texture> Emissive::display() const {
	return emissive;
}

void Emissive::for_each(const std::function<void(std::weak_ptr<Texture>&)>& f) {
	f(emissive);
}

} // namespace Materials

bool operator!=(const Materials::Lambertian& a, const Materials::Lambertian& b) {
	return a.albedo.lock() != b.albedo.lock();
}

bool operator!=(const Materials::Mirror& a, const Materials::Mirror& b) {
	return a.reflectance.lock() != b.reflectance.lock();
}

bool operator!=(const Materials::Refract& a, const Materials::Refract& b) {
	return a.transmittance.lock() != b.transmittance.lock() || a.ior != b.ior;
}

bool operator!=(const Materials::Glass& a, const Materials::Glass& b) {
	return a.reflectance.lock() != b.reflectance.lock() ||
	       a.transmittance.lock() != b.transmittance.lock() || a.ior != b.ior;
}

bool operator!=(const Materials::Emissive& a, const Materials::Emissive& b) {
	return a.emissive.lock() != b.emissive.lock();
}

bool operator!=(const Material& a, const Material& b) {
	if (a.material.index() != b.material.index()) return false;
	return std::visit(
		[&](const auto& material) {
			return material != std::get<std::decay_t<decltype(material)>>(b.material);
		},
		a.material);
}
