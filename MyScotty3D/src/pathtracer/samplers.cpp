
#include "samplers.h"
#include "../util/rand.h"

constexpr bool IMPORTANCE_SAMPLING =  true;

namespace Samplers {

Vec2 Rect::sample(RNG &rng) const {
	//A3T1 - step 2 - supersampling

    // Return a point selected uniformly at random from the rectangle [0,size.x)x[0,size.y)
    // Useful function: rng.unit()

    return Vec2{rng.unit()*size.x, rng.unit()*size.y};
}

Vec2 Jittered::sample(RNG &rng,int num_sum)const{
		int num_per_row = int(sqrt(num_sum));
		float row_size = size.x*1.f/num_per_row;
		float col_size = size.y*1.f/num_per_row;
		int i = index/num_per_row;
		int j = index%num_per_row;
		float sample_x = (i+rng.unit()) * row_size;
		float sample_y = (j+rng.unit()) * col_size;
		return Vec2{sample_x,sample_y};

}

Vec2 MultiJittered::sample(RNG &rng,int num_sum)const{
		int num_per_row = int(sqrt(num_sum));
		float row_size = size.x*1.f/num_per_row;
		float col_size = size.y*1.f/num_per_row;
		int i = index/num_per_row;
		int j = index%num_per_row;
		float sample_x = i * row_size + (j + rng.unit()) * row_size / num_per_row;
		float sample_y = j * col_size + (i + rng.unit()) * col_size / num_per_row;
		return Vec2{sample_x,sample_y};

}



float Rect::pdf(Vec2 at) const {
	if (at.x < 0.0f || at.x > size.x || at.y < 0.0f || at.y > size.y) return 0.0f;
	return 1.0f / (size.x * size.y);
}

Vec3 Point::sample(RNG &rng) const {
	return point;
}

float Point::pdf(Vec3 at) const {
	return at == point ? 1.0f : 0.0f;
}

Vec3 Triangle::sample(RNG &rng) const {
	float u = std::sqrt(rng.unit());
	float v = rng.unit();
	float a = u * (1.0f - v);
	float b = u * v;
	return a * v0 + b * v1 + (1.0f - a - b) * v2;
}

float Triangle::pdf(Vec3 at) const {
	float a = 0.5f * cross(v1 - v0, v2 - v0).norm();
	float u = 0.5f * cross(at - v1, at - v2).norm() / a;
	float v = 0.5f * cross(at - v2, at - v0).norm() / a;
	float w = 1.0f - u - v;
	if (u < 0.0f || v < 0.0f || w < 0.0f) return 0.0f;
	if (u > 1.0f || v > 1.0f || w > 1.0f) return 0.0f;
	return 1.0f / a;
}

Vec3 Hemisphere::Uniform::sample(RNG &rng) const {

	float Xi1 = rng.unit();
	float Xi2 = rng.unit();

	float theta = std::acos(Xi1);

	float phi = 2.0f * PI_F * Xi2;

	float xs = std::sin(theta) * std::cos(phi);
	float ys = std::cos(theta);
	float zs = std::sin(theta) * std::sin(phi);

	return Vec3(xs, ys, zs);
}

float Hemisphere::Uniform::pdf(Vec3 dir) const {
	if (dir.y < 0.0f) return 0.0f;
	return 1.0f / (2.0f * PI_F);
}

Vec3 Hemisphere::Cosine::sample(RNG &rng) const {

	float phi = rng.unit() * 2.0f * PI_F;
	float cos_t = std::sqrt(rng.unit());

	float sin_t = std::sqrt(1 - cos_t * cos_t);
	float x = std::cos(phi) * sin_t;
	float z = std::sin(phi) * sin_t;
	float y = cos_t;

	return Vec3(x, y, z);
}

float Hemisphere::Cosine::pdf(Vec3 dir) const {
	if (dir.y < 0.0f) return 0.0f;
	return dir.y / PI_F;
}

Vec3 Sphere::Uniform::sample(RNG &rng) const {
	//A3T7 - sphere sampler

    // Generate a uniformly random point on the unit sphere.
    // Tip: start with Hemisphere::Uniform
	float theta = PI_F*rng.unit();
	float alpha = 2.f*PI_F*rng.unit();
	float y = cosf(theta);
	float z = sinf(theta) * sinf(alpha);
	float x  = sinf(theta) * cosf(alpha);
    return Vec3{x,y,z};
	// auto sample = hemi.sample(rng);
	// printf("%f,%f,%f",sample.x,sample.y,sample.z);
	// return hemi.sample(rng);
}

float Sphere::Uniform::pdf(Vec3 dir) const {
	return 1.0f / (4.0f * PI_F);
}

Sphere::Image::Image(const HDR_Image& image) {
    //A3T7 - image sampler init

    // Set up importance sampling data structures for a spherical environment map image.
    // You may make use of the _pdf, _cdf, and total members, or create your own.

    const auto [_w, _h] = image.dimension();
    w = _w;
    h = _h;

	_pdf = std::vector<float>();
	_cdf = std::vector<float>();

	float sum_luma = 0.f;
	for (uint32_t i = 0; i < w; i++) {
        for (uint32_t j = 0; j < h; j++) {
			float luma = image.at(i,j).luma();
			sum_luma+=luma;
		}
	}

	float sum_j = 0.f;
	for (uint32_t i = 0; i < h; i++) {
        for (uint32_t j = 0; j < w; j++) {
			//  is actually the bottom left of the HDR image, not the top left 
			// float theta;
			// theta = (i+0.5f) * 1.0f * (PI_F/(h *1.0f));
			// float sin_theta = sinf(theta);
			//according to Piazza answer, no jacobian when generating, but need for pdf!
			//float jacobin = (w*h)/(2.f*PI_F*PI_F*sin_theta);
			float j_one = image.at(j,i).luma()/ sum_luma;
			sum_j+=j_one;
			// printf("%f,%f,%f,%f\n",j_one,sum_j,sin_theta,image.at(j,i).luma());
			_pdf.push_back(j_one);
			_cdf.push_back(sum_j);
		}
	}

	// normalize
	for(uint32_t i=0;i<w*h;i++){
		_pdf[i]/=sum_j;
		
		_cdf[i]/=sum_j;
		// printf("%f\n",_cdf[i]);
	}
}

Vec3 Sphere::Image::sample(RNG &rng) const {
	if(!IMPORTANCE_SAMPLING) {
		// Step 1: Uniform sampling
		// Declare a uniform sampler and return its sample
		Uniform sampler;
    	return sampler.sample(rng);
	} else {
		// Step 2: Importance sampling
		// Use your importance sampling data structure to generate a sample direction.
		// Tip: std::upper_bound
    	// return Vec3{1.f};
		float choose_cdf = rng.unit();
		auto cdf_it = std::upper_bound(_cdf.begin(),_cdf.end(),choose_cdf);
		int cdf_idx = 0;
		if(cdf_it==_cdf.end()){
			cdf_idx = int(_cdf.size()-1);
		}
		cdf_idx = int(cdf_it - _cdf.begin());

		int index = cdf_idx;

		// printf("%d\n",index);
		int x = index % w;
        int y = index / w;

		// printf("%d,%d\n",x,y);

        float phi = (2.0f * PI_F * x) / w;
        float theta = PI_F * y / h;

		// printf("%f,%f\n",phi,theta);
		float s_x = cos(phi)*sin(theta)*1.f;
		// need to use - here!!! for 0-h vs pi-0
		float s_y = -cos(theta)*1.f;
		float s_z = sin(phi)*sin(theta)*1.f;
        return Vec3(s_x,s_y ,s_z );
	}
}

float Sphere::Image::pdf(Vec3 dir) const {
    if(!IMPORTANCE_SAMPLING) {
		// Step 1: Uniform sampling
		// Declare a uniform sampler and return its pdf
    	return 1.f/(4.f*PI_F);
	} else {
		// A3T7 - image sampler importance sampling pdf
		// What is the PDF of this distribution at a particular direction?
    	// return 0.f;
		// get back
		float theta = acos(-dir.y);
		// float sin_theta = sin(theta);
		int y = int((theta/PI_F) * h);

		float phi = atan2f(dir.z,dir.x);
		//!!! phi<0:
		if(phi<0){
			phi+=2*PI_F;
		}

		int x = int((phi/(2.f*PI_F)) * w);
		// printf("%f,%f\n",phi,theta);
		// printf("%d,%d\n",x,y);
		
		

		y = std::max(0,y);
		y = std::min(int(h-1),y);

		x = std::max(0,x);
		x= std::min(int(w-1),x);

		// according to Piazza answer, no jacobian when generating, but need for pdf!
		float sin_theta = sinf(theta);
		float jacobin = (1.f*w*h)/(2.f*PI_F*PI_F*sin_theta);

		int idx = y*w + x;
		// printf("%d\n",idx);
		return _pdf.at(idx) * jacobin;

	}
	
}

} // namespace Samplers
