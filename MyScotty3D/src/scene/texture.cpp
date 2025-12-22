
#include "texture.h"

#include <iostream>

namespace Textures {


Spectrum sample_nearest(HDR_Image const &image, Vec2 uv) {
	//clamp texture coordinates, convert to [0,w]x[0,h] pixel space:
	float x = image.w * std::clamp(uv.x, 0.0f, 1.0f);
	float y = image.h * std::clamp(uv.y, 0.0f, 1.0f);

	//the pixel with the nearest center is the pixel that contains (x,y):
	int32_t ix = int32_t(std::floor(x));
	int32_t iy = int32_t(std::floor(y));

	//texture coordinates of (1,1) map to (w,h), and need to be reduced:
	ix = std::min(ix, int32_t(image.w) - 1);
	iy = std::min(iy, int32_t(image.h) - 1);

	return image.at(ix, iy);
}

Spectrum sample_bilinear(HDR_Image const &image, Vec2 uv) {
	// A1T6: sample_bilinear
	//TODO: implement bilinear sampling strategy on texture 'image'
	float x = image.w * std::clamp(uv.x, 0.0f, 1.0f);
	float y = image.h * std::clamp(uv.y, 0.0f, 1.0f);

	//!!! use -0.5 here for 1.5 is the center
	x-=0.5;
	y-=0.5;
	// not always floor -> depends
	//the pixel with the nearest center is the pixel that contains (x,y):
	int32_t ix = int32_t(std::floor(x));

	int32_t iy = int32_t(std::floor(y));

	//texture coordinates of (1,1) map to (w,h), and need to be reduced:
	ix = std::min(ix, int32_t(image.w) - 1);
	ix = std::max(ix,0);


	iy = std::min(iy, int32_t(image.h) - 1);
	iy = std::max(iy,0);

	float delta_x,delta_y; 
	delta_x= x - ix*1.0f;
	delta_y = y - iy*1.0f;

	//corner cases: x,y=0, !
	if(delta_x<0) delta_x=0;
	if(delta_y<0) delta_y=0;

	// corner cases: x,y out of range!
	if(ix==int32_t(image.w) - 1) delta_x=0;
	if(iy==int32_t(image.h) - 1) delta_y=0;


	Spectrum t_00 = image.at(ix, iy);
	Spectrum t_10,t_01,t_11;
	// corner cases
	if(delta_x!=0 && (ix+1)<int32_t(image.w))
		t_10 = image.at(ix+1, iy);
	else
		t_10 = image.at(ix, iy);
	
	if(delta_y!=0 && (iy+1)<int32_t(image.h))
		t_01 = image.at(ix, iy+1);
	else
		t_01 = image.at(ix, iy);

	if(delta_y!=0 && delta_x!=0 && ix+1<int32_t(image.w)  && iy+1<int32_t(image.h))
		t_11 = image.at(ix+1, iy+1);
	else if(delta_y==0 && delta_x!=0 && ix+1<int32_t(image.w))
		t_11 = image.at(ix+1, iy);
	else if(delta_y!=0 && delta_x==0 && iy+1<int32_t(image.h))
		t_11 = image.at(ix, iy+1);
	else
		t_11 = image.at(ix,iy);
		// return t_00;
	
	Spectrum res,tx,ty;
	tx = (1.0f-delta_x) * t_00 + delta_x * t_10;
	ty = (1.0f-delta_x) * t_01 + delta_x * t_11;
	res =  (1.0f-delta_y)*tx + delta_y*ty;

	return res;
	// return sample_nearest(image, uv); //placeholder so image doesn't look blank
}


Spectrum sample_trilinear(HDR_Image const &base, std::vector< HDR_Image > const &levels, Vec2 uv, float lod) {
	// A1T6: sample_trilinear
	int ilod = (int)std::floor(lod);
	float delta_d = lod-1.0f*ilod;
	Spectrum td,td_next;
	Vec2 px_to_uv;
	// std::cout<<ilod;
	if(ilod==0){
		px_to_uv = Vec2(1.0f / base.w, 1.0f / base.h);
		td = sample_bilinear(base, px_to_uv*uv);
	}
	else{
		px_to_uv = Vec2(1.0f / levels[ilod-1].w*1.0f, 1.0f / levels[ilod-1].h*1.0f);
		// std::cout<<px_to_uv.y;
		td = sample_bilinear(levels[ilod-1], px_to_uv*uv);
	}
	if(delta_d!=0){
		px_to_uv = Vec2(1.0f / levels[ilod].w, 1.0f / levels[ilod].h);
		td_next = sample_bilinear(levels[ilod], px_to_uv*uv);
	}else{
		// std::cout<<(delta_d!=0)<<"herhe";
		return td;
	}

	return (1-delta_d) * td + delta_d*td_next;

	// return sample_nearest(base, uv); //placeholder so image doesn't look blank
}

/*
 * generate_mipmap- generate mipmap levels from a base image.
 *  base: the base image
 *  levels: pointer to vector of levels to fill (must not be null)
 *
 * generates a stack of levels [1,n] of sizes w_i, h_i, where:
 *   w_i = max(1, floor(w_{i-1})/2)
 *   h_i = max(1, floor(h_{i-1})/2)
 *  with:
 *   w_0 = base.w
 *   h_0 = base.h
 *  and n is the smalles n such that w_n = h_n = 1
 *
 * each level should be calculated by downsampling a blurred version
 * of the previous level to remove high-frequency detail.
 *
 */
void generate_mipmap(HDR_Image const &base, std::vector< HDR_Image > *levels_) {
	assert(levels_);
	auto &levels = *levels_;


	{ // allocate sublevels sufficient to scale base image all the way to 1x1:
		int32_t num_levels = static_cast<int32_t>(std::log2(std::max(base.w, base.h)));
		assert(num_levels >= 0);

		levels.clear();
		levels.reserve(num_levels);

		uint32_t width = base.w;
		uint32_t height = base.h;
		for (int32_t i = 0; i < num_levels; ++i) {
			assert(!(width == 1 && height == 1)); //would have stopped before this if num_levels was computed correctly

			width = std::max(1u, width / 2u);
			height = std::max(1u, height / 2u);

			levels.emplace_back(width, height);
		}
		assert(width == 1 && height == 1);
		assert(levels.size() == uint32_t(num_levels));
	}

	//now fill in the levels using a helper:
	//downsample:
	// fill in dst to represent the low-frequency component of src
	auto downsample = [](HDR_Image const &src, HDR_Image &dst) {
		//dst is half the size of src in each dimension:
		assert(std::max(1u, src.w / 2u) == dst.w);
		assert(std::max(1u, src.h / 2u) == dst.h);

		// A1T6: generate
		//TODO: Write code to fill the levels of the mipmap hierarchy by downsampling
		for(uint32_t i=0;i<dst.w;i++){
			for(uint32_t j=0;j<dst.h;j++){
					// dst.at(i,j) = (src.at(i*2,j*2) + src.at(i*2+1,j*2+1) +
					// src.at(i*2+1,j*2) + src.at(i*2,j*2+1))/4.0;
					Vec2 px_to_uv(1.0f / src.w, 1.0f / src.h);
					Spectrum got = Textures::sample_bilinear(src, px_to_uv 
					* Vec2(i*2.0f+1.0f, j*2.0f+1.0f));
					dst.at(i,j) = got;
			}
		}
		//Be aware that the alignment of the samples in dst and src will be different depending on whether the image is even or odd.

	};

	std::cout << "Regenerating mipmap (" << levels.size() << " levels): [" << base.w << "x" << base.h << "]";
	std::cout.flush();
	for (uint32_t i = 0; i < levels.size(); ++i) {
		HDR_Image const &src = (i == 0 ? base : levels[i-1]);
		HDR_Image &dst = levels[i];
		std::cout << " -> [" << dst.w << "x" << dst.h << "]"; std::cout.flush();

		downsample(src, dst);
	}
	std::cout << std::endl;
	
}

Image::Image(Sampler sampler_, HDR_Image const &image_) {
	sampler = sampler_;
	image = image_.copy();
	update_mipmap();
}

Spectrum Image::evaluate(Vec2 uv, float lod) const {
	if (sampler == Sampler::nearest) {
		return sample_nearest(image, uv);
	} else if (sampler == Sampler::bilinear) {
		return sample_bilinear(image, uv);
	} else {
		return sample_trilinear(image, levels, uv, lod);
	}
}

void Image::update_mipmap() {
	if (sampler == Sampler::trilinear) {
		generate_mipmap(image, &levels);
	} else {
		levels.clear();
	}
}

GL::Tex2D Image::to_gl() const {
	return image.to_gl(1.0f);
}

void Image::make_valid() {
	update_mipmap();
}

Spectrum Constant::evaluate(Vec2 uv, float lod) const {
	return color * scale;
}

} // namespace Textures
bool operator!=(const Textures::Constant& a, const Textures::Constant& b) {
	return a.color != b.color || a.scale != b.scale;
}

bool operator!=(const Textures::Image& a, const Textures::Image& b) {
	return a.image != b.image;
}

bool operator!=(const Texture& a, const Texture& b) {
	if (a.texture.index() != b.texture.index()) return false;
	return std::visit(
		[&](const auto& data) { return data != std::get<std::decay_t<decltype(data)>>(b.texture); },
		a.texture);
}
