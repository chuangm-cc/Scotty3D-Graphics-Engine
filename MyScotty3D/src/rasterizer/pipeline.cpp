// clang-format off
#include "pipeline.h"

#include <iostream>

#include "../lib/log.h"
#include "../lib/mathlib.h"
#include "framebuffer.h"
#include "sample_pattern.h"
template<PrimitiveType primitive_type, class Program, uint32_t flags>
void Pipeline<primitive_type, Program, flags>::run(std::vector<Vertex> const& vertices,
                                                   typename Program::Parameters const& parameters,
                                                   Framebuffer* framebuffer_) {
	// Framebuffer must be non-null:
	assert(framebuffer_);
	auto& framebuffer = *framebuffer_;

	// A1T7: sample loop
	// TODO: update this function to rasterize to *all* sample locations in the framebuffer.
	//  	 This will probably involve inserting a loop of the form:
	// 		 	std::vector< Vec3 > const &samples = framebuffer.sample_pattern.centers_and_weights;
	//      	for (uint32_t s = 0; s < samples.size(); ++s) { ... }
	//   	 around some subset of the code.
	// 		 You will also need to transform the input and output of the rasterize_* functions to
	// 	     account for the fact they deal with pixels centered at (0.5,0.5).
	

	std::vector<ShadedVertex> shaded_vertices;
	// shaded_vertices.reserve(vertices.size());

	shaded_vertices.reserve(vertices.size());
	//--------------------------
	// shade vertices:
	for (auto const& v : vertices) {
		ShadedVertex sv;
		Program::shade_vertex(parameters, v.attributes, &sv.clip_position, &sv.attributes);
		shaded_vertices.emplace_back(sv);
	}

	// added here
	// defination
	// rest in below
	std::vector< Vec3 > const &samples = framebuffer.sample_pattern.centers_and_weights;

	//--------------------------
	// assemble + clip + homogeneous divide vertices:
	std::vector<ClippedVertex> clipped_vertices;

	// reserve some space to avoid reallocations later:
	if constexpr (primitive_type == PrimitiveType::Lines) {
		// clipping lines can never produce more than one vertex per input vertex:
		clipped_vertices.reserve(shaded_vertices.size());
	} else if constexpr (primitive_type == PrimitiveType::Triangles) {
		// clipping triangles can produce up to 8 vertices per input vertex:
		clipped_vertices.reserve(shaded_vertices.size() * 8);
	}
	// clang-format off

	//coefficients to map from clip coordinates to framebuffer (i.e., "viewport") coordinates:
	//x: [-1,1] -> [0,width]
	//y: [-1,1] -> [0,height]
	//z: [-1,1] -> [0,1] (OpenGL-style depth range)
	Vec3 const clip_to_fb_scale = Vec3{
		framebuffer.width / 2.0f,
		framebuffer.height / 2.0f,
		0.5f
	};
	Vec3 const clip_to_fb_offset = Vec3{
		0.5f * framebuffer.width,
		0.5f * framebuffer.height,
		0.5f
	};

	// helper used to put output of clipping functions into clipped_vertices:
	auto emit_vertex = [&](ShadedVertex const& sv) {
		ClippedVertex cv;
		float inv_w = 1.0f / sv.clip_position.w;
		cv.fb_position = clip_to_fb_scale * inv_w * sv.clip_position.xyz() + clip_to_fb_offset;
		cv.inv_w = inv_w;
		cv.attributes = sv.attributes;
		clipped_vertices.emplace_back(cv);
	};

	// actually do clipping:
	if constexpr (primitive_type == PrimitiveType::Lines) {
		for (uint32_t i = 0; i + 1 < shaded_vertices.size(); i += 2) {
			clip_line(shaded_vertices[i], shaded_vertices[i + 1], emit_vertex);
		}
	} else if constexpr (primitive_type == PrimitiveType::Triangles) {
		for (uint32_t i = 0; i + 2 < shaded_vertices.size(); i += 3) {
			clip_triangle(shaded_vertices[i], shaded_vertices[i + 1], shaded_vertices[i + 2], emit_vertex);
		}
	} else {
		static_assert(primitive_type == PrimitiveType::Lines, "Unsupported primitive type.");
	}

	//--------------------------
	// rasterize primitives:

	std::vector<Fragment> fragments;

	// helper used to put output of rasterization functions into fragments:
	auto emit_fragment = [&](Fragment const& f) { fragments.emplace_back(f); };

	// added here
	for (uint32_t s = 0; s < samples.size(); ++s) { 
		// samples[s].
		float sample_one_x = samples[s].x;
		float sample_one_y = samples[s].y;
		// actually do rasterization:
		// std::cout<<sample_one_x<<sample_one_y;
		if constexpr (primitive_type == PrimitiveType::Lines) {
			for (uint32_t i = 0; i + 1 < clipped_vertices.size(); i += 2) {

				clipped_vertices[i].fb_position.x = -0.5f+
				clipped_vertices[i].fb_position.x + sample_one_x;

				clipped_vertices[i].fb_position.y = -0.5f+
				clipped_vertices[i].fb_position.y + sample_one_y;

				clipped_vertices[i+1].fb_position.x = -0.5f+
				clipped_vertices[i+1].fb_position.x + sample_one_x;

				clipped_vertices[i+1].fb_position.y = -0.5f+
				clipped_vertices[i+1].fb_position.y + sample_one_y;

				rasterize_line(clipped_vertices[i], clipped_vertices[i + 1], emit_fragment);
			}
		} else if constexpr (primitive_type == PrimitiveType::Triangles) {
			for (uint32_t i = 0; i + 2 < clipped_vertices.size(); i += 3) {

				clipped_vertices[i].fb_position.x = -0.5f+
				clipped_vertices[i].fb_position.x + sample_one_x;

				clipped_vertices[i].fb_position.y = -0.5f+
				clipped_vertices[i].fb_position.y + sample_one_y;

				clipped_vertices[i+1].fb_position.x = -0.5f+
				clipped_vertices[i+1].fb_position.x + sample_one_x;

				clipped_vertices[i+1].fb_position.y = -0.5f+
				clipped_vertices[i+1].fb_position.y + sample_one_y;

				clipped_vertices[i+2].fb_position.x = -0.5f+
				clipped_vertices[i+2].fb_position.x + sample_one_x;

				clipped_vertices[i+2].fb_position.y = -0.5f+
				clipped_vertices[i+2].fb_position.y + sample_one_y;

				rasterize_triangle(clipped_vertices[i], clipped_vertices[i + 1], clipped_vertices[i + 2], emit_fragment);
			}
		} else {
			static_assert(primitive_type == PrimitiveType::Lines, "Unsupported primitive type.");
		}

	
		//--------------------------
		// depth test + shade + blend fragments:
		uint32_t out_of_range = 0; // check if rasterization produced fragments outside framebuffer 
								// (indicates something is wrong with clipping)
		for (auto const& f : fragments) {

			// fragment location (in pixels):
			int32_t x = (int32_t)std::floor(f.fb_position.x);
			int32_t y = (int32_t)std::floor(f.fb_position.y);

			// if clipping is working properly, this condition shouldn't be needed;
			// however, it prevents crashes while you are working on your clipping functions,
			// so we suggest leaving it in place:
			if (x < 0 || (uint32_t)x >= framebuffer.width || 
				y < 0 || (uint32_t)y >= framebuffer.height) {
				++out_of_range;
				continue;
			}

			// local names that refer to destination sample in framebuffer:
			float& fb_depth = framebuffer.depth_at(x, y, s);
			Spectrum& fb_color = framebuffer.color_at(x, y, s);


			// depth test:
			if constexpr ((flags & PipelineMask_Depth) == Pipeline_Depth_Always) {
				// "Always" means the depth test always passes.
			} else if constexpr ((flags & PipelineMask_Depth) == Pipeline_Depth_Never) {
				// "Never" means the depth test never passes.
				continue; //discard this fragment
			} else if constexpr ((flags & PipelineMask_Depth) == Pipeline_Depth_Less) {
				// "Less" means the depth test passes when the new fragment has depth less than the stored depth.
				// A1T4: Depth_Less
				// TODO: implement depth test! We want to only emit fragments that have a depth less than the stored depth, hence "Depth_Less".
				if(f.fb_position.z>=fb_depth) continue;
			} else {
				static_assert((flags & PipelineMask_Depth) <= Pipeline_Depth_Always, "Unknown depth test flag.");
			}

			// if depth test passes, and depth writes aren't disabled, write depth to depth buffer:
			if constexpr (!(flags & Pipeline_DepthWriteDisableBit)) {
				fb_depth = f.fb_position.z;
			}

			// shade fragment:
			ShadedFragment sf;
			sf.fb_position = f.fb_position;
			Program::shade_fragment(parameters, f.attributes, f.derivatives, &sf.color, &sf.opacity);

			// write color to framebuffer if color writes aren't disabled:
			if constexpr (!(flags & Pipeline_ColorWriteDisableBit)) {
				// blend fragment:
				if constexpr ((flags & PipelineMask_Blend) == Pipeline_Blend_Replace) {
					fb_color = sf.color;
				} else if constexpr ((flags & PipelineMask_Blend) == Pipeline_Blend_Add) {
					// A1T4: Blend_Add
					// TODO: framebuffer color should have fragment color multiplied by fragment opacity added to it.
					fb_color += sf.color*sf.opacity; //<-- replace this line
				} else if constexpr ((flags & PipelineMask_Blend) == Pipeline_Blend_Over) {
					// A1T4: Blend_Over
					// TODO: set framebuffer color to the result of "over" blending (also called "alpha blending") the fragment color over the framebuffer color, using the fragment's opacity
					// 		 You may assume that the framebuffer color has its alpha premultiplied already, and you just want to compute the resulting composite color
					fb_color = sf.color*sf.opacity + (1-sf.opacity)*fb_color; //<-- replace this line
				} else {
					static_assert((flags & PipelineMask_Blend) <= Pipeline_Blend_Over, "Unknown blending flag.");
				}
			}
		}
		
		if (out_of_range > 0) {
			if constexpr (primitive_type == PrimitiveType::Lines) {
				warn("Produced %d fragments outside framebuffer; this indicates something is likely "
					"wrong with the clip_line function.",
					out_of_range);
			} else if constexpr (primitive_type == PrimitiveType::Triangles) {
				warn("Produced %d fragments outside framebuffer; this indicates something is likely "
					"wrong with the clip_triangle function.",
					out_of_range);
			}
		}
	}
	

}

// -------------------------------------------------------------------------
// clipping functions

// helper to interpolate between vertices:
template<PrimitiveType p, class P, uint32_t F>
auto Pipeline<p, P, F>::lerp(ShadedVertex const& a, ShadedVertex const& b, float t) -> ShadedVertex {
	ShadedVertex ret;
	ret.clip_position = (b.clip_position - a.clip_position) * t + a.clip_position;
	for (uint32_t i = 0; i < ret.attributes.size(); ++i) {
		ret.attributes[i] = (b.attributes[i] - a.attributes[i]) * t + a.attributes[i];
	}
	return ret;
}

/*
 * clip_line - clip line to portion with -w <= x,y,z <= w, emit vertices of clipped line (if non-empty)
 *  	va, vb: endpoints of line
 *  	emit_vertex: call to produce truncated line
 *
 * If clipping shortens the line, attributes of the shortened line should respect the pipeline's interpolation mode.
 * 
 * If no portion of the line remains after clipping, emit_vertex will not be called.
 *
 * The clipped line should have the same direction as the full line.
 */
template<PrimitiveType p, class P, uint32_t flags>
void Pipeline<p, P, flags>::clip_line(ShadedVertex const& va, ShadedVertex const& vb,
                                      std::function<void(ShadedVertex const&)> const& emit_vertex) {
	// Determine portion of line over which:
	// 		pt = (b-a) * t + a
	//  	-pt.w <= pt.x <= pt.w
	//  	-pt.w <= pt.y <= pt.w
	//  	-pt.w <= pt.z <= pt.w
	// ... as a range [min_t, max_t]:

	float min_t = 0.0f;
	float max_t = 1.0f;

	// want to set range of t for a bunch of equations like:
	//    a.x + t * ba.x <= a.w + t * ba.w
	// so here's a helper:
	auto clip_range = [&min_t, &max_t](float l, float dl, float r, float dr) {
		// restrict range such that:
		// l + t * dl <= r + t * dr
		// re-arranging:
		//  l - r <= t * (dr - dl)
		if (dr == dl) {
			// want: l - r <= 0
			if (l - r > 0.0f) {
				// works for none of range, so make range empty:
				min_t = 1.0f;
				max_t = 0.0f;
			}
		} else if (dr > dl) {
			// since dr - dl is positive:
			// want: (l - r) / (dr - dl) <= t
			min_t = std::max(min_t, (l - r) / (dr - dl));
		} else { // dr < dl
			// since dr - dl is negative:
			// want: (l - r) / (dr - dl) >= t
			max_t = std::min(max_t, (l - r) / (dr - dl));
		}
	};

	// local names for clip positions and their difference:
	Vec4 const& a = va.clip_position;
	Vec4 const& b = vb.clip_position;
	Vec4 const ba = b - a;

	// -a.w - t * ba.w <= a.x + t * ba.x <= a.w + t * ba.w
	clip_range(-a.w, -ba.w, a.x, ba.x);
	clip_range(a.x, ba.x, a.w, ba.w);
	// -a.w - t * ba.w <= a.y + t * ba.y <= a.w + t * ba.w
	clip_range(-a.w, -ba.w, a.y, ba.y);
	clip_range(a.y, ba.y, a.w, ba.w);
	// -a.w - t * ba.w <= a.z + t * ba.z <= a.w + t * ba.w
	clip_range(-a.w, -ba.w, a.z, ba.z);
	clip_range(a.z, ba.z, a.w, ba.w);

	if (min_t < max_t) {
		if (min_t == 0.0f) {
			emit_vertex(va);
		} else {
			ShadedVertex out = lerp(va, vb, min_t);
			// don't interpolate attributes if in flat shading mode:
			if constexpr ((flags & PipelineMask_Interp) == Pipeline_Interp_Flat) {
				out.attributes = va.attributes;
			}
			emit_vertex(out);
		}
		if (max_t == 1.0f) {
			emit_vertex(vb);
		} else {
			ShadedVertex out = lerp(va, vb, max_t);
			// don't interpolate attributes if in flat shading mode:
			if constexpr ((flags & PipelineMask_Interp) == Pipeline_Interp_Flat) {
				out.attributes = va.attributes;
			}
			emit_vertex(out);
		}
	}
}

/*
 * clip_triangle - clip triangle to portion with -w <= x,y,z <= w, emit resulting shape as triangles (if non-empty)
 *  	va, vb, vc: vertices of triangle
 *  	emit_vertex: call to produce clipped triangles (three calls per triangle)
 *
 * If clipping truncates the triangle, attributes of the new vertices should respect the pipeline's interpolation mode.
 * 
 * If no portion of the triangle remains after clipping, emit_vertex will not be called.
 *
 * The clipped triangle(s) should have the same winding order as the full triangle.
 */
template<PrimitiveType p, class P, uint32_t flags>
void Pipeline<p, P, flags>::clip_triangle(
	ShadedVertex const& va, ShadedVertex const& vb, ShadedVertex const& vc,
	std::function<void(ShadedVertex const&)> const& emit_vertex) {
	// A1EC: clip_triangle
	// TODO: correct code!
	// print(w);
	emit_vertex(va);
	emit_vertex(vb);
	emit_vertex(vc);
}

// -------------------------------------------------------------------------
// rasterization functions

/*
 * rasterize_line:
 * calls emit_fragment( frag ) for every pixel "covered" by the line (va.fb_position.xy, vb.fb_position.xy).
 *
 *    a pixel (x,y) is "covered" by the line if it exits the inscribed diamond:
 * 
 *        (x+0.5,y+1)
 *        /        \
 *    (x,y+0.5)  (x+1,y+0.5)
 *        \        /
 *         (x+0.5,y)
 *
 *    to avoid ambiguity, we consider diamonds to contain their left and bottom points
 *    but not their top and right points. 
 * 
 * 	  since 45 degree lines breaks this rule, our rule in general is to rasterize the line as if its
 *    endpoints va and vb were at va + (e, e^2) and vb + (e, e^2) where no smaller nonzero e produces 
 *    a different rasterization result. 
 *    We will not explicitly check for 45 degree lines along the diamond edges (this will be extra credit),
 *    but you should be able to handle 45 degree lines in every other case (such as starting from pixel centers)
 *
 * for each such diamond, pass Fragment frag to emit_fragment, with:
 *  - frag.fb_position.xy set to the center (x+0.5,y+0.5)
 *  - frag.fb_position.z interpolated linearly between va.fb_position.z and vb.fb_position.z
 *  - frag.attributes set to va.attributes (line will only be used in Interp_Flat mode)
 *  - frag.derivatives set to all (0,0)
 *
 * when interpolating the depth (z) for the fragments, you may use any depth the line takes within the pixel
 * (i.e., you don't need to interpolate to, say, the closest point to the pixel center)
 *
 * If you wish to work in fixed point, check framebuffer.h for useful information about the framebuffer's dimensions.
 */
template<PrimitiveType p, class P, uint32_t flags>
void Pipeline<p, P, flags>::rasterize_line(
	ClippedVertex const& va, ClippedVertex const& vb,
	std::function<void(Fragment const&)> const& emit_fragment) {
	if constexpr ((flags & PipelineMask_Interp) != Pipeline_Interp_Flat) {
		assert(0 && "rasterize_line should only be invoked in flat interpolation mode.");
	}
	// A1T2: rasterize_line

	// TODO: Check out the block comment above this function for more information on how to fill in
	// this function!
	// The OpenGL specification section 3.5 may also come in handy.
	// do same logic in slide
	float x1 = va.fb_position.x;
	float y1 = va.fb_position.y;
	float x2 = vb.fb_position.x;
	float y2 = vb.fb_position.y;
	// check for delta
	//i=0 -> x, a.x<b.x
	float delta_x = std::abs(x2-x1);
	float delta_y = std::abs(y2-y1);
	bool x_y_flag = true; // true: delta x > delta t; false: delta y > delta x
	if(delta_x<delta_y){
		// i->1,j->0
		// same as swap x,y
		float tmp = x1;
		x1=y1;
		y1=tmp;
		tmp = x2;
		x2=y2;
		y2=tmp;
		x_y_flag = false;
	}

	// check ai,bi
	if(x1>x2){
		//swap a,b
		float tmp = x1;
		x1=x2;
		x2=tmp;
		tmp = y1;
		y1=y2;
		y2=tmp;
	}
	int t1 = (int)std::floor(x1);
	int t2 = (int)std::floor(x2);
	
	// for u from t1 to t2
	// end point considered independently
	for(int u = t1+1;u<t2;u++){
		float w = ((u+0.5f)-x1)/(x2-x1);
		float v = ((float)w*(y2-y1) + y1)*1.0f;
		// give value
		Fragment mid;
		if(x_y_flag){
			mid.fb_position.x=(float)std::floor(u)+0.5f;
			mid.fb_position.y=std::floor(v)+0.5f;
		}else{
			mid.fb_position.y=(float)std::floor(u)+0.5f;
			mid.fb_position.x=std::floor(v)+0.5f;
		}
		//= (va.fb_position + vb.fb_position) / 2.0f;
		// std::cout<<v<<" "<<u<<",";
		mid.attributes = va.attributes;
		mid.derivatives.fill(Vec2(0.0f, 0.0f));
		emit_fragment(mid);
	}


	// sepcial case:last one
	float end_diamond_x,end_diamond_y;
	// check which is end point
	end_diamond_x = (float)std::floor(vb.fb_position.x);
	end_diamond_y = (float)std::floor(vb.fb_position.y);
	Vec3 end1 = Vec3(end_diamond_x+0.5f,end_diamond_y,0.0f);
	Vec3 end2 = Vec3(end_diamond_x+1.0f,end_diamond_y+0.5f,0.0f);
	Vec3 end3 = Vec3(end_diamond_x+0.5f,end_diamond_y+1.0f,0.0f);
	Vec3 end4 = Vec3(end_diamond_x,end_diamond_y+0.5f,0.0f);
	int intersect_count=0;
	// check cross
	if(check_line_intersect(va.fb_position,vb.fb_position,end1,end2))
		intersect_count++;
	if(check_line_intersect(va.fb_position,vb.fb_position,end2,end3))
		intersect_count++;
	if(check_line_intersect(va.fb_position,vb.fb_position,end3,end4))
		intersect_count++;
	// special case, explicitly check on the endpoints for diamond-exit rule
	if(intersect_count>=1 ){
		//corner cases: not exit
		// also for 45 degrees
		if(!(vb.fb_position.x==end_diamond_x+0.5f) || 
		!(vb.fb_position.y==end_diamond_y) || 
		!(std::abs((vb.fb_position.y-va.fb_position.y)/(vb.fb_position.x-va.fb_position.x))<0.5)){

			// std::cout<<"here";
			// last one
			Fragment end_f;
			end_f.fb_position.x=end_diamond_x+0.5f;
			end_f.fb_position.y=end_diamond_y+0.5f;
			//= (va.fb_position + vb.fb_position) / 2.0f;
			end_f.attributes = va.attributes;
			end_f.derivatives.fill(Vec2(0.0f, 0.0f));
			emit_fragment(end_f);
			// std::cout<<end_diamond_x<<" "<<end_diamond_y<<",";
			// std::cout<<intersect_count<<" ";
		}
	}
		
	// As a placeholder, draw a point in the middle of the line:
	//(remove this code once you have a real implementation)

	// sepcial case:last one
	float end_diamond_x2,end_diamond_y2;
	end_diamond_x2 = (float)std::floor(va.fb_position.x);
	end_diamond_y2 = (float)std::floor(va.fb_position.y);
	end1 = Vec3(end_diamond_x2+0.5f,end_diamond_y2,0.0f);
	end2 = Vec3(end_diamond_x2+1.0f,end_diamond_y2+0.5f,0.0f);
	end3 = Vec3(end_diamond_x2+0.5f,end_diamond_y2+1.0f,0.0f);
	end4 = Vec3(end_diamond_x2,end_diamond_y2+0.5f,0.0f);
	intersect_count=0;
	// avoid duplicate
	if(((int)end_diamond_x2==(int)end_diamond_x)
	 && ((int)end_diamond_y2==(int)end_diamond_y)){
		return;
	}
	if(check_line_intersect(va.fb_position,vb.fb_position,end1,end2))
		intersect_count++;
	if(check_line_intersect(va.fb_position,vb.fb_position,end2,end3))
		intersect_count++;
	if(check_line_intersect(va.fb_position,vb.fb_position,end3,end4))
		intersect_count++;
	// special case, explicitly check on the endpoints for diamond-exit rule
	if(intersect_count>=1){
		// return;
		// std::cout<<"here";
		if(!(va.fb_position.x==end_diamond_x+0.5f) || 
		!(va.fb_position.y==end_diamond_y) || 
		!(std::abs((vb.fb_position.y-va.fb_position.y)/(vb.fb_position.x-va.fb_position.x))<0.5)){
			Fragment end_f2;
			end_f2.fb_position.x=end_diamond_x2+0.5f;
			end_f2.fb_position.y=end_diamond_y2+0.5f;
			//= (va.fb_position + vb.fb_position) / 2.0f;
			end_f2.attributes = va.attributes;
			end_f2.derivatives.fill(Vec2(0.0f, 0.0f));
			emit_fragment(end_f2);
			// std::cout<<end_diamond_x2<<" "<<end_diamond_y2<<",";
		}
	}
	
}

template<PrimitiveType p, class P, uint32_t flags>
bool Pipeline<p, P, flags>::check_point_intersect(Vec3 A, Vec3 B, Vec3 C) {
	if((C.y-A.y) * (B.x-A.x) > (B.y-A.y) * (C.x-A.x))
		return true;
    return false;
}

// Return true if line segments AB and CD intersect
template<PrimitiveType p, class P, uint32_t flags>
bool Pipeline<p, P, flags>::check_line_intersect(Vec3 A, Vec3 B, Vec3 C, Vec3 D) {
	if((check_point_intersect(A,C,D) != check_point_intersect(B,C,D)) &&
	(check_point_intersect(A,B,C) != check_point_intersect(A,B,D)))
		return true;
	return false;
}
/*
 * rasterize_triangle(a,b,c,emit) calls 'emit(frag)' at every location
 *  	(x+0.5,y+0.5) (where x,y are integers) covered by triangle (a,b,c).
 *
 * The emitted fragment should have:
 * - frag.fb_position.xy = (x+0.5, y+0.5)
 * - frag.fb_position.z = linearly interpolated fb_position.z from a,b,c (NOTE: does not depend on Interp mode!)
 * - frag.attributes = depends on Interp_* flag in flags:
 *   - if Interp_Flat: copy from va.attributes
 *   - if Interp_Smooth: interpolate as if (a,b,c) is a 2D triangle flat on the screen
 *   - if Interp_Correct: use perspective-correct interpolation
 * - frag.derivatives = derivatives w.r.t. fb_position.x and fb_position.y of the first frag.derivatives.size() attributes.
 *
 * Notes on derivatives:
 * 	The derivatives are partial derivatives w.r.t. screen locations. That is:
 *    derivatives[i].x = d/d(fb_position.x) attributes[i]
 *    derivatives[i].y = d/d(fb_position.y) attributes[i]
 *  You may compute these derivatives analytically or numerically.
 *
 *  See section 8.12.1 "Derivative Functions" of the GLSL 4.20 specification for some inspiration. (*HOWEVER*, the spec is solving a harder problem, and also nothing in the spec is binding on your implementation)
 *
 *  One approach is to rasterize blocks of four fragments and use forward and backward differences to compute derivatives.
 *  To assist you in this approach, keep in mind that the framebuffer size is *guaranteed* to be even. (see framebuffer.h)
 *
 * Notes on coverage:
 *  If two triangles are on opposite sides of the same edge, and a
 *  fragment center lies on that edge, rasterize_triangle should
 *  make sure that exactly one of the triangles emits that fragment.
 *  (Otherwise, speckles or cracks can appear in the final render.)
 * 
 *  For degenerate (co-linear) triangles, you may consider them to not be on any side of an edge.
 * 	Thus, even if two degnerate triangles share an edge that contains a fragment center, you don't need to emit it.
 *  You will not lose points for doing something reasonable when handling this case
 *
 *  This is pretty tricky to get exactly right!
 *
 */
template<PrimitiveType p, class P, uint32_t flags>
void Pipeline<p, P, flags>::rasterize_triangle(
	ClippedVertex const& va, ClippedVertex const& vb, ClippedVertex const& vc,
	std::function<void(Fragment const&)> const& emit_fragment) {
	// NOTE: it is okay to restructure this function to allow these tasks to use the
	//  same code paths. Be aware, however, that all of them need to remain working!
	//  (e.g., if you break Flat while implementing Correct, you won't get points
	//   for Flat.)
	if constexpr ((flags & PipelineMask_Interp) == Pipeline_Interp_Flat) {
		// A1T3: flat triangles
		// TODO: rasterize triangle (see block comment above this function).
		// Vec3 v0v1 = va.fb_position - vb.fb_position;
		// tests to samples that lie within a screen-space bounding box of the triangle
		float min_x = std::min(va.fb_position.x,vb.fb_position.x);
		min_x = std::min(vc.fb_position.x,min_x);
		float max_x = std::max(va.fb_position.x,vb.fb_position.x);
		max_x = std::max(vc.fb_position.x,max_x);
		float min_y = std::min(va.fb_position.y,vb.fb_position.y);
		min_y = std::min(vc.fb_position.y,min_y);
		float max_y = std::max(va.fb_position.y,vb.fb_position.y);
		max_y = std::max(vc.fb_position.y,max_y);
		int x_start = (int)std::floor(min_x);
		int x_end = (int)std::floor(max_x);
		int y_start = (int)std::floor(min_y);
		int y_end = (int)std::floor(max_y);

		Vec3 a = va.fb_position;
		Vec3 b = vb.fb_position;
		Vec3 c = vc.fb_position;

		for(int i=x_start;i<=x_end;i++){
			for(int j=y_start;j<=y_end;j++){
				Vec3 q = Vec3(i*1.0f+0.5f,j*1.0f+0.5f,0.0f);
				// check if inside and if on edge: use top-left rules
				if(!check_inside(a,b,c,q) && !check_on_edge_top_left(a,b,c,q,min_x,max_y))
					continue;
				Fragment mid;
				mid.fb_position.x=i+0.5f;
				mid.fb_position.y=j+0.5f;
				mid.fb_position.z=va.fb_position.z;
				mid.attributes = va.attributes;
				mid.derivatives.fill(Vec2(0.0f, 0.0f));
				emit_fragment(mid);
			}
		}
	} else if constexpr ((flags & PipelineMask_Interp) == Pipeline_Interp_Smooth) {
		// A1T5: screen-space smooth triangles
		// TODO: rasterize triangle (see block comment above this function).
		float min_x = std::min(va.fb_position.x,vb.fb_position.x);
		min_x = std::min(vc.fb_position.x,min_x);
		float max_x = std::max(va.fb_position.x,vb.fb_position.x);
		max_x = std::max(vc.fb_position.x,max_x);
		float min_y = std::min(va.fb_position.y,vb.fb_position.y);
		min_y = std::min(vc.fb_position.y,min_y);
		float max_y = std::max(va.fb_position.y,vb.fb_position.y);
		max_y = std::max(vc.fb_position.y,max_y);
		int x_start = (int)std::floor(min_x);
		int x_end = (int)std::floor(max_x);
		int y_start = (int)std::floor(min_y);
		int y_end = (int)std::floor(max_y);

		Vec3 a = va.fb_position;
		Vec3 b = vb.fb_position;
		Vec3 c = vc.fb_position;

		for(int i=x_start;i<=x_end;i++){
			for(int j=y_start;j<=y_end;j++){
				Vec3 q = Vec3(i*1.0f+0.5f,j*1.0f+0.5f,0.0f);
				float wa=0;
				float wb=0;
				float wc=0;
				if(!check_inside_and_ratio(a,b,c,q,wa,wb,wc)&& !check_on_edge_top_left(a,b,c,q,min_x,max_y))
					continue;
				Fragment mid;
				mid.fb_position.x=i+0.5f;
				mid.fb_position.y=j+0.5f;

				mid.fb_position.z=va.fb_position.z * wa +
				vb.fb_position.z * wb+vc.fb_position.z * wc;

				for(uint32_t k=0;k<va.attributes.size();k++){
					mid.attributes[k]= va.attributes[k]
					* wa +vb.attributes[k] * wb +
					vc.attributes[k] * wc;
				}

				// for derivates:
				// x direction:x+0.1
				// or x-0.1 if x+0.1 not work
				float wa_x,wb_x,wc_x;
				bool derivate_x_flag = true;
				// delta_x:0.1
				Vec3 qx = Vec3(i*1.0f+0.4f,j*1.0f+0.5f,0.0f);
				// get_ratio_outOfRange(a,b,c,qx,wa_x,wb_x,wc_x);
				check_inside_and_ratio(a,b,c,qx,wa_x,wb_x,wc_x);

				if(wa_x+wb_x+wc_x!=1){
					qx = Vec3(i*1.0f+0.6f,j*1.0f+0.5f,0.0f);
					check_inside_and_ratio(a,b,c,qx,wa_x,wb_x,wc_x);
					derivate_x_flag =false;
				}
				Fragment mid_x;
				for(uint32_t k=0;k<va.attributes.size();k++){
					mid_x.attributes[k]= va.attributes[k]
					* wa_x +vb.attributes[k] * wb_x +
					vc.attributes[k] * wc_x;
				}





				// std::cout<<wa<<wb<<wc<<std::endl;
				// y direction:y+0.1
				float wa_y,wb_y,wc_y;
				bool derivate_y_flag = true;
				// delat y:0.1
				Vec3 qy = Vec3(i*1.0f+0.5f,j*1.0f+0.4f,0.0f);
				// get_ratio_outOfRange(a,b,c,qy,wa_y,wb_y,wc_y);
				check_inside_and_ratio(a,b,c,qy,wa_y,wb_y,wc_y);
				if(wa_y+wb_y+wc_y!=1){
					qy = Vec3(i*1.0f+0.5f,j*1.0f+0.6f,0.0f);
					check_inside_and_ratio(a,b,c,qy,wa_y,wb_y,wc_y);
					derivate_y_flag =false;
				}
				Fragment mid_y;
				for(uint32_t k=0;k<va.attributes.size();k++){
					mid_y.attributes[k]= va.attributes[k]
					* wa_y +vb.attributes[k] * wb_y +
					vc.attributes[k] * wc_y;
				}
				

				// std::cout<<wa_y<<wb_y<<wc_y<<std::endl;
				// mid.derivatives.fill(Vec2(0.0f, 0.0f));
				for(uint32_t k=0;k<mid.derivatives.size();k++){
					float derivate_x;
					float derivate_y;
					if(derivate_x_flag)
						derivate_x = -(mid_x.attributes[k]-mid.attributes[k])*10;
					else
						derivate_x = (mid_x.attributes[k]-mid.attributes[k])*10;
					if(derivate_y_flag)
						derivate_y = -(mid_y.attributes[k]-mid.attributes[k])*10;
					else
						derivate_y = (mid_y.attributes[k]-mid.attributes[k])*10;
					mid.derivatives[k] = Vec2(derivate_x, derivate_y);					
				}
				
				emit_fragment(mid);
			}
		}
		
		// Pipeline<PrimitiveType::Lines, P, (flags & ~PipelineMask_Interp) | Pipeline_Interp_Flat>::rasterize_triangle(va, vb, vc, emit_fragment);
	} else if constexpr ((flags & PipelineMask_Interp) == Pipeline_Interp_Correct) {
		// A1T5: perspective correct triangles
		// TODO: rasterize triangle (block comment above this function).
		float min_x = std::min(va.fb_position.x,vb.fb_position.x);
		min_x = std::min(vc.fb_position.x,min_x);
		float max_x = std::max(va.fb_position.x,vb.fb_position.x);
		max_x = std::max(vc.fb_position.x,max_x);
		float min_y = std::min(va.fb_position.y,vb.fb_position.y);
		min_y = std::min(vc.fb_position.y,min_y);
		float max_y = std::max(va.fb_position.y,vb.fb_position.y);
		max_y = std::max(vc.fb_position.y,max_y);
		int x_start = (int)std::floor(min_x);
		int x_end = (int)std::floor(max_x);
		int y_start = (int)std::floor(min_y);
		int y_end = (int)std::floor(max_y);

		Vec3 a = va.fb_position;
		Vec3 b = vb.fb_position;
		Vec3 c = vc.fb_position;

		for(int i=x_start;i<=x_end;i++){
			for(int j=y_start;j<=y_end;j++){
				Vec3 q = Vec3(i*1.0f+0.5f,j*1.0f+0.5f,0.0f);
				float wa=0;
				float wb=0;
				float wc=0;
				if(!check_inside_and_ratio(a,b,c,q,wa,wb,wc)&& !check_on_edge_top_left(a,b,c,q,min_x,max_y))
					continue;
				Fragment mid;
				mid.fb_position.x=i+0.5f;
				mid.fb_position.y=j+0.5f;
				mid.fb_position.z=va.fb_position.z;
				mid.attributes = va.attributes;


				// for z
				float z_af = va.inv_w;//va.fb_position.z;
				float z_bf = vb.inv_w;//vb.fb_position.z;
				float z_cf = vc.inv_w;//vc.fb_position.z;

				double sum = z_af* wa +z_bf * wb +z_cf * wc;

				for(uint32_t k=0;k<va.attributes.size();k++){
					double res= va.attributes[k]
					* wa * z_af +vb.attributes[k] * wb *z_bf +
					vc.attributes[k] * wc*z_cf;
					mid.attributes[k]=(float)(res/sum);
				}

				mid.fb_position.z=va.fb_position.z * wa * z_af +
				vb.fb_position.z * wb * z_bf+vc.fb_position.z * wc* z_cf;

				mid.fb_position.z/=(float)sum;

				// for derivates:
				// x direction:x+0.1 or x-0.1
				bool derivate_x_flag = true; // true for x-0.1,false for x+0.1
				float wa_x,wb_x,wc_x;
				// deltax:0.1
				Vec3 qx = Vec3(i*1.0f+0.4f,j*1.0f+0.5f,0.0f);
				check_inside_and_ratio(a,b,c,qx,wa_x,wb_x,wc_x);

				if(wa_x+wb_x+wc_x!=1){
					qx = Vec3(i*1.0f+0.6f,j*1.0f+0.5f,0.0f);
					check_inside_and_ratio(a,b,c,qx,wa_x,wb_x,wc_x);
					derivate_x_flag =false;
				}

				double sum_x = z_af* wa_x +z_bf * wb_x +z_cf * wc_x;
				Fragment mid_x;
				for(uint32_t k=0;k<va.attributes.size();k++){
					double res= va.attributes[k]
					* wa_x * z_af +vb.attributes[k] * wb_x *z_bf +
					vc.attributes[k] * wc_x*z_cf;
					mid_x.attributes[k]=(float)(res/sum_x);
				}



				// y direction:y+0.1
				float wa_y,wb_y,wc_y;
				bool derivate_y_flag = true;
				// deltay:0.1
				Vec3 qy = Vec3(i*1.0f+0.5f,j*1.0f+0.4f,0.0f);
				check_inside_and_ratio(a,b,c,qy,wa_y,wb_y,wc_y);
				if(wa_y+wb_y+wc_y!=1){
					qy = Vec3(i*1.0f+0.5f,j*1.0f+0.6f,0.0f);
					check_inside_and_ratio(a,b,c,qy,wa_y,wb_y,wc_y);
					derivate_y_flag =false;
				}

				Fragment mid_y;
				double sum_y = z_af* wa_y +z_bf * wb_y +z_cf * wc_y;
				for(uint32_t k=0;k<va.attributes.size();k++){
					double res= va.attributes[k]
					* wa_y * z_af +vb.attributes[k] * wb_y *z_bf +
					vc.attributes[k] * wc_y*z_cf;
					mid_y.attributes[k]=(float)(res/sum_y);
				}


				// mid.derivatives.fill(Vec2(0.0f, 0.0f));
				for(uint32_t k=0;k<mid.derivatives.size();k++){
					float derivate_x;
					float derivate_y;
					if(derivate_x_flag)
						derivate_x = -(mid_x.attributes[k]-mid.attributes[k])*10;
					else
						derivate_x = (mid_x.attributes[k]-mid.attributes[k])*10;
					if(derivate_y_flag)
						derivate_y = -(mid_y.attributes[k]-mid.attributes[k])*10;
					else
						derivate_y = (mid_y.attributes[k]-mid.attributes[k])*10;
					mid.derivatives[k] = Vec2(derivate_x,derivate_y);
					
				}
				emit_fragment(mid);
			}
		}
		// As a placeholder, here's code that calls the Screen-space interpolation function:
		//(remove this and replace it with a real solution)
		// Pipeline<PrimitiveType::Lines, P, (flags & ~PipelineMask_Interp) | Pipeline_Interp_Smooth>::rasterize_triangle(va, vb, vc, emit_fragment);
	}
}
template<PrimitiveType p, class P, uint32_t flags>
float Pipeline<p, P, flags>::cross_product_2d(Vec3 a,Vec3 b){
	return a.x * b.y - b.x * a.y;
}
template<PrimitiveType p, class P, uint32_t flags>
bool Pipeline<p, P, flags>::check_inside(Vec3 a,Vec3 b,Vec3 c, Vec3 q){
	Vec3 ac = c-a; 
	Vec3 ab = b-a; 
	Vec3 aq = q-a; 

	Vec3 cb = b-c; 
	Vec3 ca = a-c; 
	Vec3 cq = q-c; 

	Vec3 ba = a-b; 
	Vec3 bc = c-b; 
	Vec3 bq = q-b; 
	if( ((cross_product_2d(ac,ab)*cross_product_2d(ac,aq)<0) &&
	(cross_product_2d(cb,ca)*cross_product_2d(cb,cq)<0) &&
	(cross_product_2d(ba,bc)*cross_product_2d(ba,bq)<0)
	) ||
	((cross_product_2d(ac,ab)*cross_product_2d(ac,aq)>0) &&
	(cross_product_2d(cb,ca)*cross_product_2d(cb,cq)>0) &&
	(cross_product_2d(ba,bc)*cross_product_2d(ba,bq)>0)
	)
	){
		return true;
	}
	return false;
}

template<PrimitiveType p, class P, uint32_t flags>
bool Pipeline<p, P, flags>::check_on_edge_top_left(Vec3 a,Vec3 b,Vec3 c, Vec3 q,float min_x
,float max_y){
	Vec3 ac = c-a; 
	Vec3 ab = b-a; 
	Vec3 aq = q-a; 

	Vec3 cb = b-c; 
	Vec3 ca = a-c; 
	Vec3 cq = q-c; 

	Vec3 ba = a-b; 
	Vec3 bc = c-b; 
	Vec3 bq = q-b; 
	if((cross_product_2d(ac,aq)==0)){
		// check on edge
		if(((q.x>=a.x && q.x<=c.x) || (q.x<=a.x && q.x>=c.x)) &&
		((q.y>=a.y && q.y<=c.y) || (q.y<=a.y && q.y>=c.y))){
			// top left rule!
			// if(((a.x==min_x) || (c.x==min_x)) && ((a.y==max_y) || (c.y==max_y)))
			// 	return true;
			//A top edge, is an edge that is exactly horizontal and is above the other edges.
			if(c.y==a.y){
				if(b.y<a.y)
					return true;
			}else{
				// A left edge, is an edge that is not exactly horizontal 
				// and is on the left side of the triangle. A triangle can have one or two left edges
				if(c.x==min_x){
					if(a.y>c.y)
						return true;
				}
				if(a.x==min_x){
					if(c.y>a.y)
						return true;
				}
			}
		}
	}
	if((cross_product_2d(ab,aq)==0)){
		// check on edge
		if(((q.x>=a.x && q.x<=b.x) || (q.x<=a.x && q.x>=b.x))&&
		((q.y>=a.y && q.y<=b.y) || (q.y<=a.y && q.y>=b.y))){
			// top left rule!
			// if(((a.x==min_x) || (b.x==min_x)) && ((a.y==max_y) || (b.y==max_y)))
			// 	return true;
			//A top edge, is an edge that is exactly horizontal and is above the other edges.
			if(b.y==a.y){
				if(c.y<b.y)
					return true;
			}else{
				// A left edge, is an edge that is not exactly horizontal 
				// and is on the left side of the triangle. A triangle can have one or two left edges
				if(b.x==min_x){
					return true;
				}
				if(a.x==min_x){
					return true;
				}
			}
		}
	}
	// std::cout<<q.x<<" "<<q.y<<" "<<cross_product_2d(bc,bq)<<", ";
	if((cross_product_2d(bc,bq)==0)){
		// check on edge
		// firstly: inside edge
		if(((q.x>=b.x && q.x<=c.x) || (q.x<=b.x && q.x>=c.x))&&
		((q.y>=b.y && q.y<=c.y) || (q.y<=b.y && q.y>=c.y))){

			//A top edge, is an edge that is exactly horizontal and is above the other edges.
			if(b.y==c.y){
				if(a.y<b.y)
					return true;
			}else{
				// A left edge, is an edge that is not exactly horizontal 
				// and is on the left side of the triangle. A triangle can have one or two left edges
				if(b.x==min_x){
					if(c.y>b.y)
						return true;
				}
				if(c.x==min_x){
					return true;
				}
			}
			
		}
	}
	return false;
	
}

template<PrimitiveType p, class P, uint32_t flags>
bool Pipeline<p, P, flags>::check_inside_and_ratio(Vec3 a,Vec3 b,Vec3 c, Vec3 q,
	float& wa,float& wb,float& wc){
	Vec3 ac = c-a; 
	Vec3 ab = b-a; 
	Vec3 aq = q-a; 

	Vec3 cb = b-c; 
	Vec3 ca = a-c; 
	Vec3 cq = q-c; 

	Vec3 ba = a-b; 
	Vec3 bc = c-b; 
	Vec3 bq = q-b; 

	float whole_area = std::abs(cross_product_2d(ac,ab));
	wa = std::abs(cross_product_2d(bq,cq)/whole_area);
	wb = std::abs(cross_product_2d(aq,cq)/whole_area);
	wc = std::abs(cross_product_2d(bq,aq)/whole_area);
	if( ((cross_product_2d(ac,ab)*cross_product_2d(ac,aq)<0) &&
	(cross_product_2d(cb,ca)*cross_product_2d(cb,cq)<0) &&
	(cross_product_2d(ba,bc)*cross_product_2d(ba,bq)<0)
	) ||
	((cross_product_2d(ac,ab)*cross_product_2d(ac,aq)>0) &&
	(cross_product_2d(cb,ca)*cross_product_2d(cb,cq)>0) &&
	(cross_product_2d(ba,bc)*cross_product_2d(ba,bq)>0)
	)
	){
		return true;
	}
	return false;
}

template<PrimitiveType p, class P, uint32_t flags>
void Pipeline<p, P, flags>::get_ratio_outOfRange(Vec3 a,Vec3 b,Vec3 c, Vec3 q,
	float& wa,float& wb,float& wc){
	Vec3 ac = c-a; 
	Vec3 ab = b-a; 
	Vec3 aq = q-a; 

	Vec3 cb = b-c; 
	Vec3 ca = a-c; 
	Vec3 cq = q-c; 

	Vec3 ba = a-b; 
	Vec3 bc = c-b; 
	Vec3 bq = q-b; 

	float wa_a,wb_a,wc_a;
	float whole_area = std::abs(cross_product_2d(ac,ab));
	wa_a = cross_product_2d(bq,cq)/whole_area;
	wb_a = cross_product_2d(aq,cq)/whole_area;
	wc_a = cross_product_2d(bq,aq)/whole_area;
	if(wa_a>0){
		wa = wa_a;
	}else{
		wa = std::abs(wa_a)+whole_area;
	}
	if(wb_a<0){
		wb = std::abs(wb_a);
	}else{
		wb = std::abs(wb_a)+whole_area;
	}
	if(wc_a<0){
		wc = std::abs(wc_a);
	}else{
		wc = std::abs(wc_a)+whole_area;
	}
	whole_area = wa+wb+wc;
	wa = wa/whole_area;
	wb = wb/whole_area;
	wc = wc/whole_area;
}
//-------------------------------------------------------------------------
// compile instantiations for all programs and blending and testing types:

#include "programs.h"

template struct Pipeline<PrimitiveType::Triangles, Programs::Lambertian,
                         Pipeline_Blend_Replace | Pipeline_Depth_Always | Pipeline_Interp_Flat>;
template struct Pipeline<PrimitiveType::Triangles, Programs::Lambertian,
                         Pipeline_Blend_Replace | Pipeline_Depth_Always | Pipeline_Interp_Smooth>;
template struct Pipeline<PrimitiveType::Triangles, Programs::Lambertian,
                         Pipeline_Blend_Replace | Pipeline_Depth_Always | Pipeline_Interp_Correct>;
template struct Pipeline<PrimitiveType::Triangles, Programs::Lambertian,
                         Pipeline_Blend_Replace | Pipeline_Depth_Never | Pipeline_Interp_Flat>;
template struct Pipeline<PrimitiveType::Triangles, Programs::Lambertian,
                         Pipeline_Blend_Replace | Pipeline_Depth_Never | Pipeline_Interp_Smooth>;
template struct Pipeline<PrimitiveType::Triangles, Programs::Lambertian,
                         Pipeline_Blend_Replace | Pipeline_Depth_Never | Pipeline_Interp_Correct>;
template struct Pipeline<PrimitiveType::Triangles, Programs::Lambertian,
                         Pipeline_Blend_Replace | Pipeline_Depth_Less | Pipeline_Interp_Flat>;
template struct Pipeline<PrimitiveType::Triangles, Programs::Lambertian,
                         Pipeline_Blend_Replace | Pipeline_Depth_Less | Pipeline_Interp_Smooth>;
template struct Pipeline<PrimitiveType::Triangles, Programs::Lambertian,
                         Pipeline_Blend_Replace | Pipeline_Depth_Less | Pipeline_Interp_Correct>;
template struct Pipeline<PrimitiveType::Triangles, Programs::Lambertian,
                         Pipeline_Blend_Add | Pipeline_Depth_Always | Pipeline_Interp_Flat>;
template struct Pipeline<PrimitiveType::Triangles, Programs::Lambertian,
                         Pipeline_Blend_Add | Pipeline_Depth_Always | Pipeline_Interp_Smooth>;
template struct Pipeline<PrimitiveType::Triangles, Programs::Lambertian,
                         Pipeline_Blend_Add | Pipeline_Depth_Always | Pipeline_Interp_Correct>;
template struct Pipeline<PrimitiveType::Triangles, Programs::Lambertian,
                         Pipeline_Blend_Add | Pipeline_Depth_Never | Pipeline_Interp_Flat>;
template struct Pipeline<PrimitiveType::Triangles, Programs::Lambertian,
                         Pipeline_Blend_Add | Pipeline_Depth_Never | Pipeline_Interp_Smooth>;
template struct Pipeline<PrimitiveType::Triangles, Programs::Lambertian,
                         Pipeline_Blend_Add | Pipeline_Depth_Never | Pipeline_Interp_Correct>;
template struct Pipeline<PrimitiveType::Triangles, Programs::Lambertian,
                         Pipeline_Blend_Add | Pipeline_Depth_Less | Pipeline_Interp_Flat>;
template struct Pipeline<PrimitiveType::Triangles, Programs::Lambertian,
                         Pipeline_Blend_Add | Pipeline_Depth_Less | Pipeline_Interp_Smooth>;
template struct Pipeline<PrimitiveType::Triangles, Programs::Lambertian,
                         Pipeline_Blend_Add | Pipeline_Depth_Less | Pipeline_Interp_Correct>;
template struct Pipeline<PrimitiveType::Triangles, Programs::Lambertian,
                         Pipeline_Blend_Over | Pipeline_Depth_Always | Pipeline_Interp_Flat>;
template struct Pipeline<PrimitiveType::Triangles, Programs::Lambertian,
                         Pipeline_Blend_Over | Pipeline_Depth_Always | Pipeline_Interp_Smooth>;
template struct Pipeline<PrimitiveType::Triangles, Programs::Lambertian,
                         Pipeline_Blend_Over | Pipeline_Depth_Always | Pipeline_Interp_Correct>;
template struct Pipeline<PrimitiveType::Triangles, Programs::Lambertian,
                         Pipeline_Blend_Over | Pipeline_Depth_Never | Pipeline_Interp_Flat>;
template struct Pipeline<PrimitiveType::Triangles, Programs::Lambertian,
                         Pipeline_Blend_Over | Pipeline_Depth_Never | Pipeline_Interp_Smooth>;
template struct Pipeline<PrimitiveType::Triangles, Programs::Lambertian,
                         Pipeline_Blend_Over | Pipeline_Depth_Never | Pipeline_Interp_Correct>;
template struct Pipeline<PrimitiveType::Triangles, Programs::Lambertian,
                         Pipeline_Blend_Over | Pipeline_Depth_Less | Pipeline_Interp_Flat>;
template struct Pipeline<PrimitiveType::Triangles, Programs::Lambertian,
                         Pipeline_Blend_Over | Pipeline_Depth_Less | Pipeline_Interp_Smooth>;
template struct Pipeline<PrimitiveType::Triangles, Programs::Lambertian,
                         Pipeline_Blend_Over | Pipeline_Depth_Less | Pipeline_Interp_Correct>;
template struct Pipeline<PrimitiveType::Lines, Programs::Lambertian,
                         Pipeline_Blend_Replace | Pipeline_Depth_Always | Pipeline_Interp_Flat>;
template struct Pipeline<PrimitiveType::Lines, Programs::Lambertian,
                         Pipeline_Blend_Replace | Pipeline_Depth_Never | Pipeline_Interp_Flat>;
template struct Pipeline<PrimitiveType::Lines, Programs::Lambertian,
                         Pipeline_Blend_Replace | Pipeline_Depth_Less | Pipeline_Interp_Flat>;
template struct Pipeline<PrimitiveType::Lines, Programs::Lambertian,
                         Pipeline_Blend_Add | Pipeline_Depth_Always | Pipeline_Interp_Flat>;
template struct Pipeline<PrimitiveType::Lines, Programs::Lambertian,
                         Pipeline_Blend_Add | Pipeline_Depth_Never | Pipeline_Interp_Flat>;
template struct Pipeline<PrimitiveType::Lines, Programs::Lambertian,
                         Pipeline_Blend_Add | Pipeline_Depth_Less | Pipeline_Interp_Flat>;
template struct Pipeline<PrimitiveType::Lines, Programs::Lambertian,
                         Pipeline_Blend_Over | Pipeline_Depth_Always | Pipeline_Interp_Flat>;
template struct Pipeline<PrimitiveType::Lines, Programs::Lambertian,
                         Pipeline_Blend_Over | Pipeline_Depth_Never | Pipeline_Interp_Flat>;
template struct Pipeline<PrimitiveType::Lines, Programs::Lambertian,
                         Pipeline_Blend_Over | Pipeline_Depth_Less | Pipeline_Interp_Flat>;