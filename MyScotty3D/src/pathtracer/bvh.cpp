
#include "bvh.h"
#include "aggregate.h"
#include "instance.h"
#include "tri_mesh.h"

#include <stack>

namespace PT {

struct BVHBuildData {
	BVHBuildData(size_t start, size_t range, size_t dst) : start(start), range(range), node(dst) {
	}
	size_t start; ///< start index into the primitive array
	size_t range; ///< range of index into the primitive array
	size_t node;  ///< address to update
};

struct SAHBucketData {
	BBox bb;          ///< bbox of all primitives
	size_t num_prims; ///< number of primitives in the bucket
};

template<typename Primitive>
void BVH<Primitive>::build(std::vector<Primitive>&& prims, size_t max_leaf_size) {
	//A3T3 - build a bvh


	// Keep these
    nodes.clear();
    primitives = std::move(prims);

    // Construct a BVH from the given vector of primitives and maximum leaf
    // size configuration.
	

	//TODO
	recursion(max_leaf_size,0);

}

template<typename Primitive>
void BVH<Primitive>::recursion(size_t max_leaf_size, size_t node_idx) {
	// root node
	if(node_idx==0){
		BBox box;
		for(const Primitive &p : primitives){
			box.enclose(p.bbox());
		}
		//box, start,size,l,r = 0
		// printf("%f\n",box.surface_area());
		new_node(box,0,n_primitives(),0,0);
	}

	// recurse until few p left
	if(nodes.at(node_idx).size<=max_leaf_size){
		return;
	}

	// start
	size_t start_idx = nodes.at(node_idx).start;
	size_t size = nodes.at(node_idx).size;
	auto it_start = primitives.begin() + start_idx;
	// end
	// printf("%zu,%zu,%zu\n",node_idx,start_idx,start_idx+size);
	auto it_end = primitives.begin() + start_idx + size;
	// sort(v.begin() + l, v.begin() + r + 1); 
	// first, deal with x axis


	// get cost
	// divide with each one
	BBox left_b;
	BBox right_b;
	BBox all_b;

	size_t left_count = 0;
	size_t right_count = 0;
	// size_t all_count = 0;
	float min_cost = FLT_MAX;
	size_t min_idx = 0;

	BBox best_l_box;
	BBox best_r_box;

	for(int axis_all = 0;axis_all<3;axis_all++){
		if(axis_all==0){
			// printf("x axis\n");
			// printf("%f,%f,%f",primitives[0].bbox().center().x,primitives[1].bbox().center().x,
			// primitives[2].bbox().center().x);
			std::sort(it_start, it_end, [](const Primitive &e1, const Primitive &e2){ 
				// increasing
				float e11 = e1.bbox().center().x;
				float e22 = e2.bbox().center().x;
				return e11<e22;
			});
			// check
			// printf("%f,%f,%f\n",primitives[0].bbox().center().x,primitives[1].bbox().center().x,
			// primitives[2].bbox().center().x);
		}else if(axis_all == 1){
			// printf("y axis\n");
			std::sort(it_start, it_end, [](const Primitive &e1, const Primitive &e2){ 
				// increasing
				float e11 = e1.bbox().center().y;
				float e22 = e2.bbox().center().y;
				return e11<e22;
			});
			// check
			// printf("%f,%f,%f\n",primitives[0].bbox().center().y,primitives[1].bbox().center().y,
			// primitives[2].bbox().center().y);
		}else{
			// printf("z axis\n");
			std::sort(it_start, it_end, [](const Primitive &e1, const Primitive &e2){ 
				// increasing
				float e11 = e1.bbox().center().z;
				float e22 = e2.bbox().center().z;
				return e11<e22;
			});
			// check
			// printf("%f,%f,%f\n",primitives[0].bbox().center().z,primitives[1].bbox().center().z,
			// primitives[2].bbox().center().z);
		}

		// get cost
		// printf("%f,%zd\n",min_cost,min_idx);

		if(size<=32){
			for(size_t i = 0;i<size;i++){
				for(size_t j = start_idx;j<start_idx+i;j++){
					BBox one = primitives[j].bbox();
					left_b.enclose(one);
					left_count++;
					all_b.enclose(one);
				}
				for(size_t j = start_idx+i;j<start_idx+size;j++){
					BBox one = primitives[j].bbox();
					right_b.enclose(one);
					right_count++;
					all_b.enclose(one);
				}
				// printf("num in box:%zu,%zu\n",left_count,right_count);

				float cost = left_count * left_b.surface_area()/all_b.surface_area() +
				right_count * right_b.surface_area()/all_b.surface_area();

				if(cost<min_cost){
					min_cost = cost;
					min_idx = i;
					best_l_box = BBox();
					best_l_box.enclose(left_b);
					best_r_box = BBox();
					best_r_box.enclose(right_b);
				}
				// printf("box:%f,%f\n",left_b.surface_area(),right_b.surface_area());
				// printf("cost%f\n",cost);
				// reinitial
				right_count = 0;
				left_count = 0;
				left_b = BBox();
				right_b = BBox();
				all_b = BBox();
			}
			// printf("cost and index:%f,%zd\n",min_cost,min_idx);
		}else{
			// divide into 32 bins
			size_t bin_num = size_t(std::floor(size/32));
			for(size_t ii = 0;ii<32;ii++){
			// for(size_t i = 1;i<size;i++){
				size_t i = bin_num * ii;
				for(size_t j = start_idx;j<start_idx+i;j++){
					BBox one = primitives[j].bbox();
					left_b.enclose(one);
					left_count++;
					all_b.enclose(one);
				}
				for(size_t j = start_idx+i;j<start_idx+size;j++){
					BBox one = primitives[j].bbox();
					right_b.enclose(one);
					right_count++;
					all_b.enclose(one);
				}
				// printf("num in box:%zu,%zu\n",left_count,right_count);

				float cost = left_count * left_b.surface_area()/all_b.surface_area() +
				right_count * right_b.surface_area()/all_b.surface_area();

				if(cost<min_cost){
					min_cost = cost;
					min_idx = i;
					best_l_box = BBox();
					best_l_box.enclose(left_b);
					best_r_box = BBox();
					best_r_box.enclose(right_b);
				}
				// printf("box:%f,%f\n",left_b.surface_area(),right_b.surface_area());
				// printf("cost%f\n",cost);
				// reinitial
				right_count = 0;
				left_count = 0;
				left_b = BBox();
				right_b = BBox();
				all_b = BBox();
			}
			// printf("cost and index:%f,%zd\n",min_cost,min_idx);
		}
	}
	
	// printf("final box:%f,%f\n",best_l_box.surface_area(),best_r_box.surface_area());
	// printf("final:%f,%zd\n",min_cost,min_idx);

	// !!!! index important
	size_t node_l = new_node(best_l_box, start_idx, min_idx, 0, 0);
	nodes.at(node_idx).l = node_l;
	
	size_t node_r=new_node(best_r_box, min_idx+start_idx, size-min_idx, 0, 0);
    nodes.at(node_idx).r = node_r;

	recursion(max_leaf_size,node_l);
	recursion(max_leaf_size,node_r);
}

// template<typename Primitive>
// void BVH<Primitive>::iteration_axis() {
// }

template<typename Primitive> Trace BVH<Primitive>::hit(const Ray& ray) const {
	//A3T3 - traverse your BVH

    // Implement ray - BVH intersection test. A ray intersects
    // with a BVH aggregate if and only if it intersects a primitive in
    // the BVH that is not an aggregate.

    // The starter code simply iterates through all the primitives.
    // Again, remember you can use hit() on any Primitive value.

	//TODO: replace this code with a more efficient traversal:
    Trace ret;
	// ret.distance = FLT_MAX;
    // for(const Primitive& prim : primitives) {
    //     Trace hit = prim.hit(ray);
    //     ret = Trace::min(ret, hit);
    // }
	// printf("size:%zd",nodes.size());
	// empty
	if(nodes.size()==0){
		return ret;
	}
	hit_recursion(ray,0,ret);
	return ret;
	// printf("final,%f,%d",ret.distance,ret.hit);
	// return ret;
}

// template<typename Primitive> Trace BVH<Primitive>::hit_recursion(const Ray& ray, size_t node_idx,Trace ret) const {

// 	// needed????
// 	// Vec2 curr_times=Vec2(0.0f,FLT_MAX);
// 	// bool if_hit = nodes.at(node_idx).bbox.hit(ray,curr_times);
// 	// if(!if_hit){
// 	// 	return ret;
// 	// }
// 	// printf("node:%zd\n",node_idx);
// 	// printf("%zd,%zd\n",nodes.at(node_idx).size,nodes.at(node_idx).start);

// 	if(nodes.at(node_idx).is_leaf()) {
// 		// printf("leaf\n");
//        	size_t start_idx = nodes.at(node_idx).start;
// 		size_t size = nodes.at(node_idx).size;	
//         for(size_t i = start_idx; i < start_idx+size; i++) {
//             Trace hit = primitives[i].hit(ray);
// 			if(hit.hit)
//             	ret = Trace::min(ret, hit);
// 			// printf("idx:%zd",i);
// 			// printf("%f,%d\n",hit.distance,hit.hit);
// 			// printf("%f,%d\n",ret.distance,ret.hit);
//         }
// 		return ret;
//     } else {
//         size_t left_idx = nodes.at(node_idx).l;
//         size_t right_idx = nodes.at(node_idx).r;

// 		// correct?
// 		Vec2 left_times = Vec2(0.0f,FLT_MAX);
// 		Vec2 right_times = Vec2(0.0f,FLT_MAX);

// 		bool left_h = nodes.at(left_idx).bbox.hit(ray,left_times);
// 		bool right_h = nodes.at(right_idx).bbox.hit(ray,right_times);

// 		// if(!left_h && !right_h){
// 		// 	return ret;
// 		// }else if(left_h){
// 		// 	ret = Trace::min(ret, hit_recursion(ray,left_idx,ret));
// 		// 	return ret;
// 		// }else if(right_h){
// 		// 	ret = Trace::min(ret, hit_recursion(ray,right_idx,ret));
// 		// 	return ret;
// 		// }
// 		// printf("b%f\n",left_times.x);

// 		size_t hit_first = (left_times.x <= right_times.x)?left_idx:right_idx;
// 		size_t hit_second = (left_times.x <= right_times.x)?right_idx:left_idx;


// 		// bool hit_f_ = (left_times.x <= right_times.x)?left_h:right_h;
// 		// bool hit_s_ = (left_times.x <= right_times.x)?right_h:left_h;

// 		float hit_second_t = (left_times.x <= right_times.x)?right_times.x:left_times.x;

// 		// BBox b
// 		(void) hit_first;
// 		(void) hit_second;
// 		(void) left_h;
// 		(void) right_h;
		
// 		(void) hit_second_t;

// 		// if(hit_f_){
// 		ret = Trace::min(ret, hit_recursion(ray,hit_first,ret));
// 		// }
// 		// printf("%f,%f",ret.distance,hit_second_t);

// 		// if(hit_s_){
// 		if(!ret.hit){
// 			ret = Trace::min(ret,hit_recursion(ray,hit_second,ret));
// 		}
// 		else if(hit_second_t <= ret.distance) {
// 			ret = Trace::min(ret,hit_recursion(ray,hit_second,ret));
// 		}
// 		// }


// 		// iteration all, for check if building is right
// 		// ret = Trace::min(ret, hit_recursion(ray,left_idx,ret));
// 		// ret = Trace::min(ret, hit_recursion(ray,right_idx,ret));

// 		return ret;

//     }

    
// }
template<typename Primitive> void BVH<Primitive>::hit_recursion(const Ray& ray, size_t node_idx,Trace& ret) const {
    
	// Vec2 curr_times=Vec2(0.0f,FLT_MAX);
	// bool if_hit = nodes.at(node_idx).bbox.hit(ray,curr_times);
	// // printf("%f,%d\n",nodes.at(node_idx).bbox.surface_area(),if_hit);
	// if(!if_hit){
	// 	return;
	// }
	
	if(nodes[node_idx].is_leaf()) {
        size_t start_idx = nodes[node_idx].start;
        size_t size = nodes[node_idx].size;
        for(size_t i = start_idx; i < start_idx+size; i++) {
			// Vec2 curr_times = ray.dist_bounds;
            Trace hit = primitives[i].hit(ray);
			ret = Trace::min(ret, hit);
        }
        return;
    } else {
        size_t left_idx = nodes[node_idx].l;
        size_t right_idx = nodes[node_idx].r;
		
        Vec2 left_times = ray.dist_bounds;
        Vec2 right_times = ray.dist_bounds;
        bool left_h = nodes[left_idx].bbox.hit(ray,left_times);
        bool right_h = nodes[right_idx].bbox.hit(ray,right_times);
		
		(void) left_h;
		(void) right_h;
		// printf("%zd",n_primitives());
		// if (!left_h && !right_h) return;
		// else{ 
		size_t hit_first = (left_times.x <= right_times.x)?left_idx:right_idx;
		size_t hit_second = (left_times.x <= right_times.x)?right_idx:left_idx;
		
		float hit_second_t = (left_times.x <= right_times.x)?right_times.x:left_times.x;
		(void)hit_second_t;
		hit_recursion(ray,hit_first,ret);
		// if(!ret.hit || (ret.hit && hit_second_t <= ret.distance)) {
		hit_recursion(ray,hit_second,ret);
			// }
			// return;
		// }
		// else if (left_h) {
        // 	hit_recursion(ray,left_idx,ret);
		// 	return;
		// }
		// else if (right_h) {
        // 	hit_recursion(ray,right_idx,ret);
		// 	return;
		// }


		// hit_recursion(ray,left_idx,ret);
		// hit_recursion(ray,right_idx,ret);
		// return;

    }
}

template<typename Primitive>
BVH<Primitive>::BVH(std::vector<Primitive>&& prims, size_t max_leaf_size) {
	build(std::move(prims), max_leaf_size);
}

template<typename Primitive> std::vector<Primitive> BVH<Primitive>::destructure() {
	nodes.clear();
	return std::move(primitives);
}

template<typename Primitive>
template<typename P>
typename std::enable_if<std::is_copy_assignable_v<P>, BVH<P>>::type BVH<Primitive>::copy() const {
	BVH<Primitive> ret;
	ret.nodes = nodes;
	ret.primitives = primitives;
	ret.root_idx = root_idx;
	return ret;
}

template<typename Primitive> Vec3 BVH<Primitive>::sample(RNG &rng, Vec3 from) const {
	if (primitives.empty()) return {};
	int32_t n = rng.integer(0, static_cast<int32_t>(primitives.size()));
	return primitives[n].sample(rng, from);
}

template<typename Primitive>
float BVH<Primitive>::pdf(Ray ray, const Mat4& T, const Mat4& iT) const {
	if (primitives.empty()) return 0.0f;
	float ret = 0.0f;
	for (auto& prim : primitives) ret += prim.pdf(ray, T, iT);
	return ret / primitives.size();
}

template<typename Primitive> void BVH<Primitive>::clear() {
	nodes.clear();
	primitives.clear();
}

template<typename Primitive> bool BVH<Primitive>::Node::is_leaf() const {
	// A node is a leaf if l == r, since all interior nodes must have distinct children
	return l == r;
}

template<typename Primitive>
size_t BVH<Primitive>::new_node(BBox box, size_t start, size_t size, size_t l, size_t r) {
	Node n;
	n.bbox = box;
	n.start = start;
	n.size = size;
	n.l = l;
	n.r = r;
	nodes.push_back(n);
	return nodes.size() - 1;
}
 
template<typename Primitive> BBox BVH<Primitive>::bbox() const {
	if(nodes.empty()) return BBox{Vec3{0.0f}, Vec3{0.0f}};
	return nodes[root_idx].bbox;
}

template<typename Primitive> size_t BVH<Primitive>::n_primitives() const {
	return primitives.size();
}

template<typename Primitive>
uint32_t BVH<Primitive>::visualize(GL::Lines& lines, GL::Lines& active, uint32_t level,
                                   const Mat4& trans) const {

	std::stack<std::pair<size_t, uint32_t>> tstack;
	tstack.push({root_idx, 0u});
	uint32_t max_level = 0u;

	if (nodes.empty()) return max_level;

	while (!tstack.empty()) {

		auto [idx, lvl] = tstack.top();
		max_level = std::max(max_level, lvl);
		const Node& node = nodes[idx];
		tstack.pop();

		Spectrum color = lvl == level ? Spectrum(1.0f, 0.0f, 0.0f) : Spectrum(1.0f);
		GL::Lines& add = lvl == level ? active : lines;

		BBox box = node.bbox;
		box.transform(trans);
		Vec3 min = box.min, max = box.max;

		auto edge = [&](Vec3 a, Vec3 b) { add.add(a, b, color); };

		edge(min, Vec3{max.x, min.y, min.z});
		edge(min, Vec3{min.x, max.y, min.z});
		edge(min, Vec3{min.x, min.y, max.z});
		edge(max, Vec3{min.x, max.y, max.z});
		edge(max, Vec3{max.x, min.y, max.z});
		edge(max, Vec3{max.x, max.y, min.z});
		edge(Vec3{min.x, max.y, min.z}, Vec3{max.x, max.y, min.z});
		edge(Vec3{min.x, max.y, min.z}, Vec3{min.x, max.y, max.z});
		edge(Vec3{min.x, min.y, max.z}, Vec3{max.x, min.y, max.z});
		edge(Vec3{min.x, min.y, max.z}, Vec3{min.x, max.y, max.z});
		edge(Vec3{max.x, min.y, min.z}, Vec3{max.x, max.y, min.z});
		edge(Vec3{max.x, min.y, min.z}, Vec3{max.x, min.y, max.z});

		if (!node.is_leaf()) {
			tstack.push({node.l, lvl + 1});
			tstack.push({node.r, lvl + 1});
		} else {
			for (size_t i = node.start; i < node.start + node.size; i++) {
				uint32_t c = primitives[i].visualize(lines, active, level - lvl, trans);
				max_level = std::max(c + lvl, max_level);
			}
		}
	}
	return max_level;
}

template class BVH<Triangle>;
template class BVH<Instance>;
template class BVH<Aggregate>;
template BVH<Triangle> BVH<Triangle>::copy<Triangle>() const;

} // namespace PT
