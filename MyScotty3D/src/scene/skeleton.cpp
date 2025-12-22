#include <unordered_set>
#include "skeleton.h"
#include "test.h"
#include <iostream>

void Skeleton::Bone::compute_rotation_axes(Vec3 *x_, Vec3 *y_, Vec3 *z_) const {
	assert(x_ && y_ && z_);
	auto &x = *x_;
	auto &y = *y_;
	auto &z = *z_;

	//y axis points in the direction of extent:
	y = extent.unit();
	//if extent is too short to normalize nicely, point along the skeleton's 'y' axis:
	if (!y.valid()) {
		y = Vec3{0.0f, 1.0f, 0.0f};
	}

	//x gets skeleton's 'x' axis projected to be orthogonal to 'y':
	x = Vec3{1.0f, 0.0f, 0.0f};
	x = (x - dot(x,y) * y).unit();
	if (!x.valid()) {
		//if y perfectly aligns with skeleton's 'x' axis, x, gets skeleton's z axis:
		x = Vec3{0.0f, 0.0f, 1.0f};
		x = (x - dot(x,y) * y).unit(); //(this should do nothing)
	}

	//z computed from x,y:
	z = cross(x,y);

	//x,z rotated by roll:
	float cr = std::cos(roll / 180.0f * PI_F);
	float sr = std::sin(roll / 180.0f * PI_F);
	// x = cr * x + sr * -z;
	// z = cross(x,y);
	std::tie(x, z) = std::make_pair(cr * x + sr * -z, cr * z + sr * x);
}

std::vector< Mat4 > Skeleton::bind_pose() const {
	//A4T2a: bone-to-skeleton transformations in the bind pose
	//(the bind pose does not rotate by Bone::pose)

	std::vector< Mat4 > bind;
	bind.reserve(bones.size());

	//NOTE: bones is guaranteed to be ordered such that parents appear before child bones.

	for (auto const &bone : bones) {
		(void)bone; //avoid complaints about unused bone
		//placeholder -- your code should actually compute the correct transform:
		// bind.emplace_back(Mat4::I);
		auto res = Mat4::I;
		// bone.extent;
		// bone.parent;
		Bone bone_it = bone;
		while(bone_it.parent!=-1U){
			// The transform between every bone's local space to its parent
			// For bone with parent 
			auto parent = bones[bone_it.parent];
			res = Mat4::translate(parent.extent) * res;
			bone_it = parent;
		}
		// without parent:
		res = Mat4::translate(base) * res;
		bind.emplace_back(res);
		// Mat4::angle_axis()
		// The functions Mat4::angle_axis and Mat4::translate will be helpful
	}

	assert(bind.size() == bones.size()); //should have a transform for every bone.
	return bind;
}

std::vector< Mat4 > Skeleton::current_pose() const {
    //A4T2a: bone-to-skeleton transformations in the current pose

	//Similar to bind_pose(), but takes rotation from Bone::pose into account.
	// (and translation from Skeleton::base_offset!)

	//You'll probably want to write a loop similar to bind_pose().

	//Useful functions:
	//Bone::compute_rotation_axes() will tell you what axes (in local bone space) Bone::pose should rotate around.
	//Mat4::angle_axis(angle, axis) will produce a matrix that rotates angle (in degrees) around a given axis.

	// return std::vector< Mat4 >(bones.size(), Mat4::I);
	std::vector< Mat4 > cur;
	cur.reserve(bones.size());
	//NOTE: bones is guaranteed to be ordered such that parents appear before child bones.

	for (auto const &bone : bones) {
		(void)bone; //avoid complaints about unused bone
		//placeholder -- your code should actually compute the correct transform:
		// bind.emplace_back(Mat4::I);
		// bone.extent;
		// bone.parent;
		Bone bone_it = bone;
		auto res = Mat4::I;

		// get function inputs
		// aligned with the bone's extent
		// formula
		Vec3 yb;
		Vec3 x_tmp;
		float theta_b;
		Vec3 xb;
		Vec3 zb;

		while(bone_it.parent!=-1U){
			// The transform between every bone's local space to its parent
			// For bone with parent 
			auto parent = bones[bone_it.parent];
			// get rotation 
			yb = bone_it.extent.normalize();
			x_tmp = (Vec3(1,0,0) - (Vec3(1,0,0)*yb)*yb).normalize();
			theta_b = bone_it.roll;
			xb = cos(theta_b) * x_tmp - sin(theta_b)*(cross(x_tmp,yb));
			zb = cross(xb,yb);
			// R mat
			bone_it.compute_rotation_axes(&xb,&yb,&zb);
			auto R = Mat4::angle_axis(bone_it.pose.z, zb) *Mat4::angle_axis(bone_it.pose.y, yb) 
			* Mat4::angle_axis(bone_it.pose.x, xb); 
			// all
			res = Mat4::translate(parent.extent) * R * res;
			bone_it = parent;
		}
		// without parent:
		yb = bone_it.extent.normalize();
		x_tmp = (Vec3(1,0,0) - (Vec3(1,0,0)*yb)*yb).normalize();
		theta_b = bone_it.roll;
		xb = cos(theta_b) * x_tmp - sin(theta_b)*(cross(x_tmp,yb));
		zb = cross(xb,yb);
		bone_it.compute_rotation_axes(&xb,&yb,&zb);

		auto R = Mat4::angle_axis(bone_it.pose.z, zb) *Mat4::angle_axis(bone_it.pose.y, yb) 
		* Mat4::angle_axis(bone_it.pose.x, xb); 
		res = Mat4::translate(base+base_offset) * R * res;

		cur.emplace_back(res);
		// Mat4::angle_axis()
		// The functions Mat4::angle_axis and Mat4::translate will be helpful
	}

	assert(cur.size() == bones.size()); //should have a transform for every bone.
	return cur;

}

std::vector< Vec3 > Skeleton::gradient_in_current_pose() const {
    //A4T2b: IK gradient

    // Computes the gradient (partial derivative) of IK energy relative to each bone's Bone::pose, in the current pose.

	//The IK energy is the sum over all *enabled* handles of the squared distance from the tip of Handle::bone to Handle::target
	std::vector< Vec3 > gradient(bones.size(), Vec3{0.0f, 0.0f, 0.0f});

	//TODO: loop over handles and over bones in the chain leading to the handle, accumulating gradient contributions.
	//remember bone.compute_rotation_axes() -- should be useful here, too!
	for (auto const &handle : handles) {
		(void)handle; //avoid complaints about unused bone
		Vec3 h = handle.target;
		if(!handle.enabled){
			continue;
		}
		auto cur = current_pose();
		// std::cout<<cur[handle.bone];
		Vec3 p = cur[handle.bone]*bones[handle.bone].extent;

		// std::cout<<"here1"<<std::endl;
		// over all *enabled* handles
		for(BoneIndex b = handle.bone;b<bones.size();b=bones[b].parent){
			// if(b==-1U){
			// 	break;
			// }
			// std::cout<<"here"<<std::endl;
			Bone const &bone = bones[b];

			Vec3 extent_tmp = bone.extent;
			Vec3 yb = extent_tmp.normalize();
			auto x_tmp = (Vec3(1,0,0) - (Vec3(1,0,0)*yb)*yb).normalize();
			auto theta_b = bone.roll;
			Vec3 xb = cos(theta_b) * x_tmp - sin(theta_b)*(cross(x_tmp,yb));
			Vec3 zb = cross(xb,yb);
			// R mat
			bone.compute_rotation_axes(&xb,&yb,&zb);
			// only z
			// auto Rz = Mat4::angle_axis(bone.pose.z, zb);

			Mat4 xf;
			auto parent_idx = bone.parent;
			if(parent_idx!=-1U){
				auto parent = bones[parent_idx];
				xf = cur[parent_idx] * Mat4::translate(parent.extent);
			}
			else{
				xf = Mat4::translate(base+base_offset);
			}

			// cur[parent_idx] *  Mat4::translate(parent.extent) * R;// compute linear transform
			Vec3 r = xf * Vec3{0.f,0.f,0.f};


			Vec3 x = (xf * Mat4::angle_axis(bone.pose.z, zb) *Mat4::angle_axis(bone.pose.y, yb)).rotate(xb);// x
			Vec3 y = (xf * Mat4::angle_axis(bone.pose.z, zb)).rotate(yb);// y
			Vec3 z = xf.rotate(zb) ;// z

			gradient[b].x += dot(cross(x,p-r),p-h);
			gradient[b].y += dot(cross(y,p-r),p-h);
			gradient[b].z += dot(cross(z,p-r),p-h);

		}
	}
	// std::cout<<"here"<<std::endl;

	assert(gradient.size() == bones.size());
	return gradient;
}

bool Skeleton::solve_ik(uint32_t steps) {
	//A4T2b - gradient descent
	//check which handles are enabled
	
	for (auto const &handle : handles) {
		//check which handles are enabled
		if(!handle.enabled){
			continue;
		}
		auto cur = current_pose();
		//run `steps` iterations
		for(uint32_t i=0;i<steps;i++){
			//call gradient_in_current_pose() to compute d loss / d pose
			auto grad=gradient_in_current_pose();
			//add ...
			BoneIndex b = handle.bone;
			float lr = 0.5f;
			bones[b].pose.x -= lr * grad[b].x;
			bones[b].pose.y -= lr * grad[b].y;
			bones[b].pose.z -= lr * grad[b].z;
			//if at a local minimum (e.g., gradient is near-zero), return 'true'.

			Vec3 p = cur[handle.bone]*bones[handle.bone].extent;
			Vec3 h = handle.target;
			float loss = (p-h).norm_squared();

			if(loss<0.0001){
				return false;
			}
		}
		//if run through all steps, return `false`.

	}
	
	return false;
}

Vec3 Skeleton::closest_point_on_line_segment(Vec3 const &a, Vec3 const &b, Vec3 const &p) {
	//A4T3: bone weight computation (closest point helper)

    // Return the closest point to 'p' on the line segment from a to b

	//Efficiency note: you can do this without any sqrt's! (no .unit() or .norm() is needed!)
	Vec3 ab = b-a;
	Vec3 pa = p-a;
	float ab_norm = ab.norm();
	// do normalize here!
	float t = dot(pa,ab)*1.f/ab_norm;

	// std::cout<<t;
	// !! after and before end/start!!
	if(t<=0){
		return a;
	}
	if(t>=ab_norm){
		return b;
	}

	// use unin here!
    return a + t*ab.unit();
}

void Skeleton::assign_bone_weights(Halfedge_Mesh *mesh_) const {
	assert(mesh_);
	auto &mesh = *mesh_;
	(void)mesh; //avoid complaints about unused mesh

	//A4T3: bone weight computation

	//visit every vertex and **set new values** in Vertex::bone_weights (don't append to old values)
	
	for(auto & v:mesh.vertices){
		// v.bone_weights;
		Vec3 p = v.position;
		// distance from vertex to the closest point on bone
		// uint32_t idx = 0;
		
		// (don't append to old values)
		auto new_weights = v.bone_weights;
		new_weights.clear();
		float w_sum = 0.f;
		// for(auto &bone:bones){
		auto bind_positions = bind_pose();
		for(uint32_t i =0;i<bones.size();i++){
			// //be sure to use bone positions in the bind pose (not the current pose!)
			auto bone = bones[i];
			
			float r = bone.radius;
			Vec3 b = bind_positions[i] *bone.extent;
			// Vec3 b = bone.extent; //bind_positions[i] 
			Vec3 a =bind_positions[i] *Vec3();
			// Vec3 b = bone.extent;
			// Vec3 a = Vec3();

			Vec3 closed_pos = closest_point_on_line_segment(a,b,p);
			float dis = (closed_pos-p).norm();
			float w = std::max(0.f,r-dis)/r;
			// std::cout<<dis<<","<<r<<", "<<p<<std::endl;
			// **set new values** in Vertex::bone_weights
			// you should only store nonzero bone weights in Vertex::bone_weights
			// Bone_Weight 
			if(w>0.f){
				new_weights.push_back({i,w});
				w_sum += w;
			}
			// idx++;
		}
		if(w_sum>0.f){
			for(uint32_t j=0;j<new_weights.size();j++){
				new_weights[j].weight/=w_sum;
			}
		}
		v.bone_weights = new_weights;
	}

	//you should fill in the helper closest_point_on_line_segment() before working on this function

}

Indexed_Mesh Skeleton::skin(Halfedge_Mesh const &mesh, std::vector< Mat4 > const &bind, std::vector< Mat4 > const &current) {
	assert(bind.size() == current.size());


	//A4T3: linear blend skinning

	//one approach you might take is to first compute the skinned positions (at every vertex) and normals (at every corner)
	// then generate faces in the style of Indexed_Mesh::from_halfedge_mesh

	//---- step 1: figure out skinned positions ---

	std::unordered_map< Halfedge_Mesh::VertexCRef, Vec3 > skinned_positions;
	std::unordered_map< Halfedge_Mesh::HalfedgeCRef, Vec3 > skinned_normals;
	//reserve hash table space to (one hopes) avoid re-hashing:
	skinned_positions.reserve(mesh.vertices.size());
	skinned_normals.reserve(mesh.halfedges.size());

	//(you will probably want to precompute some bind-to-current transformation matrices here)
	auto bind2cur = bind;
    for (size_t i = 0; i < bind.size(); ++i) {
        bind2cur[i] = current[i] * bind[i].inverse();
    }
	for (auto vi = mesh.vertices.begin(); vi != mesh.vertices.end(); ++vi) {
		// skinned_positions.emplace(vi, vi->position); //PLACEHOLDER! Replace with code that computes the position of the vertex according to vi->position and vi->bone_weights.
		//NOTE: vertices with empty bone_weights should remain in place.
		Vec3 v_pos = vi->position;
		Vec3 skinned_position= Vec3();
		auto weights = vi->bone_weights;
		Mat4 sum_wei = Mat4();
		if(weights.size()==0){
			sum_wei = Mat4();
		}else{
			sum_wei = Mat4::Zero;
		}
        for (auto &w : weights) {
			float w_w = w.weight;
			uint32_t bone_idx = w.bone;
			auto cur = bind2cur[bone_idx];
            sum_wei += w_w * cur;
        }
		skinned_position = sum_wei * v_pos;
        skinned_positions.emplace(vi, skinned_position);

		//circulate corners at this vertex:
		auto h = vi->halfedge;
		do {
			//NOTE: could skip if h->face->boundary, since such corners don't get emitted
			if(h->face->boundary){
				continue;
			}
			// skinned_normals.emplace(h, h->corner_normal); //PLACEHOLDER! Replace with 
			// code that properly transforms the normal vector! Make sure that you normalize correctly.
			Mat4 sum_wei2 = Mat4();
			if(weights.size()==0){
				sum_wei2 = Mat4();
			}else{
				sum_wei2 = Mat4::Zero;
			}
			for (auto &w : weights) {
				float w_w = w.weight;
				uint32_t bone_idx = w.bone;
				auto cur = bind2cur[bone_idx];

				sum_wei2 += w_w * cur;
        	}
			// still need normalizing?
			auto skinned_normal = (sum_wei2 * h->corner_normal).normalize();
            skinned_normals.emplace(h, skinned_normal);
			h = h->twin->next;
		} while (h != vi->halfedge);
	}

	//---- step 2: transform into an indexed mesh ---

	//Hint: you should be able to use the code from Indexed_Mesh::from_halfedge_mesh (SplitEdges version) pretty much verbatim, you'll just need to fill in the positions and normals.
	// return Indexed_Mesh::from_halfedge_mesh(mesh, Indexed_Mesh::SplitEdges); //PLACEHOLDER! you'll probably want to copy the SplitEdges case from this function o'er here and modify it to use skinned_positions and skinned_normals.
	
	std::vector<Indexed_Mesh::Vert> verts;
	std::vector<Indexed_Mesh::Index> idxs;
	for (Halfedge_Mesh::FaceCRef f = mesh.faces.begin(); f != mesh.faces.end(); f++) {
			if (f->boundary) continue;
			//every corner gets its own copy of a vertex:
			uint32_t corners_begin = static_cast<uint32_t>(verts.size());
			Halfedge_Mesh::HalfedgeCRef h = f->halfedge;
			do {
				Indexed_Mesh::Vert vert;
				vert.pos = skinned_positions[h->vertex];
				vert.norm = skinned_normals[h];
				vert.uv = h->corner_uv;
				vert.id = f->id;
				verts.emplace_back(vert);
				h = h->next;
			} while (h != f->halfedge);
			uint32_t corners_end = static_cast<uint32_t>(verts.size());

			//divide face into a triangle fan:
			for (size_t i = corners_begin + 1; i + 1 < corners_end; i++) {
				idxs.emplace_back(corners_begin);
				idxs.emplace_back(static_cast<uint32_t>(i));
				idxs.emplace_back(static_cast<uint32_t>(i+1));
			}
	}
	auto result =  Indexed_Mesh(std::move(verts), std::move(idxs));
	return result;
}

void Skeleton::for_bones(const std::function<void(Bone&)>& f) {
	for (auto& bone : bones) {
		f(bone);
	}
}


void Skeleton::erase_bone(BoneIndex bone) {
	assert(bone < bones.size());
	//update indices in bones:
	for (uint32_t b = 0; b < bones.size(); ++b) {
		if (bones[b].parent == -1U) continue;
		if (bones[b].parent == bone) {
			assert(b > bone); //topological sort!
			//keep bone tips in the same place when deleting parent bone:
			bones[b].extent += bones[bone].extent;
			bones[b].parent = bones[bone].parent;
		} else if (bones[b].parent > bone) {
			assert(b > bones[b].parent); //topological sort!
			bones[b].parent -= 1;
		}
	}
	// erase the bone
	bones.erase(bones.begin() + bone);
	//update indices in handles (and erase any handles on this bone):
	for (uint32_t h = 0; h < handles.size(); /* later */) {
		if (handles[h].bone == bone) {
			erase_handle(h);
		} else if (handles[h].bone > bone) {
			handles[h].bone -= 1;
			++h;
		} else {
			++h;
		}
	}
}

void Skeleton::erase_handle(HandleIndex handle) {
	assert(handle < handles.size());

	//nothing internally refers to handles by index so can just delete:
	handles.erase(handles.begin() + handle);
}


Skeleton::BoneIndex Skeleton::add_bone(BoneIndex parent, Vec3 extent) {
	assert(parent == -1U || parent < bones.size());
	Bone bone;
	bone.extent = extent;
	bone.parent = parent;
	//all other parameters left as default.

	//slightly unfortunate hack:
	//(to ensure increasing IDs within an editing session, but reset on load)
	std::unordered_set< uint32_t > used;
	for (auto const &b : bones) {
		used.emplace(b.channel_id);
	}
	while (used.count(next_bone_channel_id)) ++next_bone_channel_id;
	bone.channel_id = next_bone_channel_id++;

	//all other parameters left as default.

	BoneIndex index = BoneIndex(bones.size());
	bones.emplace_back(bone);

	return index;
}

Skeleton::HandleIndex Skeleton::add_handle(BoneIndex bone, Vec3 target) {
	assert(bone < bones.size());
	Handle handle;
	handle.bone = bone;
	handle.target = target;
	//all other parameters left as default.

	//slightly unfortunate hack:
	//(to ensure increasing IDs within an editing session, but reset on load)
	std::unordered_set< uint32_t > used;
	for (auto const &h : handles) {
		used.emplace(h.channel_id);
	}
	while (used.count(next_handle_channel_id)) ++next_handle_channel_id;
	handle.channel_id = next_handle_channel_id++;

	HandleIndex index = HandleIndex(handles.size());
	handles.emplace_back(handle);

	return index;
}


Skeleton Skeleton::copy() {
	//turns out that there aren't any fancy pointer data structures to fix up here.
	return *this;
}

void Skeleton::make_valid() {
	for (uint32_t b = 0; b < bones.size(); ++b) {
		if (!(bones[b].parent == -1U || bones[b].parent < b)) {
			warn("bones[%u].parent is %u, which is not < %u; setting to -1.", b, bones[b].parent, b);
			bones[b].parent = -1U;
		}
	}
	if (bones.empty() && !handles.empty()) {
		warn("Have %u handles but no bones. Deleting handles.", uint32_t(handles.size()));
		handles.clear();
	}
	for (uint32_t h = 0; h < handles.size(); ++h) {
		if (handles[h].bone >= HandleIndex(bones.size())) {
			warn("handles[%u].bone is %u, which is not < bones.size(); setting to 0.", h, handles[h].bone);
			handles[h].bone = 0;
		}
	}
}

//-------------------------------------------------

Indexed_Mesh Skinned_Mesh::bind_mesh() const {
	return Indexed_Mesh::from_halfedge_mesh(mesh, Indexed_Mesh::SplitEdges);
}

Indexed_Mesh Skinned_Mesh::posed_mesh() const {
	return Skeleton::skin(mesh, skeleton.bind_pose(), skeleton.current_pose());
}

Skinned_Mesh Skinned_Mesh::copy() {
	return Skinned_Mesh{mesh.copy(), skeleton.copy()};
}
