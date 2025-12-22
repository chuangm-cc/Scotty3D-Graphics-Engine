
#include "halfedge.h"

#include <unordered_map>
#include <unordered_set>
#include <functional>
#include <iostream>

/******************************************************************
*********************** Local Operations **************************
******************************************************************/

/* Note on local operation return types:

    The local operations all return a std::optional<T> type. This is used so that your
    implementation can signify that it cannot perform an operation (i.e., because
    the resulting mesh does not have a valid representation).

    An optional can have two values: std::nullopt, or a value of the type it is
    parameterized on. In this way, it's similar to a pointer, but has two advantages:
    the value it holds need not be allocated elsewhere, and it provides an API that
    forces the user to check if it is null before using the value.

    In your implementation, if you have successfully performed the operation, you can
    simply return the required reference:

            ... collapse the edge ...
            return collapsed_vertex_ref;

    And if you wish to deny the operation, you can return the null optional:

            return std::nullopt;

    Note that the stubs below all reject their duties by returning the null optional.
*/


/*
 * add_face: add a standalone face to the mesh
 *  sides: number of sides
 *  radius: distance from vertices to origin
 *
 * We provide this method as an example of how to make new halfedge mesh geometry.
 */
std::optional<Halfedge_Mesh::FaceRef> Halfedge_Mesh::add_face(uint32_t sides, float radius) {
	//faces with fewer than three sides are invalid, so abort the operation:
	if (sides < 3) return std::nullopt;


	std::vector< VertexRef > face_vertices;
	//In order to make the first edge point in the +x direction, first vertex should
	// be at -90.0f - 0.5f * 360.0f / float(sides) degrees, so:
	float const start_angle = (-0.25f - 0.5f / float(sides)) * 2.0f * PI_F;
	for (uint32_t s = 0; s < sides; ++s) {
		float angle = float(s) / float(sides) * 2.0f * PI_F + start_angle;
		VertexRef v = emplace_vertex();
		v->position = radius * Vec3(std::cos(angle), std::sin(angle), 0.0f);
		face_vertices.emplace_back(v);
	}

	assert(face_vertices.size() == sides);

	//assemble the rest of the mesh parts:
	FaceRef face = emplace_face(false); //the face to return
	FaceRef boundary = emplace_face(true); //the boundary loop around the face

	std::vector< HalfedgeRef > face_halfedges; //will use later to set ->next pointers

	for (uint32_t s = 0; s < sides; ++s) {
		//will create elements for edge from a->b:
		VertexRef a = face_vertices[s];
		VertexRef b = face_vertices[(s+1)%sides];

		//h is the edge on face:
		HalfedgeRef h = emplace_halfedge();
		//t is the twin, lies on boundary:
		HalfedgeRef t = emplace_halfedge();
		//e is the edge corresponding to h,t:
		EdgeRef e = emplace_edge(false); //false: non-sharp

		//set element data to something reasonable:
		//(most ops will do this with interpolate_data(), but no data to interpolate here)
		h->corner_uv = a->position.xy() / (2.0f * radius) + 0.5f;
		h->corner_normal = Vec3(0.0f, 0.0f, 1.0f);
		t->corner_uv = b->position.xy() / (2.0f * radius) + 0.5f;
		t->corner_normal = Vec3(0.0f, 0.0f,-1.0f);

		//thing -> halfedge pointers:
		e->halfedge = h;
		a->halfedge = h;
		if (s == 0) face->halfedge = h;
		if (s + 1 == sides) boundary->halfedge = t;

		//halfedge -> thing pointers (except 'next' -- will set that later)
		h->twin = t;
		h->vertex = a;
		h->edge = e;
		h->face = face;

		t->twin = h;
		t->vertex = b;
		t->edge = e;
		t->face = boundary;

		face_halfedges.emplace_back(h);
	}

	assert(face_halfedges.size() == sides);

	for (uint32_t s = 0; s < sides; ++s) {
		face_halfedges[s]->next = face_halfedges[(s+1)%sides];
		face_halfedges[(s+1)%sides]->twin->next = face_halfedges[s]->twin;
	}

	return face;
}


/*
 * bisect_edge: split an edge without splitting the adjacent faces
 *  e: edge to split
 *
 * returns: added vertex
 *
 * We provide this as an example for how to implement local operations.
 * (and as a useful subroutine!)
 */
std::optional<Halfedge_Mesh::VertexRef> Halfedge_Mesh::bisect_edge(EdgeRef e) {
	// Phase 0: draw a picture
	//
	// before:
	//    ----h--->
	// v1 ----e--- v2
	//   <----t---
	//
	// after:
	//    --h->    --h2->
	// v1 --e-- vm --e2-- v2
	//    <-t2-    <--t--
	//

	// Phase 1: collect existing elements
	HalfedgeRef h = e->halfedge;
	HalfedgeRef t = h->twin;
	VertexRef v1 = h->vertex;
	VertexRef v2 = t->vertex;

	// Phase 2: Allocate new elements, set data
	VertexRef vm = emplace_vertex();
	vm->position = (v1->position + v2->position) / 2.0f;
	interpolate_data({v1, v2}, vm); //set bone_weights

	EdgeRef e2 = emplace_edge();
	e2->sharp = e->sharp; //copy sharpness flag

	HalfedgeRef h2 = emplace_halfedge();
	interpolate_data({h, h->next}, h2); //set corner_uv, corner_normal

	HalfedgeRef t2 = emplace_halfedge();
	interpolate_data({t, t->next}, t2); //set corner_uv, corner_normal

	// The following elements aren't necessary for the bisect_edge, but they are here to demonstrate phase 4
    FaceRef f_not_used = emplace_face();
    HalfedgeRef h_not_used = emplace_halfedge();

	// Phase 3: Reassign connectivity (careful about ordering so you don't overwrite values you may need later!)

	vm->halfedge = h2;

	e2->halfedge = h2;

	assert(e->halfedge == h); //unchanged

	//n.b. h remains on the same face so even if h->face->halfedge == h, no fixup needed (t, similarly)

	h2->twin = t;
	h2->next = h->next;
	h2->vertex = vm;
	h2->edge = e2;
	h2->face = h->face;

	t2->twin = h;
	t2->next = t->next;
	t2->vertex = vm;
	t2->edge = e;
	t2->face = t->face;
	
	h->twin = t2;
	h->next = h2;
	assert(h->vertex == v1); // unchanged
	assert(h->edge == e); // unchanged
	//h->face unchanged

	t->twin = h2;
	t->next = t2;
	assert(t->vertex == v2); // unchanged
	t->edge = e2;
	//t->face unchanged


	// Phase 4: Delete unused elements
    erase_face(f_not_used);
    erase_halfedge(h_not_used);

	// Phase 5: Return the correct iterator
	return vm;
}


/*
 * split_edge: split an edge and adjacent (non-boundary) faces
 *  e: edge to split
 *
 * returns: added vertex. vertex->halfedge should lie along e
 *
 * Note that when splitting the adjacent faces, the new edge
 * should connect to the vertex ccw from the ccw-most end of e
 * within the face.
 *
 * Do not split adjacent boundary faces.
 */
std::optional<Halfedge_Mesh::VertexRef> Halfedge_Mesh::split_edge(EdgeRef e) {
	// A2L2 (REQUIRED): split_edge
	
	auto face1 = e->halfedge->face;
	auto face2 = e->halfedge->twin->face;
	// std::cout<<e->center();
	if(face1->boundary || face2->boundary) {
		// std::cout<<"h";
		
		// create
		auto vmid = emplace_vertex();
		auto e1 = emplace_edge();
		auto e2 = emplace_edge();
		auto hn1 = emplace_halfedge();
		auto hn2 = emplace_halfedge();
		auto hn3 = emplace_halfedge();
		auto hn4 = emplace_halfedge();
		auto f2 = emplace_face();
		// auto f3 = emplace_face();

		auto h1=e->halfedge;
        if(h1->face->boundary) 
		{
			h1 = h1->twin;
		}
		face2 = h1->twin->face;
		// std::cout<<"b"<<h1->face->boundary<<std::endl;
		// std::cout<<h1->id;
		// std::cout<<h1->next->next->next->next->id;

		// collect
		auto h2=h1->next;
		auto h_last = h1;
		while(h_last->next!=h1){
			h_last=h_last->next;
		}

		auto out_last = h1->twin;
		while(out_last->next!=h1->twin){
			out_last = out_last->next;
		}
		

		

		// auto v1 = h1->vertex;
		auto v2 = h2->vertex;
		auto v3 = h2->next->vertex;
		auto h3 = h2->next;
		auto center_e_before = e->center();
		auto f1 = h1->face;

		f1->halfedge = h1;
		// std::cout<<"1"<<std::endl;
		// // std::cout<<"face"<<h3->id<<std::endl;
		// std::cout<<"face"<<f1->id<<std::endl;
		// std::cout<<"face"<<f2->id<<std::endl;
		// std::cout<<"face"<<hn1->id<<std::endl;
		// std::cout<<"face"<<face2->id<<std::endl;
		// std::cout<<"face"<<h1->next->twin->face->id<<std::endl;
		// std::cout<<"face"<<h1->next->twin->face->halfedge->next->next->id<<std::endl;

		// auto f3 = h1->face;
		// disconnect
		// std::cout<<vmid->id<<std::endl;
		// std::cout<<v1->id<<std::endl;
		// std::cout<<v2->id<<std::endl;
		// std::cout<<v3->id<<std::endl;
		// std::cout<<h2->id<<std::endl;
		// std::cout<<h1->twin->id;
		// std::cout<<h2->twin->id;
		// std::cout<<h1->twin->face->id;
		// std::cout<<hn2->id;
		// connect
		// vertex
		vmid->position = center_e_before;
		vmid->halfedge = hn2;

		// edge
		e1 ->halfedge = hn2;
		e2->halfedge = hn3;

		// face
		f1->halfedge = h1;
		f2->halfedge = hn2;

		// halfedge
		h1->next = hn4;
		h2->next = hn3;

		h2->face = f2;

		face2->halfedge = h1->twin;

		// if(h2->twin->face->boundary)
		// 	h2->twin->next = hn1;
		out_last->next = hn1;

		h1->twin->vertex = vmid;

		// std::cout<<"face"<<h2->twin->id<<std::endl;
		// std::cout<<"face"<<h3->twin->id<<std::endl;

		hn1->set_tnvef(hn2,h1->twin,v2,e1,face2);
		hn2->set_tnvef(hn1,h2,vmid,e1,f2);
		hn3->set_tnvef(hn4,hn2,v3,e2,f2);
		hn4->set_tnvef(hn3,h3,vmid,e2,f1);
		// std::cout<<v3->position.y;

		vmid->halfedge = h1->twin;

		e->halfedge= h1;

		// std::cout<<"face"<<h1->next->twin->face->halfedge->next->next->id<<std::endl;
		
		return vmid;
    }
	// std::cout<<"hh";
	
	// create
	auto vmid = emplace_vertex();

	auto e1 = emplace_edge();
	auto e2 = emplace_edge();
	auto e3 = emplace_edge();
	// auto e4 = emplace_edge();

	// new added halfedge
	auto hn1 = emplace_halfedge();
	auto hn2 = emplace_halfedge();
	auto hn3 = emplace_halfedge();
	auto hn4 = emplace_halfedge();
	auto hn5 = emplace_halfedge();
	auto hn6 = emplace_halfedge();

	//
	auto f3 = emplace_face();
	auto f4 = emplace_face();
	// auto f5 = emplace_face();


	// collect
	auto t1=e->halfedge;
	auto h1=t1->twin;

	auto t2 = t1->next;
	auto h2 = h1->next;

	auto t3 = t2->next;
	auto h3 = h2->next;

	// auto v1 = h1->vertex;
	auto v2 =  h2->vertex;
	auto v3 = t3->vertex;
	auto v4 = h3->vertex;

	// !! consider about face!
	auto f1 = h1->face;
	auto f2 = t1->face;

	f1->halfedge = h1;
	f2->halfedge = t1;


	auto t_last = t1;
	while(t_last->next!=t1){
		t_last=t_last->next;
	}

	auto h_last = h1;
	while(h_last->next!=h1){
		h_last=h_last->next;
	}

	// std::cout<<h3->id<<std::endl;
	// std::cout<<t_last->id<<std::endl;
	// std::cout<<f2->id<<std::endl;
	// std::cout<<h1->id<<std::endl;
	// std::cout<<h2->id<<std::endl;
	// std::cout<<h3->id<<std::endl;
	// std::cout<<t1->id<<std::endl;
	// std::cout<<t2->id<<std::endl;
	// std::cout<<t3->id<<std::endl;

	// std::cout<<hn1->id<<std::endl;
	// std::cout<<hn2->id<<std::endl;
	// std::cout<<hn3->id<<std::endl;
	// std::cout<<hn4->id<<std::endl;
	// std::cout<<hn5->id<<std::endl;
	// std::cout<<hn6->id<<std::endl;
	// std::cout<<t3->id<<std::endl;
	// std::cout<<e->center();
	auto center_e_before = e->center();
	//disconnect 
	v2->halfedge = h2;
	t1->vertex = vmid;

	//connect
	// vertix:
	
	// std::cout<<center_e_before;
	vmid->position = center_e_before;
	vmid->halfedge = t1;

	// edge
	e1->halfedge = hn1;
	e2->halfedge = hn3;
	e3->halfedge = hn5;

	// face
	f3->halfedge = hn2;
	f4->halfedge = hn4;
	f2->halfedge = t1;

	// hafledge
	h1->next = hn3;

	t_last->face = f3;
	t_last->next = hn5;
	// t_last->face = f3;

	t2->next = hn1;
	t3->face = f3;

	h2->next = hn4;
	h2->face = f4;

	// twin next vertex edge face
	hn1->set_tnvef(hn2,t1,v3,e1,f2);
	hn2->set_tnvef(hn1,t3,vmid,e1,f3);
	hn3->set_tnvef(hn4,h3,vmid,e2,f1);
	hn4->set_tnvef(hn3,hn6,v4,e2,f4);
	hn5->set_tnvef(hn6,hn2,v2,e3,f3);
	hn6->set_tnvef(hn5,h2,vmid,e3,f4);

	vmid->halfedge = t1;

	// std::cout<<e->center();

	// std::cout<<vmid->position;
	

	return vmid;


	// (void)e; //this line avoids 'unused parameter' warnings. You can delete it as you fill in the function.
    // return std::nullopt;


}



/*
 * inset_vertex: divide a face into triangles by placing a vertex at f->center()
 *  f: the face to add the vertex to
 *
 * returns:
 *  std::nullopt if insetting a vertex would make mesh invalid
 *  the inset vertex otherwise
 */
std::optional<Halfedge_Mesh::VertexRef> Halfedge_Mesh::inset_vertex(FaceRef f) {
	// A2Lx4 (OPTIONAL): inset vertex
	
	(void)f;
	if(f->boundary)
    	return std::nullopt;
	
	auto h_it = f->halfedge;
	auto mid = emplace_vertex();
	mid->position = f->center();


	std::vector< HalfedgeRef > h1s;
	std::vector< HalfedgeRef > h2s;
	std::vector< FaceRef > fs;
	std::vector< HalfedgeRef > h_its;
	do{
		auto f1 = emplace_face();
		auto h1 = emplace_halfedge();
		auto h2= emplace_halfedge();
		auto e1 = emplace_edge();

		h_it -> face = f1;
		f1-> halfedge = h1;
		e1->halfedge = h1;

		// twin next vertex edge face
		h1->set_tnvef(h2,h_it,mid,e1,f1);

		// next,face not known
		h2->set_tnvef(h1,h_it,h_it->vertex,e1,f1);

		h1s.push_back(h1);
		h2s.push_back(h2);
		fs.push_back(f1);

		h_its.push_back(h_it);


		// next
		h_it=h_it->next;
	}while(h_it!=f->halfedge);

	for(uint32_t i=1;i<h1s.size();i++){
		h2s[i]->next = h1s[i-1];
		h2s[i]->face = fs[i-1];
	}
	h2s[0]->next = h1s[h1s.size()-1];
	h2s[0]->face = fs[h1s.size()-1];

	for(uint32_t i=0;i<h_its.size()-1;i++){
		h_its[i]->next = h2s[i+1];
	}
	h_its[h_its.size()-1]->next = h2s[0];

	mid->halfedge = h1s[0];

	erase_face(f);
	
	return mid;
}


/* [BEVEL NOTE] Note on the beveling process:

	Each of the bevel_vertex, bevel_edge, and extrude_face functions do not represent
	a full bevel/extrude operation. Instead, they should update the _connectivity_ of
	the mesh, _not_ the positions of newly created vertices. In fact, you should set
	the positions of new vertices to be exactly the same as wherever they "started from."

	When you click on a mesh element while in bevel mode, one of those three functions
	is called. But, because you may then adjust the distance/offset of the newly
	beveled face, we need another method of updating the positions of the new vertices.

	This is where bevel_positions and extrude_positions come in: these functions are
	called repeatedly as you move your mouse, the position of which determines the
	amount / shrink parameters. These functions are also passed an array of the original
	vertex positions, stored just after the bevel/extrude call, in order starting at
	face->halfedge->vertex, and the original element normal, computed just *before* the
	bevel/extrude call.

	Finally, note that the amount, extrude, and/or shrink parameters are not relative
	values -- you should compute a particular new position from them, not a delta to
	apply.
*/

/*
 * bevel_vertex: creates a face in place of a vertex
 *  v: the vertex to bevel
 *
 * returns: reference to the new face
 *
 * see also [BEVEL NOTE] above.
 */
std::optional<Halfedge_Mesh::FaceRef> Halfedge_Mesh::bevel_vertex(VertexRef v) {
	//A2Lx5 (OPTIONAL): Bevel Vertex
	// Reminder: This function does not update the vertex positions.
	// Remember to also fill in bevel_vertex_helper (A2Lx5h)

	(void)v;
    return std::nullopt;
}

/*
 * bevel_edge: creates a face in place of an edge
 *  e: the edge to bevel
 *
 * returns: reference to the new face
 *
 * see also [BEVEL NOTE] above.
 */
std::optional<Halfedge_Mesh::FaceRef> Halfedge_Mesh::bevel_edge(EdgeRef e) {
	//A2Lx6 (OPTIONAL): Bevel Edge
	// Reminder: This function does not update the vertex positions.
	// remember to also fill in bevel_edge_helper (A2Lx6h)

	(void)e;
    return std::nullopt;
}

/*
 * extrude_face: creates a face inset into a face
 *  f: the face to inset
 *
 * returns: reference to the inner face
 *
 * see also [BEVEL NOTE] above.
 */
std::optional<Halfedge_Mesh::FaceRef> Halfedge_Mesh::extrude_face(FaceRef f) {
	//A2L4: Extrude Face
	// Reminder:This function does not update the vertex positions.
	// Remember to also fill in extrude_helper (A2L4h)

	// corner case
	if(f->boundary) {
        return std::nullopt;
    }

	// number of angle/edge of face
    auto angle_num = 0;
	auto h_start = f->halfedge;

	// vector for not created ones, store first, then give value later
	std::vector<Halfedge_Mesh::VertexRef> vs;
    std::vector<Halfedge_Mesh::HalfedgeRef> h1s;
	std::vector<Halfedge_Mesh::HalfedgeRef> h2s;
	std::vector<Halfedge_Mesh::HalfedgeRef> h3s;
	std::vector<Halfedge_Mesh::HalfedgeRef> h4s;
	std::vector<Halfedge_Mesh::FaceRef> fs;

	std::vector<Halfedge_Mesh::HalfedgeRef> h_nexts;

	do {
		// for check
		angle_num+=1;
		// std::cout<<angle_num<<std::endl;
		// collect
	
		// connection edge
		auto h1 = emplace_halfedge();
		auto h2 = emplace_halfedge();
		auto e1 = emplace_edge();
		
		// new face edge
		auto h3 = emplace_halfedge();
		auto h4 = emplace_halfedge();
		auto e2 = emplace_edge();

		// new face and vertice
		auto f1 = emplace_face();
		auto v1 = emplace_vertex();


		// collet for later operations
		vs.push_back(v1);
		h1s.push_back(h1);
		h2s.push_back(h2);
		h3s.push_back(h3);
		h4s.push_back(h4);
		fs.push_back(f1);

		h_nexts.push_back(h_start);

		// connect
		// vertix
		v1->position = h_start->vertex->position;
		v1->halfedge = h4;
		//face
		f1->halfedge = h_start;
		//edge
		e1->halfedge = h1;
		e2->halfedge = h3;

		// for halfedge
		// HalfedgeRef twin_, HalfedgeRef next_, VertexRef vertex_, EdgeRef edge_, FaceRef face_
		h1->set_tnvef(h2,h_start,v1,e1,f1);
		// next and face not created yet,assign later, for now use h_start!
		h2->set_tnvef(h1,h_start,h_start->vertex,e1,f1);
		// vertix not created yet!
		h3->set_tnvef(h4,h1,h_start->vertex,e2,f1);
		// next not created yet!
		h4->set_tnvef(h3,h1,v1,e2,f);
		// std::cout<<f1->id<<std::endl;
		// std::cout<<h2->id<<std::endl;

		// halfedge
		// next not created yet!
		h_start->face = f1;

		h_start = h_start->next;
	} while (h_start!= f->halfedge);

	
	auto total_num = angle_num;
	// std::cout<<total_num<<std::endl;
	// assigne that not deifined in last step

	// h_start = f->halfedge;
	angle_num = 0;
	do {
		// std::cout<<angle_num<<std::endl;
		// next not created ,so assigne next
		if(angle_num==0){
			h2s[angle_num]->next = h3s[total_num-1];
			h2s[angle_num]->face = fs[total_num-1];
		}
		else{
			h2s[angle_num]->next = h3s[angle_num-1];
			h2s[angle_num]->face = fs[angle_num-1];
		}
		// vertix not created yet!
		if(angle_num==total_num-1)
			h3s[angle_num]->vertex = vs[0];
		else
			h3s[angle_num]->vertex = vs[angle_num+1];
		// next not created yet!
		if(angle_num==total_num-1)
			h4s[angle_num]->next = h4s[0];
		else
			h4s[angle_num]->next = h4s[angle_num+1];

		// halfedge
		// next not created yet!
		if(angle_num==total_num-1)
			h_nexts[angle_num]->next = h2s[0];
		else
			h_nexts[angle_num]->next = h2s[angle_num+1];

		// h_start = h_start->next;
		angle_num+=1;
		
	} while (angle_num<total_num);

	// last 
	f->halfedge = h4s[0];
	// std::cout<<f->id;

	// std::cout<<f->id;



    return f;
}

/*
 * flip_edge: rotate non-boundary edge ccw inside its containing faces
 *  e: edge to flip
 *
 * if e is a boundary edge, does nothing and returns std::nullopt
 * if flipping e would create an invalid mesh, does nothing and returns std::nullopt
 *
 * otherwise returns the edge, post-rotation
 *
 * does not create or destroy mesh elements.
 */
std::optional<Halfedge_Mesh::EdgeRef> Halfedge_Mesh::flip_edge(EdgeRef e) {

	//A2L1: Flip Edge
	// edge case1:boundary
	auto face1 = e->halfedge->face;
	if(face1->boundary) {
        return std::nullopt;
    }
	auto face2 = e->halfedge->twin->face;
	if(face2->boundary) {
        return std::nullopt;
    }

	// edge case 4:
	// not enough for flip
	// std::cout<<e->halfedge->vertex->position.y;
	// std::cout<<e->halfedge->twin->vertex->position.y;
	// if(e->halfedge->face->degree()<=3 && e->halfedge->twin->face->degree()<=3){
	// 	return std::nullopt;
	// }

	
	auto v_check1 = e->halfedge->vertex;
	uint32_t d1 = 0;
	uint32_t d11 = 0;
	HalfedgeCRef  h=v_check1->halfedge;
	do {
		if(!h->edge->on_boundary())d1++;
		h = h->twin->next;
		d11++;
	} while (h != v_check1->halfedge);
	// std::cout<<d1;
	auto v_check2 = e->halfedge->twin->vertex;
	uint32_t d2 = 0;
	uint32_t d22 = 0;
	h=v_check2->halfedge;
	do {
		if(!h->edge->on_boundary())d2++;
		h = h->twin->next;
		d22++;
	} while (h != v_check2->halfedge);
	// std::cout<<d22;
	// std::cout<<d11;

	// edge case2: edge is the boundary of shape inside, in another word, the boundry of single shape
	if(d11<=2 || d22<=2){
		return std::nullopt;
	}

	// edge case3: edge between the points making up two triangles. can't be more
	// if(d1>=3 || d2>=3){
	// 	std::cout<<"here";
	// 	return std::nullopt;
	// }



	// get every vertex and haflege
	// collect
	auto h1=e->halfedge;
	auto t1=h1->twin;

	auto t2 = t1->next;
	auto h2 = h1->next;

	auto t3 = t2->next;
	auto h3 = h2->next;

	auto v1 = h1->vertex;
	auto v2 =  h2->vertex;
	auto v3 = t3->vertex;
	auto v4 = h3->vertex;

	// !! consider about face!
	auto f1 = h1->face;
	auto f2 = t1->face;

	f1->halfedge = h1;
	f2->halfedge = t1;

	auto t_last = t1;
	while(t_last->next!=t1){
		t_last=t_last->next;
	}

	auto h_last = h1;
	while(h_last->next!=h1){
		h_last=h_last->next;
	}

	// check
	// std::cout<<face2->id<<std::endl;
	// std::cout<<h1->id<<std::endl;
	// std::cout<<h2->id<<std::endl;
	// std::cout<<h3->id<<std::endl;
	// std::cout<<t1->id<<std::endl;
	// std::cout<<t2->id<<std::endl;
	// std::cout<<t3->id<<std::endl;
	// std::cout<<face1->id<<std::endl;

	// first deal with vertex -> disconnect
	v1->halfedge = t2;
	v2 ->halfedge = h2;

	// for other vertex
	v4->halfedge = t1;
	v3->halfedge = h1;

	// then for halfedge -> consider: next,vertext
	h1->vertex = v3;
	h1->next = h3;

	t1->next = t3;
	t1->vertex = v4;

	// other halfedge:
	h2->next = t1;
	// sepcial case:not consider
	// h3->next = t2;

	t2->next = h1;

	//! last haldedge
	h_last->next = t2;
	t_last->next = h2;

	// face -> for changed
	t2->face = f1;
	h2 ->face = f2;

	// std::cout<<"flip";


    return e;
}


/*
 * make_boundary: add non-boundary face to boundary
 *  face: the face to make part of the boundary
 *
 * if face ends up adjacent to other boundary faces, merge them into face
 *
 * if resulting mesh would be invalid, does nothing and returns std::nullopt
 * otherwise returns face
 */
std::optional<Halfedge_Mesh::FaceRef> Halfedge_Mesh::make_boundary(FaceRef face) {
	//A2Lx7: (OPTIONAL) make_boundary
	auto h0 = face->halfedge;
	auto it = h0;
	
	bool pre_edge_is_boundary = false;

	// for initial
	auto f_he = face->halfedge;
	auto f_he2 = face->halfedge;
	auto f2 = face;
	std::vector<Halfedge_Mesh::EdgeRef> delete_e;
	std::vector<Halfedge_Mesh::HalfedgeRef> delete_h;
	std::vector<Halfedge_Mesh::VertexRef> delete_v;

	bool if_new_face = false;
	do{
		auto it_next = it->next;
		if(!it->edge->on_boundary()){
			auto h1 = it->twin;
			auto h1_last = h1;
			auto h2 = h1->next;

			auto v0 = it->vertex;
			auto v1 = it->twin->vertex;

			do{
				h1_last=h1_last->next;
				// std::cout<<"a";
			}while(h1_last->next!=h1);

			it->next = h1_last->twin;
			h2->twin->next = it;

			pre_edge_is_boundary = false;

			v0->halfedge = it;
			v1->halfedge = h1;

			if(if_new_face){
				f_he = it;
				if_new_face = false;
				// std::cout<<"c";
			}
			else {
				if_new_face = true;
				f_he2 = it;
				// std::cout<<"d";
				f2 = h2->twin->face;
				f2->halfedge = h2->twin;
				it->face = f2;
				// std::cout<<f2->boundary;
				// std::cout<<f2->degree();
			}

		}else{
			// std::cout<<"e";
			if(pre_edge_is_boundary){
				delete_v.push_back(it->vertex);
			}
			delete_e.push_back(it->edge);
			delete_h.push_back(it);
			delete_h.push_back(it->twin);
			pre_edge_is_boundary = true;
		}
		it = it_next;
	}while(it!=h0);

	if(!if_new_face){
		face->boundary = true;
		face->halfedge = f_he;
		// std::cout<<face->id<<std::endl;
		f_he->face = face;
		auto it3 = f_he->next;
		do{
			it3->face = face;
			it3=it3->next;
		}while(it3!=f_he);

	// corner case
	}else{
		face->boundary = true;
		// face = f2;
		auto h4 = f2->halfedge;
		do{
			h4->face = face;
			h4=h4->next;
		}while(h4!=f2->halfedge);
		face->halfedge = h4;
		// face->halfedge = f_he2;
		// std::cout<<face->id;
		erase_face(f2);
	}
	// f_he->next-= face;
	// f_he->next->next->face = f>face ace;
	// std::cout<<face->id<<std::endl;
	// erase_face(face);
	for(auto h:delete_h){
		erase_halfedge(h);
	}
	for(auto e:delete_e){
		// std::cout<<e->halfedge->corner_uv.x<<std::endl;
		// erase_halfedge(e->halfedge);
		// e->halfedge->next = e->halfedge->twin;
		// e->halfedge->twin->next = e->halfedge->next;
		// e->halfedge->corner_uv=Vec2(0.0f,0.0f);
		// erase_halfedge(e->halfedge->twin);
		// e->halfedge->twin->corner_uv=Vec2(0.0f,0.0f);
		// std::cout<<e->halfedge->id<<std::endl;
		// std::cout<<e->halfedge->twin->id<<std::endl;
		erase_edge(e);
		// std::cout<<"b";
	}

	for(auto v:delete_v){
		erase_vertex(v);
	}

	// std::cout<<faces.size();
	// for(auto f:faces){
	// 	std::cout<<f.degree();
	// }
	// std::cout<<std::endl;
	// for(auto v:vertices){
	// 	std::cout<<v.degree();
	// }
	return face;
	// return std::nullopt; //TODO: actually write this code!
}

/*
 * dissolve_vertex: merge non-boundary faces adjacent to vertex, removing vertex
 *  v: vertex to merge around
 *
 * if merging would result in an invalid mesh, does nothing and returns std::nullopt
 * otherwise returns the merged face
 */
std::optional<Halfedge_Mesh::FaceRef> Halfedge_Mesh::dissolve_vertex(VertexRef v) {
	// A2Lx1 (OPTIONAL): Dissolve Vertex
	if(v->on_boundary())
    	return std::nullopt;
	auto count = 0;
	auto h_it = v->halfedge;
	// store all half-edges on the edge
	std::vector<Halfedge_Mesh::HalfedgeRef> edge_hs;
	std::vector<Halfedge_Mesh::HalfedgeRef> delete_hs;
	std::vector<Halfedge_Mesh::EdgeRef> delete_edge;
	std::vector<Halfedge_Mesh::FaceRef> delete_face;

	// for vertix
	std::vector<Halfedge_Mesh::HalfedgeRef> first_hs;
	std::vector<Halfedge_Mesh::VertexRef> vs;
	do{
		auto h_it2 = h_it;
		auto h_start = h_it;
		// store all half edges
		// std::cout<<h_it->next->id<<std::endl;
		auto count2 = 0;
		do{
			h_it2 = h_it2->next;
			edge_hs.push_back(h_it2);
			// std::cout<<count2<<std::endl;
			count2++;
		}while(h_it2->next->next!=h_start);

		first_hs.push_back(h_start->next);
		vs.push_back(h_start->next->vertex);

		delete_hs.push_back(h_it);
		delete_hs.push_back(h_it->twin);
		delete_edge.push_back(h_it->edge);
		delete_face.push_back(h_it->face);
		
		// next
		count++;
		h_it=h_it->twin->next;
		// std::cout<<count<<std::endl;
	}while(h_it!=v->halfedge);


	// get new face and do next
	auto f_new = emplace_face();
	for(uint32_t i = 0;i<edge_hs.size()-1;i++){
		edge_hs[i]->next = edge_hs[i+1];
		edge_hs[i]->face = f_new;
		// std::cout<<edge_hs[i]->id<<std::endl;
	}
	// std::cout<<edge_hs[edge_hs.size()-1]->id<<std::endl;
	edge_hs[edge_hs.size()-1]->next = edge_hs[0];
	edge_hs[edge_hs.size()-1]->face = f_new;

	// for vertices assigne
	for(uint32_t i = 0;i<vs.size();i++){
		// std::cout<<vs[i]->id<<std::endl;
		vs[i]->halfedge = first_hs[i];
		// std::cout<<vs[i]->id<<std::endl;
		// std::cout<<first_hs[i]->id<<std::endl;
	}

	// std::cout<<vs[2]->id<<std::endl;
	// std::cout<<vs[2]->halfedge->id<<std::endl;
	// std::cout<<vs[2]->halfedge->next->id<<std::endl;
	// std::cout<<vs[2]->halfedge->twin->next->twin->id<<std::endl;

	// std::cout<<vs[2]->halfedge->twin->next->id<<std::endl;
	// std::cout<<vs[1]->halfedge->twin->id<<std::endl;
	// std::cout<<vs[1]->id

	// for twin of edge
	// also consider this
	for(uint32_t i = 1;i<vs.size();i++){
		first_hs[i]->twin->next = first_hs[i-1]->twin;

	}
	first_hs[0]->twin->next = first_hs[vs.size()-1]->twin;

	for(uint32_t i = 0;i<vs.size()-1;i++){
		first_hs[i]->twin->vertex = vs[i+1];

	}
	first_hs[vs.size()-1]->twin->vertex = vs[0];

	f_new->halfedge = edge_hs[0];

	// delete
	for(auto h:delete_hs){
		erase_halfedge(h);
	}
	for(auto f:delete_face){
		erase_face(f);
	}
	for(auto e:delete_edge){
		erase_edge(e);
	}

	// std::cout<<v->id;
	erase_vertex(v);
	return f_new;


	
}

/*
 * dissolve_edge: merge the two faces on either side of an edge
 *  e: the edge to dissolve
 *
 * merging a boundary and non-boundary face produces a boundary face.
 *
 * if the result of the merge would be an invalid mesh, does nothing and returns std::nullopt
 * otherwise returns the merged face.
 */
std::optional<Halfedge_Mesh::FaceRef> Halfedge_Mesh::dissolve_edge(EdgeRef e) {
	// A2Lx2 (OPTIONAL): dissolve_edge
	if(e->on_boundary()){
		 return std::nullopt;
	}
	//Reminder: use interpolate_data() to merge corner_uv / corner_normal data
	auto h1 = e->halfedge;
	auto t1 = h1->twin;

	auto h2 = h1->next;
	auto t2=t1->next;

	auto h_last = h1;
	do{h_last=h_last->next;}while(h_last->next!=h1);

	auto t_last = t1;
	do{t_last=t_last->next;}while(t_last->next!=t1);

	auto f1 = h1->face;
	auto f2 = t1->face;

	auto v1 = h1->vertex;
	auto v2 = t1->vertex;

	// disconnect
	auto h_it2 = h2;
	do{h_it2->face = f2;
	h_it2=h_it2->next;
	}while(h_it2!=h1);

	//connect
	h_last->next = t2;
	t_last -> next = h2;

	v1->halfedge = t2;
	v2->halfedge = h2;

	f2->halfedge = t2;

	// erase

	erase_edge(e);
	erase_halfedge(t1);
	erase_halfedge(h1);
	
	erase_face(f1);

	// std::cout<<h2->id<<std::endl;
	// std::cout<<t2->id<<std::endl;
	// std::cout<<h_last->id<<std::endl;
	// std::cout<<t_last->id<<std::endl;

	return f2;
   
}

/* collapse_edge: collapse edge to a vertex at its middle
 *  e: the edge to collapse
 *
 * if collapsing the edge would result in an invalid mesh, does nothing and returns std::nullopt
 * otherwise returns the newly collapsed vertex
 */
std::optional<Halfedge_Mesh::VertexRef> Halfedge_Mesh::collapse_edge(EdgeRef e) {
	//A2L3: Collapse Edge

	//Reminder: use interpolate_data() to merge corner_uv / corner_normal data on halfedges
	// (also works for bone_weights data on vertices!)
	auto face1 = e->halfedge->face;
	auto face2 = e->halfedge->twin->face;

	if(face1->boundary || face2->boundary) {
		// std::cout<<"a";
		auto h1 = e->halfedge;
		auto h2 = h1->next;
		// std::cout<<face2->boundary;
		auto h_last = h1;
		while(h_last->next!=h1) h_last=h_last->next;
		auto midv = emplace_vertex();
		// interpolate_data()
		midv->position = e->center();
		auto v1 = h2->vertex;
		auto v2 = h1->vertex;
		midv->position = (v1->position + v2->position) / 2.0f;
		interpolate_data({v1, v2}, midv); //set bone_weights
		auto f1 = h1->face;
		auto f2 = h1->twin->face;

		bool delete_edge = false;

		auto h_tmp = h2->twin->next;

		auto t1 = h1->twin;
		auto t2 = t1->next;
		auto t_last = t1;


		while(t_last->next!=t1) t_last=t_last->next;

		if(h1==h1->next->next->next){
			delete_edge =true;
			// std::cout<<1<<std::endl;
		}

		// std::cout<<h1->twin->id<<std::endl;
		// std::cout<<h_tmp2->next->id<<std::endl;
		// std::cout<<h_last->face->id<<std::endl;
		// std::cout<<h_tmp3->id<<std::endl;
		// std::cout<<h1->twin->next->id<<std::endl;
		// disconnect
		auto h_it1 = h2;
		do{
			// std::cout<<h_it1->id<<std::endl;
			h_it1->vertex = midv;
			h_it1 = h_it1->twin->next; 
			
		}while(h_it1!=h2);

		h_it1 = h1;
		do{
			// std::cout<<h_it1->id<<std::endl;
			h_it1->vertex = midv;
			h_it1 = h_it1->twin->next; 
			
		}while(h_it1!=h1);

		erase_halfedge(h1->twin);
		erase_halfedge(h1);
		erase_edge(e);
		// h_tmp2->next = h_tmp3;
		// h1->twin->face->halfedge = h_tmp2;
		
		erase_vertex(v1);
		erase_vertex(v2);

		// h_last->next = h2;
		t_last->next = t2;

		f2->halfedge = t2;
		

		// // connect
		if(delete_edge){
			erase_edge(h2->edge);
			erase_halfedge(h2->twin);
			erase_halfedge(h2);
			
			erase_face(f1);

			h_last->next = h_tmp;

			h_last -> face = h_tmp->face;

			h_tmp->next->next->next = h_last;

			midv->halfedge = h_tmp;

			// std::cout<<h_last->id<<std::endl;
			// std::cout<< h_tmp->face->id<<std::endl;
		}
		// std::cout<<h_tmp->face->id<<std::endl;
		else{
			f1->halfedge = h2;

			h_last->next = h2;

			midv->halfedge = h_tmp;
		}
		// // midv->halfedge = h2;
		// f1->halfedge = h2;

	
		// std::cout<<h1->twin->id<<std::endl;
		// std::cout<<h1->id<<std::endl;
	

		return midv;
	}



	// collect
	auto h1 = e->halfedge;
	auto t1 = e->halfedge->twin;
	auto h2 = h1->next;
	auto t2 = t1->next;

	auto h_last = h1;
	auto t_last = t1;
	while(h_last->next!=h1) h_last=h_last->next;
	while(t_last->next!=t1) t_last=t_last->next;

	auto midv = emplace_vertex();
	// interpolate_data()
	midv->position = e->center();
	auto v1 = h2->vertex;
	auto v2 = t2->vertex;
	midv->position = (v1->position + v2->position) / 2.0f;
	interpolate_data({v1, v2}, midv); //set bone_weights

	auto f1 = h1->face;
	auto f2 = t1->face;

	// edge cases
	auto f3 = h2->twin->face;
	auto h_tmp_f3 = h2->twin;

	auto f4 = t2->twin->face;
	auto t_tmp_f4 = t2->twin;

	// find last and second last for edge case
	auto hh2 = h2->twin;
	do{h_tmp_f3=h_tmp_f3->next;hh2=h_tmp_f3->next;}
	while(h_tmp_f3->next->next!=h2->twin);

	hh2 = h2->twin->next->next;
	h_tmp_f3 = h2->twin->next;


	auto tt2 = t2->twin;
	do{t_tmp_f4=t_tmp_f4->next;tt2=t_tmp_f4->next;}
	while(t_tmp_f4->next->next!=t2->twin);




	// h2->twin->next->next;
	// auto tt22 = h2->twin->next;
	// edge case:need to delete dege
	auto delet_edge_h = false;
	auto delet_edge_t = false;
	if(h1->next->next->next == h1){
		delet_edge_h = true;
		// std::cout<<"a";
		// std::cout<<h2->twin->id;
		// std::cout<<h2->twin->next->id;
		// std::cout<<h2->twin->next->next->id;
	}

	if(t1->next->next->next == t1){
		delet_edge_t = true;
		// std::cout<<"b";
	}
	auto h_tmp = h2->twin->next;
	auto t_tmp = t2->twin->next;


	// std::cout<<f1->id<<std::endl;
	// std::cout<<h1->id<<std::endl;
	// std::cout<<v1->id<<std::endl;
	// std::cout<<v2->id<<std::endl;
	// std::cout<<h1->id<<std::endl;

	// disconnect
	// disconnect with all hafledge hwich is from v1 and v2
	auto h_it1 = h2;
	do{
		// std::cout<<h_it1->id<<std::endl;
		h_it1->vertex = midv;
		h_it1 = h_it1->twin->next; 
		
	}while(h_it1!=h2);

	auto t_it1 = t2;
	do{
		// std::cout<<t_it1->id<<std::endl;
		t_it1->vertex = midv;
		t_it1 = t_it1->twin->next; 
	}while(t_it1!=t2);


	// h2->vertex = midv;
	// t2->vertex = midv;


	erase_halfedge(h1);
	erase_edge(e);
	erase_halfedge(t1);
	erase_vertex(v1);
	erase_vertex(v2);

	// connect
	midv->halfedge = h2;
	// std::cout<<midv->id;

	// edge cases
	if(!delet_edge_t){
		
		f2->halfedge = t2;
		t_last->next = t2;
	}else{
		
		tt2 ->next = h_last;
		erase_edge(t2->edge);
		erase_halfedge(t2->twin);
		erase_halfedge(t2);
		erase_face(f2);
		t_last->next = t_tmp;
		t_last -> face = t_tmp->face;
		t_tmp->next->next->next = t_last;
		midv->halfedge = t_tmp;

		t_last->next = t_tmp_f4;
		f4->halfedge = t_tmp_f4;

		
	}

	if(!delet_edge_h){
		f1->halfedge = h2;
		h_last->next = h2;
	}else{
		// h_last ->next = h2;

		hh2 ->next = h_last;
		// h_tmp_f3->next = tt2;
		erase_edge(h2->edge);
		erase_halfedge(h2->twin);
		erase_halfedge(h2);
		erase_face(f1);
		h_last->next = h_tmp;
		h_last -> face = h_tmp->face;
		h_tmp->next->next->next = h_last;
		midv->halfedge = h_tmp;
		// std::cout<<f3->id;
		h_last->next = h_tmp_f3;
		f3->halfedge = h_tmp_f3;

	}
 

	return midv;
    // return std::nullopt;
}

/*
 * collapse_face: collapse a face to a single vertex at its center
 *  f: the face to collapse
 *
 * if collapsing the face would result in an invalid mesh, does nothing and returns std::nullopt
 * otherwise returns the newly collapsed vertex
 */
std::optional<Halfedge_Mesh::VertexRef> Halfedge_Mesh::collapse_face(FaceRef f) {
	//A2Lx3 (OPTIONAL): Collapse Face

	//Reminder: use interpolate_data() to merge corner_uv / corner_normal data on halfedges
	// (also works for bone_weights data on vertices!)
	if(f->boundary)
    	return std::nullopt;
	auto h_it = f->halfedge;
	std::vector<Halfedge_Mesh::HalfedgeRef> delete_hs;
	std::vector<Halfedge_Mesh::EdgeRef> delete_es;
	std::vector<Halfedge_Mesh::VertexRef> delete_vs;
	// auto counts = 0;

	// collect
	auto v_center = emplace_vertex();
	v_center->position = f->center();

	std::vector<uint32_t> degrees;
	auto h_tmp = f->halfedge;
	do{
		degrees.push_back(h_tmp->vertex->degree());
		h_tmp=h_tmp->next;
	}while(h_tmp!=f->halfedge);

	auto index = 0;
	do{
		// std::cout<<count<<std::endl;
		// collect
		auto t1 = h_it->twin;
		auto t2 = t1->next;
		auto t_last = t1;
		auto f1 = t1->face;
		auto e1 = h_it->edge;
		auto v1 = h_it->vertex;

		// auto t3 = t2->twin->next;

		// std::cout<<t1->id<<std::endl;
		// std::cout<<t2->next->id<<std::endl;
		// std::cout<<t2->id<<std::endl;
		// std::cout<<t2->next->next->id<<std::endl;
		
		do{
			t_last = t_last->next;
		}while(t_last->next!=t1);
		// std::cout<<t_last->id<<std::endl;

		// std::cout<<v1->id<<std::endl;

		auto h_last2 = h_it;
		do{
			h_last2 = h_last2->next;
		}while(h_last2->next!=h_it);


		// disconnect
		auto h_it2 = t2;
		
		// std::cout<<n;
		for(uint32_t i=0;i<degrees[index]-2;i++){
			h_it2->vertex = v_center;
			h_it2 = h_it2->twin->next;
		}
		index++;
		// t2->vertex = v_center;
		// t3->vertex = v_center;

		// connect
		t_last->next = t2;
		f1->halfedge = t2;

		
		// t2->vertex = v_center;

		// t3->vertex= v_center;

		delete_hs.push_back(h_it);
		delete_hs.push_back(t1);
		delete_es.push_back(e1);
		delete_vs.push_back(v1);

		v_center->halfedge = t2;

		// next
		h_it = h_it->next;
	}while(h_it!=f->halfedge);

	for(auto h:delete_hs){
		erase_halfedge(h);
	}
	for(auto e:delete_es){
		erase_edge(e);
	}
	for(auto v:delete_vs){
		erase_vertex(v);
	}

	erase_face(f);

	return v_center;
}

/*
 * weld_edges: glue two boundary edges together to make one non-boundary edge
 *  e, e2: the edges to weld
 *
 * if welding the edges would result in an invalid mesh, does nothing and returns std::nullopt
 * otherwise returns e, updated to represent the newly-welded edge
 */
std::optional<Halfedge_Mesh::EdgeRef> Halfedge_Mesh::weld_edges(EdgeRef e, EdgeRef e2) {
	//A2Lx8: Weld Edges

	//Reminder: use interpolate_data() to merge bone_weights data on vertices!
	// edge case
	auto it1 = e->halfedge;
	auto it2 = e2->halfedge;

	do{
		do{
			if(it1->edge->id == it2->edge->id){
				return std::nullopt;
			}
			it1=it1->next;
		}while(it1!=e->halfedge);
		it2=it2->next;
	}while(it2!=e2->halfedge);

	// collect
	auto hb1 = e->halfedge->face->boundary?e->halfedge:e->halfedge->twin;
	auto hb2 = e2->halfedge->face->boundary?e2->halfedge:e2->halfedge->twin;
	// std::cout<<hb2->face->boundary;
	// std::cout<<hb1->face->boundary;

	auto tb1 = hb1->twin;
	auto tb2 = hb2->twin;

	auto tb11 = tb1->next;
	auto tb22 = tb2->next;

	auto tb1_last = tb1;

	auto v0 = tb22->vertex;
	auto v1 = tb2->vertex;

	auto v2 = tb1->vertex;
	auto v3 = hb1->vertex;

	auto f1 = tb11->twin->face;
	auto f2 = tb22->twin->face;
	do{
		tb1_last  = tb1_last ->next;
	}while(tb1_last->next!=tb1);

	auto tb2_last = tb2;
	do{
		tb2_last  = tb2_last ->next;
	}while(tb2_last->next!=tb2);

	// interpolate
	interpolate_data({v1, v3}, v3);
	interpolate_data({v0, v2}, v2);
	v3->position = (v3->position+v1->position)/2.0f;
	v2->position = (v2->position+v0->position)/2.0f;

	// disconnect and connect
	tb22->twin->face = f1;
	tb2_last->twin->face = f1;

	tb11->twin->next = tb2_last->twin;
	tb22->twin->next = tb1_last->twin;

	tb2->edge = e;

	tb2->twin = tb1;
	tb1->twin = tb2;

	tb2->vertex = hb1->vertex;
	hb2->next->vertex = hb1->vertex;

	tb22->vertex = tb1->vertex;

	// erase

	erase_halfedge(hb1);
	erase_halfedge(hb2);
	erase_edge(e2);
	erase_vertex(v0);
	erase_vertex(v1);
	erase_face(f2);

	return e;


	// return std::nullopt;
}



/*
 * bevel_positions: compute new positions for the vertices of a beveled vertex/edge
 *  face: the face that was created by the bevel operation
 *  start_positions: the starting positions of the vertices
 *     start_positions[i] is the starting position of face->halfedge(->next)^i
 *  direction: direction to bevel in (unit vector)
 *  distance: how far to bevel
 *
 * push each vertex from its starting position along its outgoing edge until it has
 *  moved distance `distance` in direction `direction`. If it runs out of edge to
 *  move along, you may choose to extrapolate, clamp the distance, or do something
 *  else reasonable.
 *
 * only changes vertex positions (no connectivity changes!)
 *
 * This is called repeatedly as the user interacts, just after bevel_vertex or bevel_edge.
 * (So you can assume the local topology is set up however your bevel_* functions do it.)
 *
 * see also [BEVEL NOTE] above.
 */
void Halfedge_Mesh::bevel_positions(FaceRef face, std::vector<Vec3> const &start_positions, Vec3 direction, float distance) {
	//A2Lx5h / A2Lx6h (OPTIONAL): Bevel Positions Helper
	
	// The basic strategy here is to loop over the list of outgoing halfedges,
	// and use the preceding and next vertex position from the original mesh
	// (in the start_positions array) to compute an new vertex position.
	
}

/*
 * extrude_positions: compute new positions for the vertices of an extruded face
 *  face: the face that was created by the extrude operation
 *  move: how much to translate the face
 *  shrink: amount to linearly interpolate vertices in the face toward the face's centroid
 *    shrink of zero leaves the face where it is
 *    positive shrink makes the face smaller (at shrink of 1, face is a point)
 *    negative shrink makes the face larger
 *
 * only changes vertex positions (no connectivity changes!)
 *
 * This is called repeatedly as the user interacts, just after extrude_face.
 * (So you can assume the local topology is set up however your extrude_face function does it.)
 *
 * Using extrude face in the GUI will assume a shrink of 0 to only extrude the selected face
 * Using bevel face in the GUI will allow you to shrink and increase the size of the selected face
 * 
 * see also [BEVEL NOTE] above.
 */
void Halfedge_Mesh::extrude_positions(FaceRef face, Vec3 move, float shrink) {
	//A2L4h: Extrude Positions Helper

	//General strategy:
	// use mesh navigation to get starting positions from the surrounding faces,
	// compute the centroid from these positions + use to shrink,
	// offset by move

	// mesh navigation to get starting positions from the surrounding faces,
	auto start_h = face->halfedge;
	auto h = face->halfedge;
	auto center_position = face->center();
	do{
		auto start_position = start_h->vertex->position;
		// std::cout<<face->center();

		// shrink
		if(shrink<=1){
			// positive shrink makes the face smaller (at shrink of 1, face is a point)

			// negative shrink makes the face larger
			auto offset = start_position - center_position;
			start_h->vertex->position = center_position + offset*(1-shrink);
		} else{
			std::cout<<"shrink value invalid"<<std::endl;
		}
		start_h->vertex->position += move;

		start_h=start_h->next;
	}while(start_h!=h);
	
}

