
#include "../geometry/spline.h"

template<typename T> T Spline<T>::at(float time) const {

	// A4T1b: Evaluate a Catumull-Rom spline

	// Given a time, find the nearest positions & tangent values
	// defined by the control point map.

	// Transform them for use with cubic_unit_spline

	// Be wary of edge cases! What if time is before the first knot,
	// before the second knot, etc...

	// If there are no knots at all in the spline, interpolation should return 
	// the default value for the interpolated type. 
	// This value can be computed by simply calling the constructor for the type: T()
	size_t size = knots.size();
	if(size ==0){
		return T();
	}
	// If there is only one knot in the spline, interpolation should always return the value of that knot
	if(size==1){
		return knots.begin()->second;
	}
	// If the query time is less than or equal to the initial knot, return the initial knot's value.
	// auto init = std::lower_bound(knots.begin(), knots.end(), time); 
	auto init = knots.upper_bound(time);
	if(init==knots.begin()){
		return knots.begin()->second;
	}
	// If the query time is greater than or equal to the final knot, return the final knot's value.
	float largest = knots.rbegin()->first;
	if(time>=largest){
		return knots.rbegin()->second;
	}
	
	// Any query time between the first and last knot
	auto k2_it = knots.upper_bound(time);
	auto k1_it = k2_it;
	--k1_it;

	auto k1 = k1_it->second;
	float t1 = k1_it->first;

	auto k2 = k2_it->second;
	float t2 = k2_it->first;

	auto k0_it = k1_it;
	auto k0 = k0_it->second;
	float t0 = k0_it->first;
	//Suppose we don't have a knot "two to the left". Then we will define a "virtual" knot .
	if(k1_it==knots.begin()){
		k0 = k1-(k2-k1);
		t0 = t1-(t2-t1);
	}else{
		--k0_it;
		k0 = k0_it->second;
		t0 = k0_it->first;
	}

	// !!!! not just use end()
	// use +1 to check end()!!
	auto k3_it = k2_it;
	auto k3 = k3_it->second;
	float t3 = k3_it->first;
	// !!! use next
	if(std::next(k2_it) == knots.end()){
		k3 = k2+(k2-k1);
		t3 = t2+(t2-t1);
	}else{
		++k3_it;
		k3 = k3_it->second;
		t3 = k3_it->first;
	}

	// We then use and as the endpoints of our cubic "piece," and for tangents we use the values
	auto m0 = (k2-k0)/(t2-t0);
	auto m1 = (k3-k1)/(t3-t1);

	//time value is between 0 and 1
	// you will have to divide by the length of the current interval
	float res_tt = (time-t1)*1.f/(t2-t1);

	// how this normalization affects the value 
	// computed by the subroutine, in comparison to the values we want to return???

	return cubic_unit_spline(res_tt, k1, k2, m0*(t2-t1), m1*(t2-t1));
}

template<typename T>
T Spline<T>::cubic_unit_spline(float time, const T& position0, const T& position1,
                               const T& tangent0, const T& tangent1) {

	// A4T1a: Hermite Curve over the unit interval

	// Given time in [0,1] compute the cubic spline coefficients and use them to compute
	// the interpolated value at time 'time' based on the positions & tangents
	// Evaluate the time, its square, and its cube (for readability, you may want to make a local copy).
	float t = time;
	float t_2 = pow(t,2.f);
	float t_3 = pow(t,3.f);
	// Using these values, as well as the position and tangent values, compute the four basis functions 
	// of a cubic polynomial in Hermite form.
	float h00 = 2*t_3 - 3*t_2 + 1;
	float h10 = t_3 - 2*t_2 + t;
	float h01 = -2*t_3 + 3*t_2;
	float h11 = t_3 - t_2;
	// Finally, combine the endpoint and tangent data using the evaluated bases, and return the result.
	T res = h00*position0 + h10*tangent0 + h01*position1 + h11*tangent1;
	// Note that Spline is parameterized on type T, which allows us to create splines over
	// any type that supports the * and + operators.

	return res;
}

template class Spline<float>;
template class Spline<double>;
template class Spline<Vec4>;
template class Spline<Vec3>;
template class Spline<Vec2>;
template class Spline<Mat4>;
template class Spline<Spectrum>;
