# A3 – Path Tracing Renderer

## Overview
Build a physically-based path tracer in Scotty3D, including ray generation, intersections, BVH acceleration, path tracing with materials, and direct/environment lighting.

---

## A3T1 – Camera Rays

**Requirement (summary):**
Generate camera rays that correctly fill the view frustum.

### Result

![Your creation](A3-writeup/student/ray.png)

![Your creation](A3-writeup/student/task1.png)

---

## A3T2 – Intersection Tests

**Requirement (summary):**
Implement ray–sphere and ray–triangle intersections.

### Result

![Your creation](A3-writeup/student/a1t22.png)
![Your creation](A3-writeup/student/a1t23.png)
![Your creation](A3-writeup/student/a1t24.png)

![Your creation](A3-writeup/student/task2.png)

---

## A3T3 – BVH

**Requirement (summary):**
Build and traverse BVH for efficient ray-scene intersection.

### Result
![Your creation](A3-writeup/student/t31.png)

---

## A3T4 – Path Tracing

**Requirement (summary):**
Implement recursive path tracing integrator (Lambertian correctness).

### Result
![Your creation](A3-writeup/student/t4.png)
![Your creation](A3-writeup/student/t42.png)


---

## A3T5 – Materials

**Requirement (summary):**
Support mirror/refraction/glass materials.

### Result
![Your creation](A3-writeup/student/t5.png)

---

## A3T6 – Direct Lighting

**Requirement (summary):**
Improve quality/speed via direct lighting sampling.

### Result
![Your creation](A3-writeup/student/t6.png)
![Your creation](A3-writeup/student/t62.png)

---

## A3T7 – Environment Lighting

**Requirement (summary):**
Implement environment map importance sampling (CDF/PDF usage, Jacobian, correct parameterization).

### Result
![Your creation](A3-writeup/student/t71.png)
![Your creation](A3-writeup/student/t72.png)

---

## Final Rendered Image

![Your creation](A3-writeup/student/render.png)

Scene description:
I placed my modeled object from previous HW into the Cornell-box-like scene, replacing the two spheres, and tested different materials (mirror/glass) for the surrounding ring.

Additional renders:

![Your creation](A3-writeup/student/render2.png)

![Your creation](A3-writeup/student/render3.png)

---

Author: Chuang Ma
