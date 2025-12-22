# A1 – Software Rasterization Pipeline

## Overview
Implements a complete software rasterization pipeline in C++, mirroring the core stages of a modern GPU.

---

## A1T1 – Scene Transformations

**Requirement (summary):**
Implement `local_to_world` and `world_to_local` for hierarchical transforms.

### Result

![task1](A1-writeup/student/task1.png)

![task1](A1-writeup/student/test1.png)

---

## A1T2 – Line Rasterization

**Requirement (summary):**
Rasterize wireframe lines (diamond-exit rule), handle endpoint edge cases, pass tests.

### Result

![Rasterize cubes wireframe](A1-writeup/student/task2.1.png)

All red:

![Rasterize cubes wireframe](A1-writeup/student/task2.2.png)

With different color:

![Rasterize cubes wireframe](A1-writeup/student/task2.3.png)

Pass all tests including Piazza's:

![Rasterize cubes wireframe](A1-writeup/student/test2.png)

---

## A1T3 – Flat Triangle Rasterization

**Requirement (summary):**
Rasterize flat-shaded triangles efficiently; handle shared edges with top-left rule.

### Result

![Rasterize flat cubes](A1-writeup/student/task3.1.png)

![Rasterize flat cubes](A1-writeup/student/task3.2.png)

![Rasterize flat cubes](A1-writeup/student/test3.png)

---

## A1T4 – Depth Testing & Blending

**Requirement (summary):**
Implement depth-less testing and blending (Add, Over). Answer blend/depth questions.

### Result

Three colors (one valid configuration):
- Red: Blend Replace, Depth Always
- Green: Blend Replace, Depth Always
- Blue: Blend Replace, Depth Always

![task4](A1-writeup/student/task4.3.png)

Five colors (one valid configuration):
- Red: Blend Add, Depth Always
- Green: Blend Add, Depth Always
- Blue: Blend Add, Depth Always

![task4](A1-writeup/student/task4.4.png)

Additional results:

![task4](A1-writeup/student/task4.2.png)

![task4](A1-writeup/student/test4.png)

---

## A1T5 – Attribute Interpolation (Smooth & Perspective-Correct)

**Requirement (summary):**
Implement screen-space + perspective-correct interpolation and compute derivatives.

### Result

![Rasterize smooth and correct cubes](A1-writeup/student/task5.3.png)

![Rasterize smooth and correct cubes](A1-writeup/student/task5.2.png)

![Rasterize smooth and correct cubes](A1-writeup/student/test5.png)

---

## A1T6 – Mipmapping

**Requirement (summary):**
Implement LOD, mip generation, and texture sampling.

### Result

![Rasterize smooth and correct cubes](A1-writeup/student/task6.png)

![Rasterize smooth and correct cubes](A1-writeup/student/test6.png)

---

## A1T7 – Supersampling

**Requirement (summary):**
Support multi-sample framebuffer indexing, resolve, and a custom sample pattern.

### Result

Sample pattern:

![Sample pattern](A1-writeup/student/task7.1.png)

Explanation:
I alternate sampling points between odd/even sequences (5 points vs 4 points). Distances between points stay constant, but edge-to-boundary distances differ. This may help irregular shapes but can be less effective on highly regular images.

Final test:

![Rasterize smooth and correct cubes](A1-writeup/student/test7.png)

---

## Final Rasterized Image

![Your creation](A1-writeup/student/render2.png)

Description:
I tried to create a sci-fi teleportation point with a surrounding halo. The final result looks closer to a cylindrical space station with a ring, built from two cylinders plus an outer ring.

---

Author: Chuang Ma
