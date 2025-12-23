# MyScotty3D – Computer Graphics Systems

This repository is built on top of **Scotty3D**, a C++ computer graphics framework originally developed at Carnegie Mellon University for experimenting with modern graphics systems.

Scotty3D repository:  
[MyScotty3D](./MyScotty3D)

---

## Scotty3D Codebase Structure

The Scotty3D framework has the following core directories.

### deps/
Third-party dependencies required by Scotty3D.

These components are treated as external libraries and are not modified directly.

### media/
Scene descriptions, meshes, textures, and example assets used for testing and rendering.  
This includes `.js3d` scene files and imported geometry used across rasterization, path tracing, and animation.

### src/
The core implementation directory of Scotty3D.  
All major graphics systems are implemented or extended here, including:
- Software rasterization pipeline
- Halfedge mesh data structures and geometry processing
- BVH acceleration structures
- Path tracing renderer
- Animation, skinning, and particle simulation systems

### tests/
Unit and integration tests used to validate correctness of geometry operations, rendering algorithms, and animation behavior.

---

## Implemented Graphics Systems

On top of the Scotty3D framework, the following major graphics subsystems were implemented.  
Each folder below contains a dedicated README with task-level descriptions, implementation details, and visual results.

---

## A0 – Scotty3D Setup and Debugging  
[A0-Scotty3D-Setup-and-Debugging](./A0-Scotty3D-Setup-and-Debugging)

Establishes a reliable C++ development and debugging workflow for Scotty3D.

- Build and run Scotty3D in both GUI and CLI modes  
- Debug compilation and runtime issues in a large C++ codebase  
- Diagnose and fix a Halfedge mesh validation crash triggered from the GUI  

---

## A1 – Software Rasterizer  
[A1-Software-Rasterizer](./A1-Software-Rasterizer)

Implements a complete GPU-style rasterization pipeline in software.

- Scene graph transformations  
- Line and triangle rasterization  
- Depth testing and blending  
- Attribute interpolation (flat, smooth, perspective-correct)  
- Texture sampling with mipmapping  
- Supersampling anti-aliasing  

---

## A2 – Halfedge Mesh Editing  
[A2-Halfedge-Mesh-Editing](./A2-Halfedge-Mesh-Editing)

Implements interactive mesh editing and geometry processing using a Halfedge data structure.

- Local topology operations (flip, split, collapse, extrude, bevel)  
- Global geometry processing (triangulation, subdivision)  
- Catmull–Clark and Loop subdivision  
- Robust handling of boundary and degenerate cases  

---

## A3 – Path Tracing Renderer  
[A3-Path-Tracing-Renderer](./A3-Path-Tracing-Renderer)

Builds a physically based path tracing renderer with global illumination.

- Camera ray generation  
- Ray–primitive intersection tests  
- BVH acceleration structures  
- Recursive path tracing integrator  
- Diffuse, mirror, and glass materials  
- Direct and environment lighting with importance sampling  

---

## A4 – Animation, Skinning, and Particles  
[A4-Animation-Skinning-Particles](./A4-Animation-Skinning-Particles)

Completes the graphics pipeline with animation and simulation systems.

- Spline-based animation interpolation  
- Hierarchical skeletal animation (FK and IK)  
- Linear Blend Skinning for deformable meshes  
- Particle simulation with collision handling  
