# Katerina - A Geodesic Path Tracer

**Katerina** is a geodesic path tracer inspired by [Smallpt](https://www.kevinbeason.com/smallpt/#moreinfo), written entirely in Rust.

---

## Requirements

- **Rust**: Ensure Rust is installed on your system. You can install it from [rust-lang.org](https://www.rust-lang.org/).
- **Rayon**: Rayon is used for parallel processing and should be automatically installed when running the project.

---

## Getting Started

### IMPORTANT! Make sure you're on the 'all-broken' branch, conveniently named. TODO: Merge in. 

Run the project using the following command:  

```bash
cargo run --release  
```

**Important:** Always use the `--release` flag for optimized performance; without it, render times will be significantly longer.

### Configuration Options

You can modify the following parameters directly in the `main` function:

- **Image Resolution**: Adjust the resolution of the rendered image.
- **Number of Samples**: Configure the number of samples per pixel.

Additionally, you can set the integration type to either:  
- **Iterative** (default)  
- **Recursive** (faster due to light source sampling).  

For the recursive option, you can reduce the number of samples for quicker results.

---

## Metric Parameters

The `metric` parameters in the `main` function allow you to adjust properties related to the black hole simulation:

- **`s`**: Position of the black hole/mass center.
- **`a`**: Angular momentum (relevant for the Kerr metric; not fully implemented yet).
- **`m`**: Mass of the black hole (used for the Schwarzschild metric).

---

## Scene Description

The default scene is a simple Cornell box located in `scenes/cornell.rs`.  

Two scene variants are available:
1. **Default Scene**: Works with the iterative method.
2. **Light Surface Sampling Scene**: Required for the recursive method, where the light sphere is modified for better light surface sampling.

The appropriate scene is automatically selected based on the chosen integration type.  

Rendered images are saved as `image.ppm` in the root directory.

---

## Project Structure

### Main Files:
- **`main.rs`**  
  The entry point where rays are initialized and the `integrate` function from `integrator.rs` is invoked.

- **`integrator.rs`**  
  Handles ray tracing by testing intersections using the `world.trace_geodesic` function from `world.rs`.  
  For successful intersections, rays are scattered according to the BSDF of the surface.

- **`world.rs`**  
  Contains the `trace_geodesic` function, which marches rays through geodesic world lines using solvers from `solver.rs`.  
  The RK4 solver is currently implemented, while RKF45 is a work in progress.

- **`solver.rs`**  
  Provides functions for numerical integration of equations of motion (EOMs) and geodesic metric calculations.  
  It steps rays through the simulation and returns their successive positions and momenta.

### Additional Files:
- **`sphere.rs`**  
  Defines equations for spheres, which currently compose the entire world.

- **`tup.rs`**  
  Utility functions for vector operations.

- **`ray.rs`**  
  Defines a ray object, including its position and direction in Cartesian coordinates.

- **`geodesic.rs`**  
  Extends the ray object with additional data required for spherical coordinate systems and angular momentum (useful for the Kerr metric, not yet fully implemented).

- **`metric.rs`**  
  Stores information about the metric, including its position, mass (`m` for Schwarzschild), and angular momentum (`a` for Kerr).

- **`interval.rs`**  
  Previously used for solvers; now replaced by `solver.rs`.

---

## Notes

- The project uses the RK4 solver for geodesic calculations by default.
- Light source sampling is only available for the recursive integration type.
- The Kerr metric implementation is almost complete but will be implemented in future updates.

---

## Licence

The source is released under the MIT license, which is open and compatible with GPL. See LICENSE.txt in the distribution for more information.
