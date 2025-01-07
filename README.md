# Katerina - A Geodesic Path Tracer

Besed of a smallpt, written in rust. 

https://www.kevinbeason.com/smallpt/#moreinfo

### Rust is required to run this project

You can run this project by running the command:

cargo run --release

!!IMPORTANT!! use --release or you will be waiting a long time ;)

Rayon is required to run this project, but should be installed automatically (or when you run cargo install)

You can modify the image resolution and number of samples in the main function (right at the top of the function).

The IntegrationType can be either Iterative (default) or Recursive. The recursive function is faster, because it also uses light source sampling (try reducing the number of samples).

You can also update the values in the metric, m (under the main function). The parameters are:

s: the position of the centre of the black hole/ mass
a: the angular momentum (useful for the Kerr metric, not fully implemented yet)
m: the mass of the the black hole (and hence the Schwarszchild metric)

The current scene is a simple cornell box, under scenes>cornell.rs. There are currently two variants, the second scene is used by the recursive function, and is needed because the light sphere is originally a massive sphere and only a small portion of it is visible, thus does not work well with light surface sampling. The scenes will swap automatically when switching the IntegrationType.

Images are saves as image.ppm, in the root directory.

## Project Structure

main.rs
The entrypoint is in main.rs, where a ray is initialised and then calls the integrate function in integrator.rs.

integrator.rs
Within the integrator, the ray is tested for an intersection within the world.trace_geodesic function in world.rs, and for a sucessful intersection, the ray is then scattered off its surface, according the the surfaces BSDF. This is an iterative (or recursive process)

world.rs
The trace_geodisic function is the main function in this file, which marches the ray through the geodesic world lines using a solver within the solver.rs file, the current working solver is RK4. There is also a RKf45 method (WIP). The solver returns the rays position and direction at each step. The geodesic tracer then tests for an intersection at close proximity. If an intersection occurs, the function returns true.

solver.rs
This file contains various functions for both the metrics, and the numerical integration of the EOM's (equations of motion). They take the current ray, and step them through a timestep, returning the successive position and momentum of the ray.

sphere.rs
Defines the equations for a sphere, which is currently what the world is comprised of.

tup.rs
Useful vector functions needed.

ray.rs
A ray object, with position and direction, in cartesian form.

geodesic.rs
Think of this as a ray++ object, which contains the ray, and much more. This is primarily needed as the Kerr solution requires working in spherical coordinate systems, and parameters for the angular momentum. This is only useful for the Kerr metric and beyond, which is not fully implemented yet.

metric.rs
Contains information about the Metric, such as its position in the world, the mass (schwarzschild) and the a parameter (kerr).

interval.rs
Redundant, replaced with solver.rs

