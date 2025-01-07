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
