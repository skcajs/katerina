use crate::{
    sphere::{RflType, Sphere},
    tup::Tup,
};

pub fn cornell_box() -> Vec<Sphere> {
    vec![
        // Scene: radius, position, emission, color, material
        Sphere::new(
            1e5,
            Tup(1e5 - 49., -11.2, 81.6 - 25.),
            Tup::zeros(),
            Tup(0.75, 0.25, 0.25),
            RflType::DIFF,
        ), // Left
        Sphere::new(
            1e5,
            Tup(-1e5 + 49., -11.2, 81.6 - 25.),
            Tup::zeros(),
            Tup(0.25, 0.25, 0.75),
            RflType::DIFF,
        ), // Right
        Sphere::new(
            1e5,
            Tup(0., -11.2, 1e5 - 25.),
            Tup::zeros(),
            Tup(0.75, 0.75, 0.75), //Tup(0.25, 0.75, 0.75),
            RflType::DIFF,
        ), // Back
        Sphere::new(
            1e5,
            Tup(0., -11.2, -1e5 + 145.),
            Tup::zeros(),
            Tup::zeros(), //Tup(0.75, 0.75, 0.25),
            RflType::DIFF,
        ), // Front
        Sphere::new(
            1e5,
            Tup(0., 1e5 - 52., 81.6 - 25.),
            Tup::zeros(),
            Tup(0.75, 0.75, 0.75),
            RflType::DIFF,
        ), // Bottom
        Sphere::new(
            1e5,
            Tup(0., -1e5 + 29.6, 81.6 - 25.),
            Tup::zeros(),
            Tup(0.75, 0.75, 0.75),
            RflType::DIFF,
        ), // Top
        Sphere::new(
            16.5,
            Tup(-23., -35.5, 47.0 - 25.),
            Tup::zeros(),
            Tup(1., 1., 1.) * 0.999,
            RflType::SPEC,
        ), // Mirror
        Sphere::new(
            16.5,
            Tup(23., -35.5, 78. - 25.),
            Tup::zeros(),
            Tup(1., 1., 1.) * 0.999,
            RflType::REFR,
        ), // Glass
        Sphere::new(
            600.,
            Tup(0., 629.6 - 0.27, 81.6 - 25.),
            Tup(12., 12., 12.),
            Tup::zeros(),
            RflType::DIFF,
        ), // Light
           // Sphere::new(
           //     4.,
           //     Tup(-1., -13.2, 0.),
           //     Tup(12., 12., 12.),
           //     Tup::zeros(),
           //     RflType::DIFF,
           // ), // einstein ring
    ]
}

pub fn cornell_box_recursive() -> Vec<Sphere> {
    vec![
        // Scene: radius, position, emission, color, material
        Sphere::new(
            1e5,
            Tup(1e5 - 49., -11.2, 81.6 - 25.),
            Tup::zeros(),
            Tup(0.75, 0.25, 0.25),
            RflType::DIFF,
        ), // Left
        Sphere::new(
            1e5,
            Tup(-1e5 + 49., -11.2, 81.6 - 25.),
            Tup::zeros(),
            Tup(0.25, 0.25, 0.75),
            RflType::DIFF,
        ), // Right
        Sphere::new(
            1e5,
            Tup(0., -11.2, 1e5 - 25.),
            Tup::zeros(),
            Tup(0.75, 0.75, 0.75), //Tup(0.25, 0.75, 0.75),
            RflType::DIFF,
        ), // Back
        Sphere::new(
            1e5,
            Tup(0., -11.2, -1e5 + 170. - 25.),
            Tup::zeros(),
            Tup::zeros(), //Tup(0.75, 0.75, 0.25),
            RflType::DIFF,
        ), // Front
        Sphere::new(
            1e5,
            Tup(0., 1e5 - 52., 81.6 - 25.),
            Tup::zeros(),
            Tup(0.75, 0.75, 0.75),
            RflType::DIFF,
        ), // Bottom
        Sphere::new(
            1e5,
            Tup(0., -1e5 + 29.6, 81.6 - 25.),
            Tup::zeros(),
            Tup(0.75, 0.75, 0.75),
            RflType::DIFF,
        ), // Top
        Sphere::new(
            16.5,
            Tup(-23., -35.5, 47.0 - 25.),
            Tup::zeros(),
            Tup(1., 1., 1.) * 0.999,
            RflType::SPEC,
        ), // Mirror
        Sphere::new(
            16.5,
            Tup(23., -35.5, 78. - 25.),
            Tup::zeros(),
            Tup(1., 1., 1.) * 0.999,
            RflType::REFR,
        ), // Glass
        Sphere::new(
            600.,
            Tup(0., 629.6 - 0.27, 81.6 - 25.),
            Tup(12., 12., 12.),
            Tup::zeros(),
            RflType::DIFF,
        ), // Light
    ]
}
