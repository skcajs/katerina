mod filter;
mod integrator;
mod interval;
mod ray;
mod sampler;
mod sphere;
mod tup;
mod world;

use std::fs::File;
use std::io::Write;
use std::sync::atomic::AtomicUsize;
use std::sync::atomic::Ordering;
use std::time::Instant;

use integrator::integrate;
use integrator::IntegrationType;
use rayon::prelude::*;

use filter::tent_filter;
use ray::Ray;
use sampler::Sampler;
use tup::Tup;
use world::World;

fn clamp(x: f64) -> f64 {
    if x < 0. {
        return 0.;
    } else if x > 1. {
        return 1.;
    }
    x
}

fn to_int(x: f64) -> i32 {
    (clamp(x).powf(1. / 2.2) * 255. + 0.5) as i32
}

fn main() {
    let w = 320;
    let h = 240;
    let num_samples: isize = 5; // will be evaluated to num_samples * 4
    let cam = Ray {
        o: Tup(50. - 50., 52. - 52., 295.6 - 25.),
        d: Tup(0., -0.046, -1.).norm(),
    };

    let cx = Tup(w as f64 * 0.5135 / h as f64, 0.0, 0.0);
    let cy = (cx.cross(cam.d)).norm() * 0.5135;
    let mut data: Vec<(usize, usize, Tup)> = vec![];
    for i in (0..h).rev() {
        for j in 0..w {
            data.push((i, j, Tup(0., 0., 0.)));
        }
    }

    let world = World::new();

    let now = Instant::now();

    let progress_counter = AtomicUsize::new(0);
    let total_pixels = h * w;

    data.par_chunks_mut(100).for_each(|slice| {
        let mut sampler = Sampler::new();
        slice.into_iter().for_each(|p| {
            let y = p.0;
            let x = p.1;
            for sy in 0..2 {
                for sx in 0..2 {
                    let mut rad = Tup(0., 0., 0.);
                    rad = (0..num_samples).into_iter().fold(rad, |acc, _| {
                        let (dx, dy) = tent_filter(&mut sampler);

                        let d = cx * (((sx as f64 + 0.5 + dx) / 2. + x as f64) / w as f64 - 0.5)
                            + cy * (((sy as f64 + 0.5 + dy) / 2. + y as f64) / h as f64 - 0.5)
                            + cam.d;

                        acc + integrate(
                            &world,
                            Ray {
                                o: cam.o + d * 140.,
                                d: d.norm(),
                            },
                            0,
                            &mut sampler,
                            IntegrationType::default(),
                        ) * (1. / num_samples as f64)
                    });

                    p.2 = p.2 + Tup(clamp(rad.0), clamp(rad.1), clamp(rad.2)) * 0.25;
                }
            }
        });

        // Increment progress
        let prev_count = progress_counter.fetch_add(1, Ordering::SeqCst);

        // Print progress
        let progress = 10000. * prev_count as f64 / (total_pixels as f64);
        print!("\rRendering {0} spp {1:.2}%", num_samples, progress);
    });

    let elapsed_time = now.elapsed();
    println!(
        "\nRunning integrator took {} seconds.",
        elapsed_time.as_secs(),
    );

    let mut f = File::create("image.ppm").unwrap();
    writeln!(f, "P3\n{} {}\n255", w, h).unwrap();
    for i in 0..w * h {
        writeln!(
            f,
            "{} {} {}",
            to_int(data[i].2 .0),
            to_int(data[i].2 .1),
            to_int(data[i].2 .2)
        )
        .unwrap();
    }
}
