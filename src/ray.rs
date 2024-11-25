use nalgebra::Vector4;

use crate::geodesic::Geodesic;
use crate::helpers::{rest_frame_to_zamo, zamo_to_global};

use super::tup::Tup;

#[derive(Clone)]
pub struct Ray {
    pub o: Tup,
    pub d: Tup,
}

impl Ray {
    pub fn compute(&self, cam_props: Geodesic) -> Geodesic {
        let a = cam_props.a;
        let bl_coords = self.o.cartesian_to_boyer_lindquist(a);
        let t = 0.0;
        let r = bl_coords.0;
        let theta = bl_coords.1;
        let phi = bl_coords.2;
        let delta = (r * r) - (2.0 * r) + (a * a);
        let sigma = (r * r) + (a * a) * (f64::cos(theta) * f64::cos(theta));
        let zamo_dir: Vector4<f64> = rest_frame_to_zamo(self.d, cam_props, a);
        let photon_momentum: Vector4<f64> =
            zamo_to_global(zamo_dir, r, theta, phi, delta, sigma, a);

        Geodesic {
            t,
            r,
            theta,
            phi,
            u0: photon_momentum[0],
            u1: photon_momentum[1],
            u2: photon_momentum[2],
            u3: photon_momentum[3],
            h_step: 0.0,
            dx: 0.0,
            a,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_ray() {
        let o = Tup(0.0, 0.0, 0.0);
        let d = Tup(0.3, 0.5, 0.4);

        let r = Ray { o, d };

        assert_eq!(r.o, o);
        assert_eq!(r.d, d)
    }
}
