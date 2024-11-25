use crate::geodesic::Geodesic;
use crate::helpers::u0;
use crate::tup::Tup;

pub struct Cam {
    pub pos: Tup,
    pub dir: Tup,
}

impl Cam {
    pub fn props(&self, s: Tup, a: f64) -> Geodesic {
        let bl_coords = (self.pos - s).cartesian_to_boyer_lindquist(a);
        let r_cam = bl_coords.0;
        let delta = r_cam.powi(2) - 2. * r_cam + a.powi(2);
        let sigma = r_cam.powi(2) + a.powi(2) * f64::cos(bl_coords.1).powi(2);
        let u0_cam = u0(
            self.dir.0,
            self.dir.1,
            self.dir.2,
            bl_coords.0,
            bl_coords.1,
            delta,
            sigma,
            -1.0,
            1.,
        );
        Geodesic {
            t: 0.,
            r: r_cam,
            theta: bl_coords.1,
            phi: bl_coords.2,
            u0: u0_cam,
            u1: self.dir.0,
            u2: self.dir.1,
            u3: self.dir.2,
            h_step: 0.,
            dx: 0.,
            a,
        }
    }
}
