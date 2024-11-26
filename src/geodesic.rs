use nalgebra::{Matrix4, Vector4};

use crate::{metric::Metric, ray::Ray, tup::Tup};

pub struct Geodesic {
    pub ray: Ray, // For resolving the image
    pub m: Metric,
    pub t: f64,
    pub r: f64,
    pub theta: f64,
    pub phi: f64,
    pub u0: f64,
    pub u1: f64,
    pub u2: f64,
    pub u3: f64,
}

impl Geodesic {
    pub fn ray(pos: Tup, dir: Tup, m: Metric) -> Self {
        let bl_coords = (pos - m.s).cartesian_to_boyer_lindquist(m.a);
        let r_cam = bl_coords.0;
        let delta = r_cam.powi(2) - 2. * r_cam + m.a.powi(2);
        let sigma = r_cam.powi(2) + m.a.powi(2) * f64::cos(bl_coords.1).powi(2);
        let u0_cam = Self::u0(
            dir.0,
            dir.1,
            dir.2,
            bl_coords.0,
            bl_coords.1,
            delta,
            sigma,
            -1.0,
            1.,
        );
        Geodesic {
            ray: Ray { o: pos, d: dir },
            m,
            t: 0.,
            r: r_cam,
            theta: bl_coords.1,
            phi: bl_coords.2,
            u0: u0_cam,
            u1: dir.0,
            u2: dir.1,
            u3: dir.2,
        }
    }

    pub fn init_ray(pos: Tup, dir: Tup, props: &Geodesic) -> Geodesic {
        let a = props.m.a;
        let bl_coords = (pos - props.m.s).cartesian_to_boyer_lindquist(a);
        let t = 0.0;
        let r = bl_coords.0;
        let theta = bl_coords.1;
        let phi = bl_coords.2;
        let delta = (r * r) - (2.0 * r) + (a * a);
        let sigma = (r * r) + (a * a) * (f64::cos(theta) * f64::cos(theta));
        let zamo_dir: Vector4<f64> = Self::rest_frame_to_zamo(dir, props, a);
        let photon_momentum: Vector4<f64> =
            Self::zamo_to_global(zamo_dir, r, theta, phi, delta, sigma, a);

        Geodesic {
            ray: Ray { o: pos, d: dir },
            m: props.m.clone(),
            t,
            r,
            theta,
            phi,
            u0: photon_momentum[0],
            u1: photon_momentum[1],
            u2: photon_momentum[2],
            u3: photon_momentum[3],
        }
    }

    fn rest_frame_to_zamo(ray_dir: Tup, cam_props: &Geodesic, a: f64) -> Vector4<f64> {
        let r = cam_props.r;
        let theta = cam_props.theta;
        let phi = cam_props.phi;
        let delta = (r * r) - (2.0 * r) + (a * a);
        let sigma = r * r + a * a * f64::cos(theta) * f64::cos(theta);
        let u0 = cam_props.u0;
        let u1 = cam_props.u1;
        let u2 = cam_props.u2;
        let u3 = cam_props.u3;

        let w2 = (r * r) + (a * a);
        let lambda = (w2 * w2) - (a * a * delta * (f64::sin(theta) * f64::sin(theta)));
        let u_t = f64::sqrt((delta * sigma) / lambda) * u0;
        let u_r = f64::sqrt(sigma / delta) * u1;
        let u_theta = f64::sqrt(sigma) * u2;
        let u_phi = (f64::sin(theta) * f64::sqrt(lambda / sigma) * u3)
            - (((2.0 * a * r * f64::sin(theta)) / f64::sqrt(lambda * sigma)) * u0);
        let v_r = u_r / u_t;
        let v_theta = u_theta / u_t;
        let v_phi = u_phi / u_t;
        let v_x = (1. / f64::sqrt((r * r) + (a * a) * (f64::cos(theta) * f64::cos(theta))))
            * ((r * f64::sin(theta) * f64::cos(phi) * v_r)
                + (f64::sqrt((r * r) + (a * a)) * f64::cos(theta) * f64::cos(phi) * v_theta))
            - (f64::sin(phi) * v_phi);
        let v_y = (1. / f64::sqrt((r * r) + (a * a) * (f64::cos(theta) * f64::cos(theta))))
            * ((f64::sqrt((r * r) + (a * a)) * f64::cos(theta) * v_r)
                - (r * f64::sin(theta) * v_theta));
        let v_z = (1. / f64::sqrt((r * r) + (a * a) * (f64::cos(theta) * f64::cos(theta))))
            * ((r * f64::sin(theta) * f64::sin(phi) * v_r)
                + (f64::sqrt((r * r) + (a * a)) * f64::cos(theta) * f64::sin(phi) * v_theta))
            + (f64::cos(phi) * v_phi);
        let v_2 = v_x.powi(2) + v_y.powi(2) + v_z.powi(2);
        let gamma = 1.0 / f64::sqrt(1.0 - v_2);
        let m_cam: Matrix4<f64> = Matrix4::new(
            gamma,
            -gamma * v_x,
            -gamma * v_y,
            -gamma * v_z,
            -gamma * v_x,
            1.0 + (gamma - 1.0) * ((v_x * v_x) / (v_2)),
            (gamma - 1.0) * ((v_x * v_y) / v_2),
            (gamma - 1.0) * ((v_x * v_z) / v_2),
            -gamma * v_y,
            (gamma - 1.0) * ((v_x * v_y) / v_2),
            1.0 + (gamma - 1.0) * ((v_y * v_y) / (v_2)),
            (gamma - 1.0) * ((v_y * v_z) / v_2),
            -gamma * v_z,
            (gamma - 1.0) * ((v_x * v_z) / v_2),
            (gamma - 1.0) * ((v_y * v_z) / v_2),
            1.0 + (gamma - 1.0) * ((v_z * v_z) / (v_2)),
        );

        let mut ray4: Vector4<f64> = Vector4::new(
            f64::sqrt(ray_dir.0 * ray_dir.0 + ray_dir.1 * ray_dir.1 + ray_dir.2 * ray_dir.2),
            ray_dir.0,
            ray_dir.1,
            ray_dir.2,
        );
        ray4 = m_cam * ray4;
        Vector4::new(ray4[1], ray4[2], ray4[3], ray4[0])
    }

    fn zamo_to_global(
        ray_dir: Vector4<f64>,
        r: f64,
        theta: f64,
        phi: f64,
        delta: f64,
        sigma: f64,
        a: f64,
    ) -> Vector4<f64> {
        let ray_r = (1.0 / f64::sqrt((r * r) + (a * a) * (f64::cos(theta) * f64::cos(theta))))
            * ((r * f64::sin(theta) * f64::cos(phi) * ray_dir[0])
                + (r * f64::sin(theta) * f64::sin(phi) * ray_dir[2])
                + (f64::sqrt((r * r) + (a * a)) * f64::cos(theta) * ray_dir[1]));
        let ray_theta = (1.0 / f64::sqrt((r * r) + (a * a) * (f64::cos(theta) * f64::cos(theta))))
            * ((f64::sqrt((r * r) + (a * a)) * f64::cos(theta) * f64::cos(phi) * ray_dir[0])
                + (f64::sqrt((r * r) + (a * a)) * f64::cos(theta) * f64::sin(phi) * ray_dir[2])
                - (r * f64::sin(theta) * ray_dir[1]));
        let ray_phi = (-f64::sin(phi) * ray_dir[0]) + (f64::cos(phi) * ray_dir[2]);
        let ray_t = ray_dir[3];
        let w2 = (r * r) + (a * a);
        let lambda = (w2 * w2) - (a * a * delta * (f64::sin(theta) * f64::sin(theta)));
        let u0 = f64::sqrt(lambda / (delta * sigma)) * ray_t
            + (((2.0 * a * r) / f64::sqrt(lambda * delta * sigma)) * ray_phi);
        let u1 = f64::sqrt(delta / sigma) * ray_r;
        let u2 = (1.0 / f64::sqrt(sigma)) * ray_theta;
        let u3 = f64::sqrt(sigma / lambda) * (1.0 / f64::sin(theta)) * ray_phi;

        Vector4::new(u0, u1, u2, u3)
    }

    fn u0(
        u1: f64,
        u2: f64,
        u3: f64,
        r: f64,
        theta: f64,
        delta: f64,
        sigma: f64,
        m_y: f64,
        a: f64,
    ) -> f64 {
        let gt_phi = -((2.0 * a * r) / sigma) * (f64::sin(theta) * f64::sin(theta));
        let gtt = -(1.0 - ((2.0 * r) / sigma));
        let grr = sigma / delta;
        let g_theta_theta = sigma;
        let lambda = (((r * r) + (a * a)) * ((r * r) + (a * a)))
            - (a * a) * delta * (f64::sin(theta) * f64::sin(theta));
        let g_phi_phi = (lambda / sigma) * f64::sin(theta) * f64::sin(theta);
        let mut root_term = (((gt_phi / gtt) * (gt_phi / gtt)) * (u3 * u3))
            - ((1.0 / gtt)
                * (grr * (u1 * u1) + g_theta_theta * (u2 * u2) + g_phi_phi * (u3 * u3) - m_y));

        root_term = if root_term < 0. {
            0.
        } else {
            f64::sqrt(root_term)
        };

        ((-gt_phi / gtt) * u3) + root_term
    }
}
