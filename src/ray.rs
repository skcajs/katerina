use super::tup::Tup;

#[derive(Clone, Debug)]
pub struct Ray {
    pub o: Tup,
    pub d: Tup,
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
