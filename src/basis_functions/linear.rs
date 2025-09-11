#![allow(non_snake_case)]
#![allow(dead_code)]

/*
Linear basis functions and their derivatives for 1D elements.
*/

pub fn N1(xi: f64) -> f64 {
    return (1.0 - xi) / 2.0;
}

pub fn dN1(_xi: f64) -> f64 {
    // xi argument is just for consistency
    return -0.5;
}

pub fn N2(xi: f64) -> f64 {
    return (1.0 + xi) / 2.0;
}

pub fn dN2(_xi: f64) -> f64 {
    // xi argument is just for consistency
    return 0.5;
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_N1() {
        assert_eq!(N1(-1.0), 1.0);
        assert_eq!(N1(0.0), 0.5);
        assert_eq!(N1(1.0), 0.0);
    }

    #[test]
    fn test_dN1() {
        assert_eq!(dN1(-1.0), -0.5);
        assert_eq!(dN1(0.0), -0.5);
        assert_eq!(dN1(1.0), -0.5);
    }

    #[test]
    fn test_N2() {
        assert_eq!(N2(-1.0), 0.0);
        assert_eq!(N2(0.0), 0.5);
        assert_eq!(N2(1.0), 1.0);
    }

    #[test]
    fn test_dN2() {
        assert_eq!(dN2(-1.0), 0.5);
        assert_eq!(dN2(0.0), 0.5);
        assert_eq!(dN2(1.0), 0.5);
    }
}
