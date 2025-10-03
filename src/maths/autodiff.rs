#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DualNumber {
    pub f: f64,
    pub d: f64,
}

impl std::fmt::Display for DualNumber {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "f:{}, df:{}", self.f, self.d)
    }
}

impl DualNumber {
    pub fn new(f: f64, d: f64) -> Self {
        DualNumber { f, d }
    }
}

// Operations with other DualNumbers
impl std::ops::Add<&DualNumber> for &DualNumber {
    type Output = DualNumber;

    fn add(self, other: &DualNumber) -> DualNumber {
        DualNumber {
            f: self.f + other.f,
            d: self.d + other.d,
        }
    }
}

impl std::ops::Add for DualNumber {
    type Output = DualNumber;

    fn add(self, other: DualNumber) -> DualNumber {
        &self + &other
    }
}

impl std::ops::Sub<&DualNumber> for &DualNumber {
    type Output = DualNumber;

    fn sub(self, other: &DualNumber) -> DualNumber {
        DualNumber {
            f: self.f - other.f,
            d: self.d - other.d,
        }
    }
}

impl std::ops::Sub for DualNumber {
    type Output = DualNumber;

    fn sub(self, other: DualNumber) -> DualNumber {
        &self - &other
    }
}

impl std::ops::Mul<&DualNumber> for &DualNumber {
    type Output = DualNumber;

    fn mul(self, other: &DualNumber) -> DualNumber {
        DualNumber {
            f: self.f * other.f,
            d: self.f * other.d + self.d * other.f,
        }
    }
}

impl std::ops::Mul for DualNumber {
    type Output = DualNumber;

    fn mul(self, other: DualNumber) -> DualNumber {
        &self * &other
    }
}

impl std::ops::Div<&DualNumber> for &DualNumber {
    type Output = DualNumber;

    fn div(self, other: &DualNumber) -> DualNumber {
        DualNumber {
            f: self.f / other.f,
            d: (self.d * other.f - self.f * other.d) / (other.f * other.f),
        }
    }
}

impl std::ops::Div for DualNumber {
    type Output = DualNumber;

    fn div(self, other: DualNumber) -> DualNumber {
        &self / &other
    }
}

impl std::ops::Neg for &DualNumber {
    type Output = DualNumber;

    fn neg(self) -> DualNumber {
        DualNumber {
            f: -self.f,
            d: -self.d,
        }
    }
}

impl std::ops::Neg for DualNumber {
    type Output = DualNumber;

    fn neg(self) -> DualNumber {
        -&self
    }
}

impl DualNumber {
    pub fn sin(&self) -> DualNumber {
        DualNumber {
            f: self.f.sin(),
            d: self.d * self.f.cos(),
        }
    }

    pub fn cos(&self) -> DualNumber {
        DualNumber {
            f: self.f.cos(),
            d: -self.d * self.f.sin(),
        }
    }

    pub fn exp(&self) -> DualNumber {
        let exp_f = self.f.exp();
        DualNumber {
            f: exp_f,
            d: self.d * exp_f,
        }
    }

    pub fn ln(&self) -> DualNumber {
        DualNumber {
            f: self.f.ln(),
            d: self.d / self.f,
        }
    }

    pub fn powf(&self, n: f64) -> DualNumber {
        DualNumber {
            f: self.f.powf(n),
            d: self.d * n * self.f.powf(n - 1.0),
        }
    }

    pub fn powi(&self, n: i32) -> DualNumber {
        DualNumber {
            f: self.f.powi(n),
            d: self.d * (n as f64) * self.f.powi(n - 1),
        }
    }

    pub fn sqrt(&self) -> DualNumber {
        DualNumber {
            f: self.f.sqrt(),
            d: self.d / (2.0 * self.f.sqrt()),
        }
    }
}

// Operations with f64 on the right
impl std::ops::Add<f64> for &DualNumber {
    type Output = DualNumber;

    fn add(self, other: f64) -> DualNumber {
        DualNumber {
            f: self.f + other,
            d: self.d,
        }
    }
}

impl std::ops::Add<f64> for DualNumber {
    type Output = DualNumber;

    fn add(self, other: f64) -> DualNumber {
        &self + other
    }
}

impl std::ops::Sub<f64> for &DualNumber {
    type Output = DualNumber;

    fn sub(self, other: f64) -> DualNumber {
        DualNumber {
            f: self.f - other,
            d: self.d,
        }
    }
}

impl std::ops::Sub<f64> for DualNumber {
    type Output = DualNumber;

    fn sub(self, other: f64) -> DualNumber {
        &self - other
    }
}

impl std::ops::Mul<f64> for &DualNumber {
    type Output = DualNumber;

    fn mul(self, other: f64) -> DualNumber {
        DualNumber {
            f: self.f * other,
            d: self.d * other,
        }
    }
}

impl std::ops::Mul<f64> for DualNumber {
    type Output = DualNumber;

    fn mul(self, other: f64) -> DualNumber {
        &self * other
    }
}

impl std::ops::Div<f64> for &DualNumber {
    type Output = DualNumber;

    fn div(self, other: f64) -> DualNumber {
        DualNumber {
            f: self.f / other,
            d: self.d / other,
        }
    }
}

impl std::ops::Div<f64> for DualNumber {
    type Output = DualNumber;

    fn div(self, other: f64) -> DualNumber {
        &self / other
    }
}

// Operations with f64 on the left
impl std::ops::Add<&DualNumber> for f64 {
    type Output = DualNumber;

    fn add(self, other: &DualNumber) -> DualNumber {
        DualNumber {
            f: self + other.f,
            d: other.d,
        }
    }
}

impl std::ops::Add<DualNumber> for f64 {
    type Output = DualNumber;

    fn add(self, other: DualNumber) -> DualNumber {
        self + &other
    }
}

impl std::ops::Sub<&DualNumber> for f64 {
    type Output = DualNumber;

    fn sub(self, other: &DualNumber) -> DualNumber {
        DualNumber {
            f: self - other.f,
            d: -other.d,
        }
    }
}

impl std::ops::Sub<DualNumber> for f64 {
    type Output = DualNumber;

    fn sub(self, other: DualNumber) -> DualNumber {
        self - &other
    }
}

impl std::ops::Mul<&DualNumber> for f64 {
    type Output = DualNumber;

    fn mul(self, other: &DualNumber) -> DualNumber {
        DualNumber {
            f: self * other.f,
            d: self * other.d,
        }
    }
}

impl std::ops::Mul<DualNumber> for f64 {
    type Output = DualNumber;

    fn mul(self, other: DualNumber) -> DualNumber {
        self * &other
    }
}

impl std::ops::Div<&DualNumber> for f64 {
    type Output = DualNumber;

    fn div(self, other: &DualNumber) -> DualNumber {
        DualNumber {
            f: self / other.f,
            d: -self * other.d / (other.f * other.f),
        }
    }
}

impl std::ops::Div<DualNumber> for f64 {
    type Output = DualNumber;

    fn div(self, other: DualNumber) -> DualNumber {
        self / &other
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn addition_with_scalar() {
        let x = DualNumber::new(2.0, 1.0);
        let y = x + 3.0;
        assert_eq!(y.f, 5.0);
        assert_eq!(y.d, 1.0);
    }

    #[test]
    fn subtraction_with_scalar() {
        let x = DualNumber::new(2.0, 1.0);
        let y = x - 3.0;
        assert_eq!(y.f, -1.0);
        assert_eq!(y.d, 1.0);
    }

    #[test]
    fn multiplication_with_scalar() {
        let x = DualNumber::new(2.0, 1.0);
        let y = x * 3.0;
        assert_eq!(y.f, 6.0);
        assert_eq!(y.d, 3.0);
    }

    #[test]
    fn division_with_scalar() {
        let x = DualNumber::new(2.0, 1.0);
        let y = x / 2.0;
        assert_eq!(y.f, 1.0);
        assert_eq!(y.d, 0.5);
    }

    #[test]
    fn linear_function() {
        let x = DualNumber::new(3.0, 1.0);
        let y = 2.0 * x + 5.0;
        assert_eq!(y.f, 11.0);
        assert_eq!(y.d, 2.0);
    }

    #[test]
    fn quadratic_function_mult() {
        let x = DualNumber::new(3.0, 1.0);
        let y = x * x + 2.0 * x + 1.0; // (x + 1)^2
        assert_eq!(y.f, 16.0);
        assert_eq!(y.d, 8.0); // dy/dx = 2x + 2 = 8 at x=3
    }

    #[test]
    fn quadratic_function_powi() {
        let x = DualNumber::new(3.0, 1.0);
        let y = x.powi(2) + 2.0 * x + 1.0; // (x + 1)^2
        assert_eq!(y.f, 16.0);
        assert_eq!(y.d, 8.0); // dy/dx = 2x + 2 = 8 at x=3
    }
}
