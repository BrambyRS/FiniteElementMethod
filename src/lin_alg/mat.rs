/*
Implements a basic 2D matrix struct and some fundamental operations.
*/

#[derive(Debug)]
pub struct Matrix<T> {
    data: Vec<T>,
    dim: (usize, usize),
}

impl Matrix<f64> {
    // Basic constructor, getters and setters
    pub fn new(dim: (usize, usize)) -> Self {
        let data: Vec<f64> = vec![0.0; dim.0 * dim.1];
        return Self { data, dim };
    }

    pub fn identity(size: usize) -> Self {
        let mut mat: Matrix<f64> = Matrix::<f64>::new((size, size));
        for i in 0..size {
            mat.set(i, i, 1.0);
        }
        return mat;
    }

    pub fn set(&mut self, r: usize, c: usize, val: f64) {
        if r >= self.dim.0 || c >= self.dim.1 {
            panic!("Index out of bounds.");
        }
        self.data[r * self.dim.1 + c] = val;
    }

    pub fn get(&self, r: usize, c: usize) -> f64 {
        if r >= self.dim.0 || c >= self.dim.1 {
            panic!("Index out of bounds.");
        }
        return self.data[r * self.dim.1 + c];
    }

    pub fn get_dim(&self) -> (usize, usize) {
        return self.dim;
    }

    // Matrix specific mathematical operations
    pub fn transpose(&self) -> Self {
        let mut result: Matrix<f64> = Matrix::<f64>::new((self.dim.1, self.dim.0));
        for r in 0..self.dim.0 {
            for c in 0..self.dim.1 {
                result.set(c, r, self.get(r, c));
            }
        }
        return result;
    }

    pub fn dot_product(&self, rhs: &Self) -> f64 {
        assert_eq!(self.dim, rhs.dim);
        let mut sum: f64 = 0.0;
        for r in 0..self.dim.0 {
            for c in 0..self.dim.1 {
                sum += self.get(r, c) * rhs.get(r, c);
            }
        }
        return sum;
    }

    pub fn norm(&self) -> f64 {
        return self.dot_product(self).sqrt();
    }

    pub fn abs(&self) -> Matrix<f64> {
        let mut result: Matrix<f64> = Matrix::<f64>::new(self.dim);

        for r in 0..self.dim.0 {
            for c in 0..self.dim.1 {
                result.set(r, c, self.get(r, c).abs())
            }
        }

        return result;
    }
}

// Mathematical operations overloading
// Reference + Reference
impl std::ops::Add<&Matrix<f64>> for &Matrix<f64> {
    type Output = Matrix<f64>;

    fn add(self, rhs: &Matrix<f64>) -> Self::Output {
        assert_eq!(
            self.dim, rhs.dim,
            "Matrix dimensions must match for addition."
        );
        let mut result: Matrix<f64> = Matrix::<f64>::new(self.dim);
        for r in 0..self.dim.0 {
            for c in 0..self.dim.1 {
                result.set(r, c, self.get(r, c) + rhs.get(r, c));
            }
        }
        return result;
    }
}

// Owned + Owned
impl std::ops::Add for Matrix<f64> {
    type Output = Matrix<f64>;
    fn add(self, rhs: Self) -> Self::Output {
        &self + &rhs
    }
}

// Scalar + Reference (Scalar on the left)
impl std::ops::Add<&Matrix<f64>> for f64 {
    type Output = Matrix<f64>;

    fn add(self, rhs: &Matrix<f64>) -> Self::Output {
        let mut result: Matrix<f64> = Matrix::<f64>::new(rhs.dim);
        for r in 0..rhs.dim.0 {
            for c in 0..rhs.dim.1 {
                result.set(r, c, self + rhs.get(r, c));
            }
        }
        return result;
    }
}

// Reference + Scalar (Scalar on the right)
impl std::ops::Add<f64> for &Matrix<f64> {
    type Output = Matrix<f64>;

    fn add(self, rhs: f64) -> Self::Output {
        rhs + self
    }
}

// Reference - Reference
impl std::ops::Sub<&Matrix<f64>> for &Matrix<f64> {
    type Output = Matrix<f64>;

    fn sub(self, rhs: &Matrix<f64>) -> Self::Output {
        assert_eq!(
            self.dim, rhs.dim,
            "Matrix dimensions must match for subtraction."
        );
        let mut result: Matrix<f64> = Matrix::<f64>::new(self.dim);
        for r in 0..self.dim.0 {
            for c in 0..self.dim.1 {
                result.set(r, c, self.get(r, c) - rhs.get(r, c));
            }
        }
        return result;
    }
}

// Owned - Owned
impl std::ops::Sub for Matrix<f64> {
    type Output = Matrix<f64>;
    fn sub(self, rhs: Self) -> Self::Output {
        &self - &rhs
    }
}

// Scalar - Reference (Scalar on the left)
impl std::ops::Sub<&Matrix<f64>> for f64 {
    type Output = Matrix<f64>;

    fn sub(self, rhs: &Matrix<f64>) -> Self::Output {
        let mut result: Matrix<f64> = Matrix::<f64>::new(rhs.dim);
        for r in 0..rhs.dim.0 {
            for c in 0..rhs.dim.1 {
                result.set(r, c, self - rhs.get(r, c));
            }
        }
        return result;
    }
}

// Reference - Scalar (Scalar on the right)
impl std::ops::Sub<f64> for &Matrix<f64> {
    type Output = Matrix<f64>;

    fn sub(self, rhs: f64) -> Self::Output {
        let mut result: Matrix<f64> = Matrix::<f64>::new(self.dim);
        for r in 0..self.dim.0 {
            for c in 0..self.dim.1 {
                result.set(r, c, self.get(r, c) - rhs);
            }
        }
        return result;
    }
}

// Reference * Reference
impl std::ops::Mul<&Matrix<f64>> for &Matrix<f64> {
    type Output = Matrix<f64>;

    fn mul(self, rhs: &Matrix<f64>) -> Self::Output {
        assert!(self.dim.1 == rhs.dim.0);
        let mut result: Matrix<f64> = Matrix::<f64>::new((self.dim.0, rhs.dim.1));
        for r1 in 0..self.dim.0 {
            for c2 in 0..rhs.dim.1 {
                let mut sum = 0.0;
                for c1 in 0..self.dim.1 {
                    sum += self.get(r1, c1) * rhs.get(c1, c2);
                }
                result.set(r1, c2, sum);
            }
        }
        return result;
    }
}

// Owned * Owned (forward to reference version)
impl std::ops::Mul for Matrix<f64> {
    type Output = Matrix<f64>;
    fn mul(self, rhs: Self) -> Self::Output {
        &self * &rhs
    }
}

// Scalar * Reference (Scalar on the left)
impl std::ops::Mul<&Matrix<f64>> for f64 {
    type Output = Matrix<f64>;
    fn mul(self, rhs: &Matrix<f64>) -> Self::Output {
        let mut result: Matrix<f64> = Matrix::<f64>::new(rhs.dim);
        for r in 0..rhs.dim.0 {
            for c in 0..rhs.dim.1 {
                result.set(r, c, self * rhs.get(r, c));
            }
        }
        return result;
    }
}

// Reference * Scalar (Scalar on the right)
impl std::ops::Mul<f64> for &Matrix<f64> {
    type Output = Matrix<f64>;
    fn mul(self, rhs: f64) -> Self::Output {
        rhs * self
    }
}

// Reference == Reference
impl PartialEq for Matrix<f64> {
    fn eq(&self, other: &Self) -> bool {
        if self.dim != other.dim {
            return false;
        }
        for r in 0..self.dim.0 {
            for c in 0..self.dim.1 {
                if self.get(r, c) != other.get(r, c) {
                    return false;
                }
            }
        }
        return true;
    }
}

// Reference / Scalar
impl std::ops::Div<f64> for &Matrix<f64> {
    type Output = Matrix<f64>;
    fn div(self, rhs: f64) -> Self::Output {
        let mut result: Matrix<f64> = Matrix::<f64>::new(self.dim);
        for r in 0..self.dim.0 {
            for c in 0..self.dim.1 {
                result.set(r, c, self.get(r, c) / rhs);
            }
        }
        return result;
    }
}

// Data management operations
impl Clone for Matrix<f64> {
    fn clone(&self) -> Self {
        let mut new_data: Vec<f64> = vec![0.0; self.dim.0 * self.dim.1];
        new_data.clone_from_slice(&self.data);
        return Self {
            data: new_data,
            dim: self.dim,
        };
    }
}

// Tests
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_matrix_creation() {
        let m: Matrix<f64> = Matrix::<f64>::new((2, 3));
        assert_eq!(m.get_dim(), (2, 3));
        for r in 0..2 {
            for c in 0..3 {
                assert_eq!(m.get(r, c), 0.0);
            }
        }
    }

    #[test]
    fn test_identity_matrix() {
        let m: Matrix<f64> = Matrix::<f64>::identity(3);
        assert_eq!(m.get_dim(), (3, 3));
        for r in 0..3 {
            for c in 0..3 {
                if r == c {
                    assert_eq!(m.get(r, c), 1.0);
                } else {
                    assert_eq!(m.get(r, c), 0.0);
                }
            }
        }
    }

    #[test]
    fn test_matrix_get_set() {
        let mut m: Matrix<f64> = Matrix::<f64>::new((2, 2));
        m.set(0, 0, 1.0);
        m.set(0, 1, 2.0);
        m.set(1, 0, 3.0);
        m.set(1, 1, 4.0);
        assert_eq!(m.get(0, 0), 1.0);
        assert_eq!(m.get(0, 1), 2.0);
        assert_eq!(m.get(1, 0), 3.0);
        assert_eq!(m.get(1, 1), 4.0);
    }

    // TODO: Add tests for out of bounds

    #[test]
    fn test_matrix_transpose() {
        let mut m: Matrix<f64> = Matrix::<f64>::new((2, 3));
        m.set(0, 0, 1.0);
        m.set(0, 1, 2.0);
        m.set(0, 2, 3.0);
        m.set(1, 0, 4.0);
        m.set(1, 1, 5.0);
        m.set(1, 2, 6.0);

        let mt: Matrix<f64> = m.transpose();
        assert_eq!(mt.get_dim(), (3, 2));
        assert_eq!(mt.get(0, 0), 1.0);
        assert_eq!(mt.get(0, 1), 4.0);
        assert_eq!(mt.get(1, 0), 2.0);
        assert_eq!(mt.get(1, 1), 5.0);
        assert_eq!(mt.get(2, 0), 3.0);
        assert_eq!(mt.get(2, 1), 6.0);
    }

    #[test]
    fn test_matrix_dot_product() {
        let mut m1: Matrix<f64> = Matrix::<f64>::new((2, 2));
        m1.set(0, 0, 1.0);
        m1.set(0, 1, 2.0);
        m1.set(1, 0, 3.0);
        m1.set(1, 1, 4.0);

        let mut m2: Matrix<f64> = Matrix::<f64>::new((2, 2));
        m2.set(0, 0, 5.0);
        m2.set(0, 1, 6.0);
        m2.set(1, 0, 7.0);
        m2.set(1, 1, 8.0);

        let dp: f64 = m1.dot_product(&m2);
        assert_eq!(dp, 70.0); // 1*5 + 2*6 + 3*7 + 4*8 = 70
    }

    #[test]
    fn test_matrix_norm() {
        let mut m: Matrix<f64> = Matrix::<f64>::new((2, 2));
        m.set(0, 0, 3.0);
        m.set(0, 1, 4.0);
        m.set(1, 0, 0.0);
        m.set(1, 1, 0.0);

        let norm: f64 = m.norm();
        assert_eq!(norm, 5.0); // sqrt(3^2 + 4^2 + 0^2 + 0^2) = 5
    }

    #[test]
    fn test_matrix_abs() {
        let mut m: Matrix<f64> = Matrix::<f64>::new((2, 2));
        m.set(0, 0, 0.0);
        m.set(0, 1, 1.0);
        m.set(1, 0, -2.0);
        m.set(1, 1, 3.0);

        let mut m_expected: Matrix<f64> = Matrix::<f64>::new((2, 2));
        m_expected.set(0, 0, 0.0);
        m_expected.set(0, 1, 1.0);
        m_expected.set(1, 0, 2.0);
        m_expected.set(1, 1, 3.0);

        assert_eq!(m.abs(), m_expected);
    }

    #[test]
    fn test_matrix_add() {
        let mut m1: Matrix<f64> = Matrix::<f64>::new((2, 2));
        m1.set(0, 0, 1.0);
        m1.set(0, 1, 2.0);
        m1.set(1, 0, 3.0);
        m1.set(1, 1, 4.0);

        let mut m2: Matrix<f64> = Matrix::<f64>::new((2, 2));
        m2.set(0, 0, 5.0);
        m2.set(0, 1, 6.0);
        m2.set(1, 0, 7.0);
        m2.set(1, 1, 8.0);

        let mut expected: Matrix<f64> = Matrix::<f64>::new((2, 2));
        expected.set(0, 0, 6.0);
        expected.set(0, 1, 8.0);
        expected.set(1, 0, 10.0);
        expected.set(1, 1, 12.0);

        let result: Matrix<f64> = m1 + m2;
        assert_eq!(result, expected);
    }

    #[test]
    fn test_scalar_add() {
        let mut m1: Matrix<f64> = Matrix::<f64>::new((2, 2));
        m1.set(0, 0, 1.0);
        m1.set(0, 1, 2.0);
        m1.set(1, 0, 3.0);
        m1.set(1, 1, 4.0);

        let mut expected: Matrix<f64> = Matrix::<f64>::new((2, 2));
        expected.set(0, 0, 11.0);
        expected.set(0, 1, 12.0);
        expected.set(1, 0, 13.0);
        expected.set(1, 1, 14.0);

        let result1: Matrix<f64> = &m1 + 10.0;
        assert_eq!(result1, expected);

        let result2: Matrix<f64> = 10.0 + &m1;
        assert_eq!(result2, expected);
    }

    #[test]
    fn test_matrix_sub() {
        let mut m1: Matrix<f64> = Matrix::<f64>::new((2, 2));
        m1.set(0, 0, 5.0);
        m1.set(0, 1, 6.0);
        m1.set(1, 0, 7.0);
        m1.set(1, 1, 8.0);

        let mut m2: Matrix<f64> = Matrix::<f64>::new((2, 2));
        m2.set(0, 0, 1.0);
        m2.set(0, 1, 2.0);
        m2.set(1, 0, 3.0);
        m2.set(1, 1, 4.0);

        let mut expected: Matrix<f64> = Matrix::<f64>::new((2, 2));
        expected.set(0, 0, 4.0);
        expected.set(0, 1, 4.0);
        expected.set(1, 0, 4.0);
        expected.set(1, 1, 4.0);

        let result: Matrix<f64> = m1 - m2;
        assert_eq!(result, expected);
    }

    #[test]
    fn test_scalar_sub() {
        let mut m1: Matrix<f64> = Matrix::<f64>::new((2, 2));
        m1.set(0, 0, 1.0);
        m1.set(0, 1, 2.0);
        m1.set(1, 0, 3.0);
        m1.set(1, 1, 4.0);

        let mut expected1: Matrix<f64> = Matrix::<f64>::new((2, 2));
        expected1.set(0, 0, -9.0);
        expected1.set(0, 1, -8.0);
        expected1.set(1, 0, -7.0);
        expected1.set(1, 1, -6.0);

        let mut expected2: Matrix<f64> = Matrix::<f64>::new((2, 2));
        expected2.set(0, 0, 9.0);
        expected2.set(0, 1, 8.0);
        expected2.set(1, 0, 7.0);
        expected2.set(1, 1, 6.0);

        let result1: Matrix<f64> = &m1 - 10.0;
        assert_eq!(result1, expected1);

        let result2: Matrix<f64> = 10.0 - &m1;
        assert_eq!(result2, expected2);
    }

    #[test]
    fn test_square_mat_mul() {
        let mut m1: Matrix<f64> = Matrix::<f64>::new((2, 2));
        m1.set(0, 0, 1.0);
        m1.set(0, 1, 2.0);
        m1.set(1, 0, 3.0);
        m1.set(1, 1, 4.0);

        let mut m2: Matrix<f64> = Matrix::<f64>::new((2, 2));
        m2.set(0, 0, 5.0);
        m2.set(0, 1, 6.0);
        m2.set(1, 0, 7.0);
        m2.set(1, 1, 8.0);

        let mut expected: Matrix<f64> = Matrix::<f64>::new((2, 2));
        expected.set(0, 0, 19.0);
        expected.set(0, 1, 22.0);
        expected.set(1, 0, 43.0);
        expected.set(1, 1, 50.0);

        let result: Matrix<f64> = m1 * m2;
        assert_eq!(result, expected);
    }

    #[test]
    fn test_rect_mat_mul() {
        let mut m1: Matrix<f64> = Matrix::<f64>::new((2, 3));
        m1.set(0, 0, 1.0);
        m1.set(0, 1, 2.0);
        m1.set(0, 2, 3.0);
        m1.set(1, 0, 4.0);
        m1.set(1, 1, 5.0);
        m1.set(1, 2, 6.0);

        let mut m2: Matrix<f64> = Matrix::<f64>::new((3, 2));
        m2.set(0, 0, 7.0);
        m2.set(0, 1, 8.0);
        m2.set(1, 0, 9.0);
        m2.set(1, 1, 10.0);
        m2.set(2, 0, 11.0);
        m2.set(2, 1, 12.0);

        let mut expected: Matrix<f64> = Matrix::<f64>::new((2, 2));
        expected.set(0, 0, 58.0);
        expected.set(0, 1, 64.0);
        expected.set(1, 0, 139.0);
        expected.set(1, 1, 154.0);

        let result: Matrix<f64> = m1 * m2;
        assert_eq!(result, expected);
    }

    #[test]
    fn test_scalar_mul() {
        let mut m1: Matrix<f64> = Matrix::<f64>::new((2, 2));
        m1.set(0, 0, 1.0);
        m1.set(0, 1, 2.0);
        m1.set(1, 0, 3.0);
        m1.set(1, 1, 4.0);

        let mut expected: Matrix<f64> = Matrix::<f64>::new((2, 2));
        expected.set(0, 0, 10.0);
        expected.set(0, 1, 20.0);
        expected.set(1, 0, 30.0);
        expected.set(1, 1, 40.0);

        let result1: Matrix<f64> = &m1 * 10.0;
        assert_eq!(result1, expected);

        let result2: Matrix<f64> = 10.0 * &m1;
        assert_eq!(result2, expected);
    }

    #[test]
    fn test_scalar_div() {
        let mut m1: Matrix<f64> = Matrix::<f64>::new((2, 2));
        m1.set(0, 0, 10.0);
        m1.set(0, 1, 20.0);
        m1.set(1, 0, 30.0);
        m1.set(1, 1, 40.0);

        let mut expected: Matrix<f64> = Matrix::<f64>::new((2, 2));
        expected.set(0, 0, 1.0);
        expected.set(0, 1, 2.0);
        expected.set(1, 0, 3.0);
        expected.set(1, 1, 4.0);

        let result: Matrix<f64> = &m1 / 10.0;
        assert_eq!(result, expected);
    }

    #[test]
    fn test_matrix_clone() {
        let mut m1: Matrix<f64> = Matrix::<f64>::new((2, 2));
        m1.set(0, 0, 1.0);
        m1.set(0, 1, 2.0);
        m1.set(1, 0, 3.0);
        m1.set(1, 1, 4.0);

        let m2: Matrix<f64> = m1.clone();
        assert_eq!(m2.get_dim(), (2, 2));
        assert_eq!(m1, m2);
    }

    #[test]
    fn test_change_clone() {
        let mut m1: Matrix<f64> = Matrix::<f64>::new((2, 2));
        m1.set(0, 0, 1.0);
        m1.set(0, 1, 2.0);
        m1.set(1, 0, 3.0);
        m1.set(1, 1, 4.0);

        let mut m2: Matrix<f64> = m1.clone();
        m2.set(0, 0, 5.0);
        assert_eq!(m1.get(0, 0), 1.0); // Ensure original matrix is unchanged
        assert_eq!(m2.get(0, 0), 5.0); // Ensure cloned matrix is changed
    }
}
