/*
Implements a basic 2D matrix struct and some fundamental operations.
*/

struct Matrix<T> {
    data: Vec<T>,
    dim: (usize, usize),
}

impl Matrix<f64> {
    // Basic constructor, getters and setters
    fn new(dim: (usize, usize)) -> Self {
        let data: Vec<f64> = vec![0.0; dim.0 * dim.1];
        return Self { data, dim };
    }

    fn set(&mut self, r: usize, c: usize, val: f64) {
        if r >= self.dim.0 || c >= self.dim.1 {
            panic!("Index out of bounds.");
        }
        self.data[r * self.dim.1 + c] = val;
    }

    fn get(&self, r: usize, c: usize) -> f64 {
        if r >= self.dim.0 || c >= self.dim.1 {
            panic!("Index out of bounds.");
        }
        return self.data[r * self.dim.1 + c];
    }

    fn get_dim(&self) -> (usize, usize) {
        return self.dim;
    }

    // Matrix specific mathematical operations
    fn transpose(&self) -> Self {
        let mut result: Matrix<f64> = Matrix::<f64>::new((self.dim.1, self.dim.0));
        for r in 0..self.dim.0 {
            for c in 0..self.dim.1 {
                result.set(c, r, self.get(r, c));
            }
        }
        return result;
    }

    fn dot_product(&self, rhs: &Self) -> f64 {
        assert_eq!(self.dim, rhs.dim);
        let mut sum: f64 = 0.0;
        for r in 0..self.dim.0 {
            for c in 0..self.dim.1 {
                sum += self.get(r, c) * rhs.get(r, c);
            }
        }
        return sum;
    }
}

// Mathematical operations overloading
impl std::ops::Add for Matrix<f64> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        assert_eq!(self.dim, rhs.dim, "Matrix dimensions must match for addition.");
        let mut result: Matrix<f64> = Matrix::<f64>::new(self.dim);
        for r in 0..self.dim.0 {
            for c in 0..self.dim.1 {
                result.set(r, c, self.get(r, c) + rhs.get(r, c));
            }
        }
        return result;
    }
}

impl std::ops::Sub for Matrix<f64> {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        assert_eq!(self.dim, rhs.dim, "Matrix dimensions must match for subtraction.");
        let mut result: Matrix<f64> = Matrix::<f64>::new(self.dim);
        for r in 0..self.dim.0 {
            for c in 0..self.dim.1 {
                result.set(r, c, self.get(r, c) - rhs.get(r, c));
            }
        }
        return result;
    }
}

impl std::ops::Mul for Matrix<f64> {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
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
        let mut m: Matrix<f64> = Matrix::<f64>::new( (2, 3));
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

        let m3: Matrix<f64> = m1 + m2;
        assert_eq!(m3.get(0, 0), 6.0);
        assert_eq!(m3.get(0, 1), 8.0);
        assert_eq!(m3.get(1, 0), 10.0);
        assert_eq!(m3.get(1, 1), 12.0);
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

        let m3: Matrix<f64> = m1 - m2;
        assert_eq!(m3.get(0, 0), 4.0);
        assert_eq!(m3.get(0, 1), 4.0);
        assert_eq!(m3.get(1, 0), 4.0);
        assert_eq!(m3.get(1, 1), 4.0);
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

        let m3: Matrix<f64> = m1 * m2;
        assert_eq!(m3.get(0, 0), 19.0);
        assert_eq!(m3.get(0, 1), 22.0);
        assert_eq!(m3.get(1, 0), 43.0);
        assert_eq!(m3.get(1, 1), 50.0);
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

        let m3: Matrix<f64> = m1 * m2;
        assert_eq!(m3.get(0, 0), 58.0);
        assert_eq!(m3.get(0, 1), 64.0);
        assert_eq!(m3.get(1, 0), 139.0);
        assert_eq!(m3.get(1, 1), 154.0);
    }
}
