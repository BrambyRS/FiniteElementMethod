pub mod mat;

pub fn lin_solve(A: &mat::Matrix<f64>, b: &mat::Matrix<f64>) -> mat::Matrix<f64> {
    assert_eq!(
        A.get_dim().0,
        b.get_dim().0,
        "Incompatible matrix dimensions"
    );
    assert_eq!(b.get_dim().1, 1, "b must be a column vector");

    // TODO: Add more safety checks

    let (r, c) = A.get_dim();

    let (Q, R) = qr_decomposition(A);
    let y: mat::Matrix<f64> = &Q.transpose() * b;
    let mut x: mat::Matrix<f64> = mat::Matrix::new((A.get_dim().1, 1));
    for i in (0..r).rev() {
        let mut sum: f64 = 0.0;
        for j in (i + 1)..c {
            sum += R.get(i, j) * x.get(j, 0);
        }
        x.set(i, 0, (y.get(i, 0) - sum) / R.get(i, i));
    }

    return x;
}

fn qr_decomposition(A: &mat::Matrix<f64>) -> (mat::Matrix<f64>, mat::Matrix<f64>) {
    assert!(
        A.get_dim().0 >= A.get_dim().1,
        "Matrix A of dimension m x n must have m >= n."
    );

    // TODO: Add check to make sure A is not singular

    let (r, c) = A.get_dim();

    let mut Q: mat::Matrix<f64> = mat::Matrix::identity(r);
    let mut R: mat::Matrix<f64> = A.clone();

    // TODO: This can be optimised to avoid creating so many temporary matrices
    for i in 0..c {
        // Get the column vector from R for th k-th column from row k to r
        let rv: mat::Matrix<f64> = {
            let mut rv_temp = mat::Matrix::new((r - i, 1));
            for j in i..r {
                rv_temp.set(j - i, 0, R.get(j, i));
            }
            rv_temp
        };

        let e: mat::Matrix<f64> = {
            let mut e_temp: mat::Matrix<f64> = mat::Matrix::new((r - i, 1));
            e_temp.set(0, 0, rv.norm());
            e_temp
        };

        let sign_rv1 = if rv.get(0, 0) >= 0.0 { 1.0 } else { -1.0 };
        let v: mat::Matrix<f64> = {
            let mut v_temp: mat::Matrix<f64> = mat::Matrix::new((r - i, 1));
            for j in 0..(r - i) {
                v_temp.set(j, 0, rv.get(j, 0) + sign_rv1 * e.get(j, 0));
            }
            v_temp
        };
        let v_dot: f64 = v.dot_product(&v);
        let sub_h: mat::Matrix<f64> =
            mat::Matrix::identity(r - i) - 2.0 * (1.0 / v_dot) * &(&v * &v.transpose());
        let mut h: mat::Matrix<f64> = mat::Matrix::identity(r);
        for j in i..r {
            for k in i..r {
                h.set(j, k, sub_h.get(j - i, k - i));
            }
        }

        Q = &Q * &h.transpose();
        R = &h * &R;
    }

    return (Q, R);
}

#[cfg(test)]
mod tests {
    use super::*;

    // Helper function to check if all elements of a matrix are below a tolerance
    fn matrix_elements_below_tolerance(matrix: &mat::Matrix<f64>, tolerance: f64) -> bool {
        let (rows, cols) = matrix.get_dim();
        for r in 0..rows {
            for c in 0..cols {
                if matrix.get(r, c) >= tolerance {
                    return false;
                }
            }
        }
        true
    }

    #[test]
    fn test_qr_decomposition() {
        let mut a = mat::Matrix::new((3, 3));
        a.set(0, 0, 2.0);
        a.set(0, 1, -1.0);
        a.set(0, 2, 0.0);
        a.set(1, 0, 1.0);
        a.set(1, 1, 2.0);
        a.set(1, 2, 1.0);
        a.set(2, 0, 0.0);
        a.set(2, 1, 1.0);
        a.set(2, 2, 3.0);

        let (q, r) = qr_decomposition(&a);
        let a_reconstructed: mat::Matrix<f64> = &q * &r;

        // I also happen to know what Q and R shoud be for this example
        let mut q_expected = mat::Matrix::new((3, 3));
        q_expected.set(0, 0, -0.894427190999916);
        q_expected.set(0, 1, 0.408248290463863);
        q_expected.set(0, 2, -0.182574185835055);
        q_expected.set(1, 0, -0.447213595499958);
        q_expected.set(1, 1, -0.816496580927726);
        q_expected.set(1, 2, 0.365148371670111);
        q_expected.set(2, 0, 0.000000000000000);
        q_expected.set(2, 1, -0.408248290463863);
        q_expected.set(2, 2, -0.912870929175277);

        let mut r_expected = mat::Matrix::new((3, 3));
        r_expected.set(0, 0, -2.236067977499790);
        r_expected.set(0, 1, 0.000000000000000);
        r_expected.set(0, 2, -0.447213595499958);
        r_expected.set(1, 0, 0.000000000000000);
        r_expected.set(1, 1, -2.449489742783178);
        r_expected.set(1, 2, -2.041241452319315);
        r_expected.set(2, 0, 0.000000000000000);
        r_expected.set(2, 1, 0.000000000000000);
        r_expected.set(2, 2, -2.373464415855720);

        // Check that Q*R reconstructs the original matrix A
        // For floating point results, we need tolerance-based comparison
        let reconstruction_error = (&a - &a_reconstructed).abs();
        assert!(matrix_elements_below_tolerance(
            &reconstruction_error,
            1e-10
        ));

        // Note: Due to floating point precision, we can't use exact equality for Q and R
        // so we use matrix operations with tolerance checking
        let q_error = (&q - &q_expected).abs();
        let r_error = (&r - &r_expected).abs();
        assert!(matrix_elements_below_tolerance(&q_error, 1e-10));
        assert!(matrix_elements_below_tolerance(&r_error, 1e-10));
    }

    #[test]
    fn test_solve() {
        let mut a = mat::Matrix::new((3, 3));
        a.set(0, 0, 2.0);
        a.set(0, 1, -1.0);
        a.set(0, 2, 0.0);
        a.set(1, 0, 1.0);
        a.set(1, 1, 2.0);
        a.set(1, 2, 1.0);
        a.set(2, 0, 0.0);
        a.set(2, 1, 1.0);
        a.set(2, 2, 3.0);

        let mut b = mat::Matrix::new((3, 1));
        b.set(0, 0, 1.0);
        b.set(1, 0, 4.0);
        b.set(2, 0, 2.0);

        let x = lin_solve(&a, &b);

        let mut x_expected = mat::Matrix::new((3, 1));
        x_expected.set(0, 0, 1.153846153846154);
        x_expected.set(1, 0, 1.307692307692307);
        x_expected.set(2, 0, 0.230769230769231);

        // Note: Due to floating point precision, we can't use exact equality
        // so we use matrix operations with tolerance checking
        let solution_error = (&x - &x_expected).abs();
        assert!(matrix_elements_below_tolerance(&solution_error, 1e-10));
    }
}
