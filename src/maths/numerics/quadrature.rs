pub struct GaussLegendre {
    pub order: usize,
    pub nodes: Vec<f64>,
    pub weights: Vec<f64>,
}

impl GaussLegendre {
    pub fn new(order: usize) -> Self {
        // Get nodes and weights in bi-unit domain [-1, 1]
        // TODO: Generalise for higher orders
        let (nodes, weights) = match order {
            1 => (vec![0.0], vec![2.0]),
            2 => {
                let term: f64 = 1.0 / 3f64.sqrt();
                (vec![-term, term], vec![1.0, 1.0])
            }
            3 => {
                let term: f64 = (3f64 / 5f64).sqrt();
                (
                    vec![-term, 0.0, term],
                    vec![5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0],
                )
            }
            _ => panic!("Order not implemented"),
        };

        return Self {
            order,
            nodes,
            weights,
        };
    }

    pub fn integrate(&self, f: fn(f64) -> f64, a: f64, b: f64) -> f64 {
        assert!(
            a < b,
            "Invalid integration limits. a = {a} must be less than b = {b}."
        );

        let mut integral: f64 = 0.0;
        for i in 0..self.order {
            let x: f64 = ((b - a) * self.nodes[i] + (a + b)) * 0.5;
            integral += self.weights[i] * f(x);
        }
        return integral * (b - a) * 0.5;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_gauss_legendre_order_1() {
        let gl = GaussLegendre::new(1);
        assert_eq!(gl.nodes, vec![0.0]);
        assert_eq!(gl.weights, vec![2.0]);

        // Test integration of f(x) = 2x+1 over [2, 4]
        let expected_integral = 14.0;
        let integral = gl.integrate(|x: f64| 2.0 * x + 1.0, 2.0, 4.0);
        assert!(
            (integral - expected_integral).abs() < 1e-12,
            "Integral does not match expected value. Got {integral}, expected {expected_integral}"
        );
    }

    #[test]
    fn test_gauss_legendre_order_2() {
        let gl = GaussLegendre::new(2);
        let term = 1.0 / 3f64.sqrt();
        let expected_nodes = vec![-term, term];
        assert!(
            gl.nodes
                .iter()
                .zip(expected_nodes.iter())
                .all(|(a, b)| (a - b).abs() < 1e-12),
            "Nodes do not match expected values"
        );
        assert_eq!(gl.weights, vec![1.0, 1.0]);

        // Test integration of f(x) = x^3-3x^2+3x-1 over [-1, 2]
        let expected_integral: f64 = -3.75;
        let integral = gl.integrate(
            |x: f64| x.powi(3) - 3.0 * x.powi(2) + 3.0 * x - 1.0,
            -1.0,
            2.0,
        );
        assert!(
            (integral - expected_integral).abs() < 1e-12,
            "Integral does not match expected value. Got {integral}, expected {expected_integral}",
        );
    }

    #[test]
    fn test_gauss_legendre_order_3() {
        let gl = GaussLegendre::new(3);
        let term: f64 = (3f64 / 5f64).sqrt();
        let expected_nodes = vec![-term, 0.0, term];
        assert!(
            gl.nodes
                .iter()
                .zip(expected_nodes.iter())
                .all(|(a, b)| (a - b).abs() < 1e-12),
            "Nodes do not match expected values"
        );
        assert_eq!(gl.weights, vec![5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0]);

        // Test integration of f(x) = 6x^5-3x^3-3x^2+2 over [-3, 2]
        let expected_integral: f64 = -641.25;
        let integral = gl.integrate(
            |x: f64| 6.0 * x.powi(5) - 3.0 * x.powi(3) - 3.0 * x.powi(2) + 2.0,
            -3.0,
            2.0,
        );
        assert!(
            (integral - expected_integral).abs() < 1e-12,
            "Integral does not match expected value. Got {integral}, expected {expected_integral}",
        );
    }

    #[test]
    fn test_faulty_interval() {
        let gl = GaussLegendre::new(2);
        let result: Result<(), Box<dyn std::any::Any + Send + 'static>> =
            std::panic::catch_unwind(|| {
                gl.integrate(|x: f64| x.powi(2), 2.0, 1.0);
            });
        assert!(result.is_err(), "Expected panic for invalid interval");
    }
}
