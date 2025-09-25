pub struct GaussLegendre {
    pub nodes: Vec<f64>,
    pub weights: Vec<f64>,
}

impl GaussLegendre {
    pub fn new(interval: (f64, f64), order: usize) -> Self {
        // Get nodes and weights in bi-unit domain [-1, 1]
        // TODO: Generalise for higher orders
        let (nodes_biunit, weights) = match order {
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

        // Map nodes to desired interval [a, b]
        let a: f64 = interval.0;
        let b: f64 = interval.1;
        let mut nodes: Vec<f64> = Vec::with_capacity(nodes_biunit.len());
        for &xi in nodes_biunit.iter() {
            let x_mapped: f64 = 0.5 * (b - a) * xi + 0.5 * (a + b);
            nodes.push(x_mapped);
        }

        return Self { nodes, weights };
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_gauss_legendre_order_1() {
        let gl = GaussLegendre::new((-2.0, 4.0), 1);
        assert_eq!(gl.nodes, vec![1.0]);
        assert_eq!(gl.weights, vec![2.0]);
    }

    #[test]
    fn test_gauss_legendre_order_2() {
        let gl = GaussLegendre::new((1.0, 3.0), 2);
        let term = 1.0 / 3f64.sqrt();
        let expected_nodes = vec![2.0 - term, 2.0 + term];
        assert!(
            gl.nodes
                .iter()
                .zip(expected_nodes.iter())
                .all(|(a, b)| (a - b).abs() < 1e-12),
            "Nodes do not match expected values"
        );
        assert_eq!(gl.weights, vec![1.0, 1.0]);
    }

    #[test]
    fn test_gauss_legendre_order_3() {
        let gl = GaussLegendre::new((-1.0, 5.0), 3);
        let term: f64 = (3f64 / 5f64).sqrt();
        let expected_nodes = vec![2.0 - 3.0 * term, 2.0, 2.0 + 3.0 * term];
        assert!(
            gl.nodes
                .iter()
                .zip(expected_nodes.iter())
                .all(|(a, b)| (a - b).abs() < 1e-12),
            "Nodes do not match expected values"
        );
        assert_eq!(gl.weights, vec![5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0]);
    }
}
