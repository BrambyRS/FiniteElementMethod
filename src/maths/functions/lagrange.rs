pub fn lagrange(x: f64, order: usize, n: usize) -> f64 {
    if n > order + 1 {
        panic!("Lagrange polynomial n = {n} does not exist for order = {order}. Value n must be less than order.");
    }

    let mut nodes: Vec<f64> = Vec::with_capacity(order + 1);
    for i in 0..order + 1 {
        let node: f64 = 2.0 * (i as f64 / order as f64) - 1.0;
        nodes.push(node);
    }

    let mut numerator: f64 = 1.0;
    let mut denominator: f64 = 1.0;

    for i in 0..order + 1 {
        if i + 1 != n {
            numerator *= x - nodes[i];
            denominator *= nodes[n - 1] - nodes[i];
        }
    }

    println!("x: {x}, order: {order}, n: {n}");
    println!("Numerator: {numerator}, Denominator: {denominator}");
    return numerator / denominator;
}

#[cfg(test)]
mod tests {
    use super::*;

    fn test_kronecker_delta_property(order: usize) {
        let mut nodes: Vec<f64> = Vec::with_capacity(order + 1);
        for i in 0..order + 1 {
            let node: f64 = 2.0 * (i as f64 / order as f64) - 1.0;
            nodes.push(node);
        }

        for i in 0..order + 1 {
            for j in 0..order + 1 {
                let value = lagrange(nodes[i], order, j + 1);
                if i == j {
                    assert!(
                        (value - 1.0).abs() < 1e-12,
                        "Failed Kronecker delta property at node {}: expected 1, got {}",
                        i + 1,
                        value
                    );
                } else {
                    assert!(
                        value.abs() < 1e-12,
                        "Failed Kronecker delta property at node {}: expected 0, got {}",
                        i + 1,
                        value
                    );
                }
            }
        }
    }

    #[test]
    fn test_lagrange_linear() {
        let x: [f64; 3] = [-1.0, 0.0, 1.0];
        let order: usize = 1;
        let n: [usize; 2] = [1, 2];
        let expected_n1: [f64; 3] = [1.0, 0.5, 0.0];
        let expected_n2: [f64; 3] = [0.0, 0.5, 1.0];
        for i in 0..x.len() {
            let computed_n1 = lagrange(x[i], order, n[0]);
            let computed_n2 = lagrange(x[i], order, n[1]);
            assert!(
                (computed_n1 - expected_n1[i]).abs() < 1e-12,
                "Failed at x = {}. Expected {}, got {}",
                x[i],
                expected_n1[i],
                computed_n1
            );
            assert!(
                (computed_n2 - expected_n2[i]).abs() < 1e-12,
                "Failed at x = {}. Expected {}, got {}",
                x[i],
                expected_n2[i],
                computed_n2
            );
        }
        test_kronecker_delta_property(order);
    }
}
