struct Mesh1D {
    n_elements: usize,

    element_size: Vec<f64>,     // Length of each element
    node_coordinates: Vec<f64>, // Position of each node
}

impl Mesh1D {
    fn uniform_mesh(length: f64, n_elements: usize) -> Self {
        let element_size: Vec<f64> = vec![length / n_elements as f64; n_elements];
        let mut node_coordinates: Vec<f64> = Vec::with_capacity(n_elements + 1);
        for i in 0..=n_elements {
            node_coordinates.push(i as f64 * element_size[0]);
        }
        return Self {
            n_elements,
            element_size,
            node_coordinates,
        };
    }

    fn general_mesh(node_coordinates: Vec<f64>) -> Self {
        let n_elements = node_coordinates.len() - 1;
        let mut element_size: Vec<f64> = Vec::with_capacity(n_elements);
        for i in 0..n_elements {
            element_size.push(node_coordinates[i + 1] - node_coordinates[i]);
        }
        return Self {
            n_elements,
            element_size,
            node_coordinates,
        };
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_uniform_mesh() {
        let mesh = Mesh1D::uniform_mesh(10.0, 5);
        assert_eq!(mesh.n_elements, 5);
        assert_eq!(mesh.element_size, vec![2.0, 2.0, 2.0, 2.0, 2.0]);
        assert_eq!(mesh.node_coordinates, vec![0.0, 2.0, 4.0, 6.0, 8.0, 10.0]);
    }

    #[test]
    fn test_general_mesh() {
        let mesh = Mesh1D::general_mesh(vec![0.0, 2.0, 4.0, 6.0, 8.0, 10.0]);
        assert_eq!(mesh.n_elements, 5);
        assert_eq!(mesh.element_size, vec![2.0, 2.0, 2.0, 2.0, 2.0]);
        assert_eq!(mesh.node_coordinates, vec![0.0, 2.0, 4.0, 6.0, 8.0, 10.0]);
    }
}
