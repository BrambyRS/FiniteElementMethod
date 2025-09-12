use std::vec;

pub struct Mesh1D {
    n_elements: usize,

    element_size: Vec<f64>,     // Length of each element
    node_coordinates: Vec<f64>, // Position of each node

    pub node_values: Vec<f64>, // Value at each node

    pub elasticity: f64, // Material property (e.g., Young's modulus, thermal conductivity)
    pub area: f64,       // Cross-sectional area
}

impl Mesh1D {
    pub fn uniform_mesh(length: f64, n_elements: usize) -> Self {
        let element_size: Vec<f64> = vec![length / n_elements as f64; n_elements];
        let mut node_coordinates: Vec<f64> = Vec::with_capacity(n_elements + 1);
        let node_values: Vec<f64> = vec![0.0; n_elements + 1];
        for i in 0..=n_elements {
            node_coordinates.push(i as f64 * element_size[0]);
        }
        return Self {
            n_elements,
            element_size,
            node_coordinates,
            node_values,
            elasticity: 0.0,
            area: 0.0,
        };
    }

    pub fn general_mesh(node_coordinates: Vec<f64>) -> Self {
        let n_elements = node_coordinates.len() - 1;
        let mut element_size: Vec<f64> = Vec::with_capacity(n_elements);
        let node_values: Vec<f64> = vec![0.0; n_elements + 1];
        for i in 0..n_elements {
            element_size.push(node_coordinates[i + 1] - node_coordinates[i]);
        }
        return Self {
            n_elements,
            element_size,
            node_coordinates,
            node_values,
            elasticity: 0.0,
            area: 0.0,
        };
    }

    pub fn get_n_elements(&self) -> usize {
        return self.n_elements;
    }

    pub fn get_element_size(&self, element: usize) -> Option<f64> {
        self.element_size.get(element).copied()
    }

    pub fn get_node_coordinates(&self, node: usize) -> Option<f64> {
        self.node_coordinates.get(node).copied()
    }

    pub fn get_node_values(&self, node: usize) -> Option<f64> {
        self.node_values.get(node).copied()
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
        assert_eq!(mesh.node_values, vec![0.0; 6]);
    }

    #[test]
    fn test_general_mesh() {
        let mesh = Mesh1D::general_mesh(vec![0.0, 2.0, 4.0, 6.0, 8.0, 10.0]);
        assert_eq!(mesh.n_elements, 5);
        assert_eq!(mesh.element_size, vec![2.0, 2.0, 2.0, 2.0, 2.0]);
        assert_eq!(mesh.node_coordinates, vec![0.0, 2.0, 4.0, 6.0, 8.0, 10.0]);
        assert_eq!(mesh.node_values, vec![0.0; 6]);
    }

    #[test]
    fn test_getters() {
        let mesh = Mesh1D::uniform_mesh(10.0, 5);
        assert_eq!(mesh.get_n_elements(), 5);
        assert_eq!(mesh.get_element_size(0), Some(2.0));
        assert_eq!(mesh.get_node_coordinates(0), Some(0.0));
        assert_eq!(mesh.get_node_values(0), Some(0.0));
    }
}
