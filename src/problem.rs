/*
Only allows for problems with a Dirichlet boundary condition at the left boundary
and a Neumann boundary condition at the right boundary
*/

use crate::lin_alg;
use crate::lin_alg::mat;
use crate::meshing::Mesh1D;

pub struct Problem1D {
    mesh: Mesh1D,
    dirichlet_condition: f64,
    neumann_condition: f64,
}

impl Problem1D {
    pub fn new(mesh: Mesh1D, dirichlet_condition: f64, neumann_condition: f64) -> Self {
        return Self {
            mesh,
            dirichlet_condition,
            neumann_condition,
        };
    }

    pub fn solve(&mut self) {
        let (k, f) = self.construct_matrices();
        let d: mat::Matrix<f64> = lin_alg::lin_solve(&k, &f);

        // Apply nodal values to mesh
        for i in 0..self.mesh.get_n_elements() + 1 {
            self.mesh.node_values[i] = d.get(i, 0);
        }
    }

    fn construct_matrices(&self) -> (mat::Matrix<f64>, mat::Matrix<f64>) {
        // Stiffness matrix
        let mut k: mat::Matrix<f64> = mat::Matrix::new((
            self.mesh.get_n_elements() + 1,
            self.mesh.get_n_elements() + 1,
        ));
        // Force matrix
        let mut f: mat::Matrix<f64> = mat::Matrix::new((self.mesh.get_n_elements() + 1, 1));

        // TODO: Implement matrix assembly

        return (k, f);
    }
}
