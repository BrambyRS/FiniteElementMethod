/*
Only allows for problems with a Dirichlet boundary condition at the left boundary
and a Neumann boundary condition at the right boundary
*/

use crate::lin_alg;
use crate::lin_alg::mat;
use crate::meshing::Mesh1D;

pub enum BoundaryCondition {
    Dirichlet(f64),
    Neumann(f64),
}

pub struct Problem1D {
    pub mesh: Mesh1D,
    pub left_bc: BoundaryCondition,
    pub right_bc: BoundaryCondition,
}

impl Problem1D {
    pub fn new(mesh: Mesh1D, left_bc: BoundaryCondition, right_bc: BoundaryCondition) -> Self {
        return Self {
            mesh,
            left_bc,
            right_bc,
        };
    }

    pub fn solve(&mut self) {
        let (k, f) = self.construct_matrices();
        let d: mat::Matrix<f64> = lin_alg::lin_solve(&k, &f);

        // Apply nodal values to mesh
        // Length of d depends on the boundary conditions
        let n_el = self.mesh.get_n_elements();
        let mut node_idx: usize = 0;

        match self.left_bc {
            BoundaryCondition::Dirichlet(val) => {
                self.mesh.node_values[0] = val;
                node_idx += 1;
            }
            BoundaryCondition::Neumann(_) => {}
        }

        for i in 0..d.get_dim().0 {
            self.mesh.node_values[i + node_idx] = d.get(i, 0);
        }

        match self.right_bc {
            BoundaryCondition::Dirichlet(val) => {
                self.mesh.node_values[n_el] = val;
            }
            BoundaryCondition::Neumann(_) => {}
        }
    }

    fn construct_matrices(&self) -> (mat::Matrix<f64>, mat::Matrix<f64>) {
        let int_n1: f64 = 1.0;
        let int_n2: f64 = 1.0;
        let int_dn1: f64 = -0.5;
        let int_dn2: f64 = 0.5;

        let n_el: usize = self.mesh.get_n_elements();
        // Matrix sizes depend on the degrees of freedom, where Dirichlet BC reduces by 1
        let mut n_dof: usize = n_el + 1;
        match self.left_bc {
            BoundaryCondition::Dirichlet(_) => n_dof -= 1,
            _ => {}
        }
        match self.right_bc {
            BoundaryCondition::Dirichlet(_) => n_dof -= 1,
            _ => {}
        }

        let mut k: mat::Matrix<f64> = mat::Matrix::new((n_dof, n_dof));
        let mut f: mat::Matrix<f64> = mat::Matrix::new((n_dof, 1));

        let k_fac: f64 = 4.0 * self.mesh.elasticity * self.mesh.area;
        let f_fac: f64 = self.mesh.internal_force * self.mesh.area;

        // Apply left boundary condition
        let k_fac_local: f64 = match self.mesh.get_element_size(0) {
            Some(size) => k_fac / size,
            None => panic!("Element size not found"),
        };

        let f_fac_local: f64 = match self.mesh.get_element_size(0) {
            Some(size) => f_fac * size,
            None => panic!("Element size not found"),
        };

        match self.left_bc {
            BoundaryCondition::Dirichlet(val) => {
                k.set(0, 0, k.get(0, 0) + k_fac_local * int_dn2 * int_dn2);
                f.set(
                    0,
                    0,
                    f.get(0, 0) + f_fac_local * int_n2 - k_fac_local * int_dn1 * int_dn2 * val,
                );
            }
            BoundaryCondition::Neumann(val) => {
                f.set(0, 0, f.get(0, 0) + val * self.mesh.area);
            }
        }

        for i in 0..(n_dof - 1) {
            let mut k_local: mat::Matrix<f64> = mat::Matrix::new((2, 2));
            let mut f_local: mat::Matrix<f64> = mat::Matrix::new((2, 1));

            let k_fac_element: f64 = match self.mesh.get_element_size(i) {
                Some(size) => k_fac / size,
                None => panic!("Element size not found"),
            };

            let f_fac_element: f64 = match self.mesh.get_element_size(i) {
                Some(size) => f_fac * size,
                None => panic!("Element size not found"),
            };

            k_local.set(0, 0, k_fac_element * int_dn1 * int_dn1);
            k_local.set(0, 1, k_fac_element * int_dn1 * int_dn2);
            k_local.set(1, 0, k_fac_element * int_dn2 * int_dn1);
            k_local.set(1, 1, k_fac_element * int_dn2 * int_dn2);

            f_local.set(0, 0, f_fac_element * int_n1);
            f_local.set(1, 0, f_fac_element * int_n2);

            // Assemble into global matrices
            k.set(i, i, k.get(i, i) + k_local.get(0, 0));
            k.set(i, i + 1, k.get(i, i + 1) + k_local.get(0, 1));
            k.set(i + 1, i, k.get(i + 1, i) + k_local.get(1, 0));
            k.set(i + 1, i + 1, k.get(i + 1, i + 1) + k_local.get(1, 1));

            f.set(i, 0, f.get(i, 0) + f_local.get(0, 0));
            f.set(i + 1, 0, f.get(i + 1, 0) + f_local.get(1, 0));
        }

        // Adjust for right boundary condition
        let k_fac_local: f64 = match self.mesh.get_element_size(n_el - 1) {
            Some(size) => k_fac / size,
            None => panic!("Element size not found"),
        };

        let f_fac_local: f64 = match self.mesh.get_element_size(n_el - 1) {
            Some(size) => f_fac * size,
            None => panic!("Element size not found"),
        };

        match self.right_bc {
            BoundaryCondition::Dirichlet(val) => {
                k.set(
                    n_dof - 1,
                    n_dof - 1,
                    k.get(n_dof - 1, n_dof - 1) + k_fac_local * int_dn1 * int_dn1,
                );
                f.set(
                    n_dof - 1,
                    0,
                    f.get(n_dof - 1, 0) + f_fac_local * int_n1
                        - k_fac_local * int_dn1 * int_dn2 * val,
                );
            }
            BoundaryCondition::Neumann(val) => {
                f.set(n_dof - 1, 0, f.get(n_dof - 1, 0) + val * self.mesh.area);
            }
        }

        return (k, f);
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::meshing;

    #[test]
    fn test_bc_enum() {
        let dirichlet = BoundaryCondition::Dirichlet(100.0);
        let neumann = BoundaryCondition::Neumann(50.0);

        match dirichlet {
            BoundaryCondition::Dirichlet(val) => assert_eq!(val, 100.0),
            _ => panic!("Expected Dirichlet condition"),
        }

        match neumann {
            BoundaryCondition::Neumann(val) => assert_eq!(val, 50.0),
            _ => panic!("Expected Neumann condition"),
        }
    }

    #[test]
    fn test_problem_construction_dirichlet_neumann() {
        let mut mesh: meshing::Mesh1D = meshing::Mesh1D::uniform_mesh(5.0, 5); // To ensure elements are of size 1

        // Set all material properties to 1 for simplicity
        mesh.elasticity = 1.0;
        mesh.area = 1.0;
        mesh.internal_force = 1.0;

        let left_bc = BoundaryCondition::Dirichlet(273.15);
        let right_bc = BoundaryCondition::Neumann(100.0);

        let problem: Problem1D = Problem1D::new(mesh, left_bc, right_bc);
        // assert_eq!(problem.dirichlet_condition, 273.15);
        // assert_eq!(problem.neumann_condition, 100.0);
        assert_eq!(problem.mesh.get_n_elements(), 5);

        let mut k_expected: mat::Matrix<f64> = mat::Matrix::new((5, 5));
        k_expected.set(0, 0, 2.0);
        k_expected.set(0, 1, -1.0);
        k_expected.set(1, 0, -1.0);
        k_expected.set(1, 1, 2.0);
        k_expected.set(1, 2, -1.0);
        k_expected.set(2, 1, -1.0);
        k_expected.set(2, 2, 2.0);
        k_expected.set(2, 3, -1.0);
        k_expected.set(3, 2, -1.0);
        k_expected.set(3, 3, 2.0);
        k_expected.set(3, 4, -1.0);
        k_expected.set(4, 3, -1.0);
        k_expected.set(4, 4, 1.0);

        let mut f_expected: mat::Matrix<f64> = mat::Matrix::new((5, 1));
        f_expected.set(0, 0, 273.15 + 2.0);
        f_expected.set(1, 0, 2.0);
        f_expected.set(2, 0, 2.0);
        f_expected.set(3, 0, 2.0);
        f_expected.set(4, 0, 1.0 + 100.0);

        let (k, f) = problem.construct_matrices();
        assert_eq!(k, k_expected);
        assert_eq!(f, f_expected);
    }

    #[test]
    fn test_problem_construction_neumann_dirichlet() {
        let mut mesh: meshing::Mesh1D = meshing::Mesh1D::uniform_mesh(5.0, 5); // To ensure elements are of size 1

        // Set all material properties to 1 for simplicity
        mesh.elasticity = 1.0;
        mesh.area = 1.0;
        mesh.internal_force = 1.0;

        let left_bc = BoundaryCondition::Neumann(50.0);
        let right_bc = BoundaryCondition::Dirichlet(373.15);

        let problem: Problem1D = Problem1D::new(mesh, left_bc, right_bc);
        // assert_eq!(problem.dirichlet_condition, 273.15);
        // assert_eq!(problem.neumann_condition, 100.0);
        assert_eq!(problem.mesh.get_n_elements(), 5);

        let mut k_expected: mat::Matrix<f64> = mat::Matrix::new((5, 5));
        k_expected.set(0, 0, 1.0);
        k_expected.set(0, 1, -1.0);
        k_expected.set(1, 0, -1.0);
        k_expected.set(1, 1, 2.0);
        k_expected.set(1, 2, -1.0);
        k_expected.set(2, 1, -1.0);
        k_expected.set(2, 2, 2.0);
        k_expected.set(2, 3, -1.0);
        k_expected.set(3, 2, -1.0);
        k_expected.set(3, 3, 2.0);
        k_expected.set(3, 4, -1.0);
        k_expected.set(4, 3, -1.0);
        k_expected.set(4, 4, 2.0);

        let mut f_expected: mat::Matrix<f64> = mat::Matrix::new((5, 1));
        f_expected.set(0, 0, 50.0 + 1.0);
        f_expected.set(1, 0, 2.0);
        f_expected.set(2, 0, 2.0);
        f_expected.set(3, 0, 2.0);
        f_expected.set(4, 0, 2.0 + 373.15);

        let (k, f) = problem.construct_matrices();
        assert_eq!(k, k_expected);
        assert_eq!(f, f_expected);
    }

    #[test]
    fn test_problem_construction_pure_dirichlet() {
        let mut mesh: meshing::Mesh1D = meshing::Mesh1D::uniform_mesh(5.0, 5); // To ensure elements are of size 1

        // Set all material properties to 1 for simplicity
        mesh.elasticity = 1.0;
        mesh.area = 1.0;
        mesh.internal_force = 1.0;

        let left_bc = BoundaryCondition::Dirichlet(273.15);
        let right_bc = BoundaryCondition::Dirichlet(373.15);

        let problem: Problem1D = Problem1D::new(mesh, left_bc, right_bc);
        // assert_eq!(problem.dirichlet_condition, 273.15);
        // assert_eq!(problem.neumann_condition, 100.0);
        assert_eq!(problem.mesh.get_n_elements(), 5);

        let mut k_expected: mat::Matrix<f64> = mat::Matrix::new((4, 4));
        k_expected.set(0, 0, 2.0);
        k_expected.set(0, 1, -1.0);
        k_expected.set(1, 0, -1.0);
        k_expected.set(1, 1, 2.0);
        k_expected.set(1, 2, -1.0);
        k_expected.set(2, 1, -1.0);
        k_expected.set(2, 2, 2.0);
        k_expected.set(2, 3, -1.0);
        k_expected.set(3, 2, -1.0);
        k_expected.set(3, 3, 2.0);

        let mut f_expected: mat::Matrix<f64> = mat::Matrix::new((4, 1));
        f_expected.set(0, 0, 273.15 + 2.0);
        f_expected.set(1, 0, 2.0);
        f_expected.set(2, 0, 2.0);
        f_expected.set(3, 0, 373.15 + 2.0);

        let (k, f) = problem.construct_matrices();
        assert_eq!(k, k_expected);
        assert_eq!(f, f_expected);
    }

    #[test]
    fn test_problem_construction_pure_neumann() {
        let mut mesh: meshing::Mesh1D = meshing::Mesh1D::uniform_mesh(5.0, 5); // To ensure elements are of size 1

        // Set all material properties to 1 for simplicity
        mesh.elasticity = 1.0;
        mesh.area = 1.0;
        mesh.internal_force = 1.0;

        let left_bc = BoundaryCondition::Neumann(50.0);
        let right_bc = BoundaryCondition::Neumann(100.0);

        let problem: Problem1D = Problem1D::new(mesh, left_bc, right_bc);
        // assert_eq!(problem.dirichlet_condition, 273.15);
        // assert_eq!(problem.neumann_condition, 100.0);
        assert_eq!(problem.mesh.get_n_elements(), 5);

        let mut k_expected: mat::Matrix<f64> = mat::Matrix::new((6, 6));
        k_expected.set(0, 0, 1.0);
        k_expected.set(0, 1, -1.0);
        k_expected.set(1, 0, -1.0);
        k_expected.set(1, 1, 2.0);
        k_expected.set(1, 2, -1.0);
        k_expected.set(2, 1, -1.0);
        k_expected.set(2, 2, 2.0);
        k_expected.set(2, 3, -1.0);
        k_expected.set(3, 2, -1.0);
        k_expected.set(3, 3, 2.0);
        k_expected.set(3, 4, -1.0);
        k_expected.set(4, 3, -1.0);
        k_expected.set(4, 4, 2.0);
        k_expected.set(4, 5, -1.0);
        k_expected.set(5, 4, -1.0);
        k_expected.set(5, 5, 1.0);

        let mut f_expected: mat::Matrix<f64> = mat::Matrix::new((6, 1));
        f_expected.set(0, 0, 50.0 + 1.0);
        f_expected.set(1, 0, 2.0);
        f_expected.set(2, 0, 2.0);
        f_expected.set(3, 0, 2.0);
        f_expected.set(4, 0, 2.0);
        f_expected.set(5, 0, 1.0 + 100.0);

        let (k, f) = problem.construct_matrices();
        assert_eq!(k, k_expected);
        assert_eq!(f, f_expected);
    }
}
