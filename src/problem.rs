/*
Only allows for problems with a Dirichlet boundary condition at the left boundary
and a Neumann boundary condition at the right boundary
*/

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

    pub fn solve(&self) {
        // TODO: Implement solve
    }
}
