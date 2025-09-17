mod basis_functions;
mod lin_alg;
mod meshing;
mod problem;

fn kelvin_to_celsius(k: f64) -> f64 {
    k - 273.15
}

fn main() {
    println!("Finite Element Method 1D Solver");

    let length: f64 = 0.1; // meters
    let n_elements: usize = 100;
    let thermal_conductivity: f64 = 45.0; // W/m.K for steel
    let area = 0.01; // m^2
    let left_bc: problem::BoundaryCondition = problem::BoundaryCondition::Dirichlet(373.15); // Temperature at left end in K
    let right_bc: problem::BoundaryCondition = problem::BoundaryCondition::Dirichlet(-1000.0); // Heat flux at right end in W/m^2

    println!("Simulating heat conduction in a {length} m rod with {n_elements} elements...");

    let mut mesh: meshing::Mesh1D = meshing::Mesh1D::uniform_mesh(length, n_elements);
    mesh.elasticity = thermal_conductivity;
    mesh.area = area;
    mesh.internal_force = -25.0; // W/m^3, volumetric heat generation/loss

    let mut problem: problem::Problem1D = problem::Problem1D::new(mesh, left_bc, right_bc);
    problem.solve();

    println!("Nodal temperatures (°C):");
    for i in 0..=n_elements {
        let temp_c: f64 = match problem.mesh.get_node_values(i) {
            Some(t) => kelvin_to_celsius(t),
            None => panic!("Node value not found"),
        };
        let coord = match problem.mesh.get_node_coordinates(i) {
            Some(c) => c,
            None => panic!("Node coordinate not found"),
        };
        println!("Node {i}: x = {coord:.4} m, T = {temp_c:.2} °C");
    }
}
