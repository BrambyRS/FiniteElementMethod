import numpy as np

def problem_assembly_prototype(E, A, u0, f, t, h):
    """
    Python equivalent of the MATLAB ProblemAssemblyPrototype function.

    Assembles the global stiffness matrix K and force vector F for a 1D finite element problem.

    Parameters:
    -----------
    E : float
        Young's modulus (elasticity)
    A : float
        Cross-sectional area
    u0 : float
        Dirichlet boundary condition value at left boundary
    f : float
        Distributed force per unit volume
    t : float
        Neumann boundary condition value at right boundary (traction)
    h : array-like
        Element sizes (lengths)

    Returns:
    --------
    K : numpy.ndarray
        Global stiffness matrix (n_el x n_el)
    F : numpy.ndarray
        Global force vector (n_el x 1)
    """
    n_el = len(h)
    K = np.zeros((n_el, n_el))
    F = np.zeros((n_el, 1))

    # Shape function integrals and derivatives
    iN1 = 1.0  # Integral of N1 over the xi range [-1, 1]
    iN2 = 1.0  # Integral of N2 over the xi range [-1, 1]
    dN1 = -0.5  # Derivative of N1 with respect to xi
    dN2 = 0.5   # Derivative of N2 with respect to xi

    # Loop through elements for assembly
    for i in range(n_el - 1):
        # Local stiffness matrix
        K_local = 4 * E * A / h[i] * np.array([
            [dN1 * dN1, dN1 * dN2],
            [dN2 * dN1, dN2 * dN2]
        ])

        # Local force vector
        F_local = f * A * h[i] * np.array([[iN1], [iN2]])

        # Assembly into global matrices
        K[i:i+2, i:i+2] += K_local
        F[i:i+2] += F_local

    # Apply Dirichlet boundary condition at left boundary
    K[0, 0] += 4 * E * A / h[0] * dN2 * dN2
    F[0] += h[0] - 4 * E * A / h[0] * dN1 * dN2 * u0

    # Apply Neumann boundary condition at right boundary
    F[-1] += t * A

    return K, F


def main():
    """
    Example usage and test of the problem_assembly_prototype function.
    """
    # Example parameters
    E = 1.0      # Young's modulus
    A = 1.0      # Cross-sectional area
    u0 = 273.15  # Dirichlet BC (temperature at left boundary)
    f = 1.0      # Distributed force
    t = 100.0    # Neumann BC (heat flux at right boundary)
    h = [1.0, 1.0, 1.0, 1.0, 1.0]  # Element sizes (5 elements of size 1)

    # Assemble matrices
    K, F = problem_assembly_prototype(E, A, u0, f, t, h)

    print("Global Stiffness Matrix K:")
    print(K)
    print("\nGlobal Force Vector F:")
    print(F)

    # Verify the expected structure for this example
    print("\nExpected K matrix structure (for verification):")
    expected_K = np.array([
        [2.0, -1.0, 0.0, 0.0, 0.0],
        [-1.0, 2.0, -1.0, 0.0, 0.0],
        [0.0, -1.0, 2.0, -1.0, 0.0],
        [0.0, 0.0, -1.0, 2.0, -1.0],
        [0.0, 0.0, 0.0, -1.0, 1.0]
    ])
    print(expected_K)

    print("\nExpected F vector structure (for verification):")
    expected_F = np.array([[273.15 + 2.0], [2.0], [2.0], [2.0], [1.0 + 100.0]])
    print(expected_F)


if __name__ == "__main__":
    main()
