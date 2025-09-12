#![allow(non_snake_case)]
#![allow(dead_code)]

//! Linear basis functions and their derivatives for 1D elements.
//!
//! These functions implement the standard linear Lagrange basis functions
//! on the reference element ξ ∈ [-1, 1].

/// First linear basis function N₁(ξ) = (1 - ξ)/2.
///
/// This function equals 1 at ξ = -1 and 0 at ξ = 1, providing
/// linear interpolation from the left node of a 1D element.
///
/// # Arguments
///
/// * `xi` - Natural coordinate in the range [-1, 1]
///
/// # Returns
///
/// Value of the basis function at the given point
///
/// # Examples
///
/// ```
/// use finite_element_method::basis_functions::linear::N1;
/// assert_eq!(N1(-1.0), 1.0);  // Left node
/// assert_eq!(N1(1.0), 0.0);   // Right node
/// ```
pub fn N1(xi: f64) -> f64 {
    return (1.0 - xi) / 2.0;
}

/// Derivative of the first linear basis function dN₁/dξ = -1/2.
///
/// The derivative is constant since N₁ is linear in ξ.
///
/// # Arguments
///
/// * `_xi` - Natural coordinate (unused, kept for API consistency)
///
/// # Returns
///
/// Constant derivative value -0.5
pub fn dN1(_xi: f64) -> f64 {
    // xi argument is just for consistency
    return -0.5;
}

/// Second linear basis function N₂(ξ) = (1 + ξ)/2.
///
/// This function equals 0 at ξ = -1 and 1 at ξ = 1, providing
/// linear interpolation from the right node of a 1D element.
///
/// # Arguments
///
/// * `xi` - Natural coordinate in the range [-1, 1]
///
/// # Returns
///
/// Value of the basis function at the given point
///
/// # Examples
///
/// ```
/// use finite_element_method::basis_functions::linear::N2;
/// assert_eq!(N2(-1.0), 0.0);  // Left node
/// assert_eq!(N2(1.0), 1.0);   // Right node
pub fn N2(xi: f64) -> f64 {
    return (1.0 + xi) / 2.0;
}

/// Derivative of the second linear basis function dN₂/dξ = 1/2.
///
/// The derivative is constant since N₂ is linear in ξ.
///
/// # Arguments
///
/// * `_xi` - Natural coordinate (unused, kept for API consistency)
///
/// # Returns
///
/// Constant derivative value 0.5
pub fn dN2(_xi: f64) -> f64 {
    // xi argument is just for consistency
    return 0.5;
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_N1() {
        assert_eq!(N1(-1.0), 1.0);
        assert_eq!(N1(0.0), 0.5);
        assert_eq!(N1(1.0), 0.0);
    }

    #[test]
    fn test_dN1() {
        assert_eq!(dN1(-1.0), -0.5);
        assert_eq!(dN1(0.0), -0.5);
        assert_eq!(dN1(1.0), -0.5);
    }

    #[test]
    fn test_N2() {
        assert_eq!(N2(-1.0), 0.0);
        assert_eq!(N2(0.0), 0.5);
        assert_eq!(N2(1.0), 1.0);
    }

    #[test]
    fn test_dN2() {
        assert_eq!(dN2(-1.0), 0.5);
        assert_eq!(dN2(0.0), 0.5);
        assert_eq!(dN2(1.0), 0.5);
    }
}
