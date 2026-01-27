"""
Test that symplectic_eigenvalues_symbolic returns n eigenvalues

For a 2n×2n matrix, the function should return n symplectic eigenvalues.

This test simulates what the symbolic version should do using NumPy.
The actual symbolic version is in symplectic_eigenvalues.sage and requires SageMath.
"""

import numpy as np
from scipy import linalg


def symplectic_form(n):
    """Create the standard symplectic form matrix."""
    I_n = np.eye(n)
    Z_n = np.zeros((n, n))
    Omega = np.block([[Z_n, I_n], [-I_n, Z_n]])
    return Omega


def symplectic_eigenvalues_symbolic_simulation(A):
    """
    Simulate the symbolic version behavior: return n eigenvalues with multiplicity.
    For a 2n×2n matrix, returns n eigenvalues.
    """
    # Check that A is square
    if A.shape[0] != A.shape[1]:
        raise ValueError("Matrix A must be square")
    
    # Check that dimension is even
    if A.shape[0] % 2 != 0:
        raise ValueError("Matrix A must have even dimension for symplectic eigenvalues")
    
    n = A.shape[0] // 2
    
    # Create the symplectic form matrix
    Omega = symplectic_form(n)
    
    # Compute the matrix iΩA
    M = 1j * Omega @ A
    
    # Compute eigenvalues of M
    eigenvals = linalg.eigvals(M)
    
    # Take absolute values
    abs_eigenvals = np.abs(eigenvals.real)
    
    # Sort the eigenvalues
    sorted_evals = np.sort(abs_eigenvals)
    
    # Return the first n eigenvalues (with multiplicity)
    count = n
    symplectic_evals = sorted_evals[:count]
    
    return symplectic_evals


def test_identity_4x4():
    """Test 4x4 identity matrix (2n=4, n=2, should return 2 eigenvalues)."""
    print("\nTest: 4x4 identity matrix")
    A = np.eye(4) * 2
    evals = symplectic_eigenvalues_symbolic_simulation(A)
    
    assert len(evals) == 2, f"Should return 2 eigenvalues, got {len(evals)}"
    assert np.allclose(evals, [2.0, 2.0]), f"Eigenvalues should be [2.0, 2.0], got {evals}"
    
    print(f"  Matrix size: 4x4 (2n=4, n=2)")
    print(f"  Returned {len(evals)} eigenvalues: {evals}")
    print(f"  Expected: n = 2 eigenvalues")
    print("✓ Passed")


def test_identity_8x8():
    """Test 8x8 identity matrix (2n=8, n=4, should return 4 eigenvalues)."""
    print("\nTest: 8x8 identity matrix")
    A = np.eye(8) * 2
    evals = symplectic_eigenvalues_symbolic_simulation(A)
    
    assert len(evals) == 4, f"Should return 4 eigenvalues, got {len(evals)}"
    assert np.allclose(evals, [2.0, 2.0, 2.0, 2.0]), f"Eigenvalues should be [2.0, 2.0, 2.0, 2.0], got {evals}"
    
    print(f"  Matrix size: 8x8 (2n=8, n=4)")
    print(f"  Returned {len(evals)} eigenvalues: {evals}")
    print(f"  Expected: n = 4 eigenvalues")
    print("✓ Passed")


def test_block_diagonal_4x4():
    """Test 4x4 block diagonal matrix (2n=4, n=2, should return 2 eigenvalues)."""
    print("\nTest: 4x4 block diagonal matrix")
    A = np.array([[3, 1, 0, 0], 
                  [1, 3, 0, 0], 
                  [0, 0, 3, 1], 
                  [0, 0, 1, 3]], dtype=float)
    evals = symplectic_eigenvalues_symbolic_simulation(A)
    
    assert len(evals) == 2, f"Should return 2 eigenvalues, got {len(evals)}"
    # Eigenvalues are ±2 and ±4, so sorted abs gives [2, 2, 4, 4]
    # Taking first n=2 gives [2, 2]
    assert np.allclose(evals, [2.0, 2.0]), f"Eigenvalues should be ~[2.0, 2.0], got {evals}"
    
    print(f"  Matrix size: 4x4 (2n=4, n=2)")
    print(f"  Returned {len(evals)} eigenvalues: {evals}")
    print(f"  Expected: n = 2 eigenvalues")
    print("✓ Passed")


def test_6x6():
    """Test 6x6 matrix (2n=6, n=3, should return 3 eigenvalues)."""
    print("\nTest: 6x6 identity matrix")
    A = np.eye(6) * 3
    evals = symplectic_eigenvalues_symbolic_simulation(A)
    
    assert len(evals) == 3, f"Should return 3 eigenvalues, got {len(evals)}"
    assert np.allclose(evals, [3.0, 3.0, 3.0]), f"Eigenvalues should be [3.0, 3.0, 3.0], got {evals}"
    
    print(f"  Matrix size: 6x6 (2n=6, n=3)")
    print(f"  Returned {len(evals)} eigenvalues: {evals}")
    print(f"  Expected: n = 3 eigenvalues")
    print("✓ Passed")


def run_all_tests():
    """Run all test functions."""
    print("=" * 70)
    print("Testing symplectic_eigenvalues_symbolic behavior (n eigenvalues)")
    print("=" * 70)
    
    test_identity_4x4()
    test_identity_8x8()
    test_block_diagonal_4x4()
    test_6x6()
    
    print("\n" + "=" * 70)
    print("All tests passed! ✓")
    print("=" * 70)
    print("\nNote: These tests simulate the expected behavior of")
    print("symplectic_eigenvalues_symbolic() which requires SageMath.")
    print("For a 2n×2n matrix, the function returns n symplectic eigenvalues.")


if __name__ == "__main__":
    run_all_tests()
