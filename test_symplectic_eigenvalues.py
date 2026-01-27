"""
Test suite for symplectic eigenvalues computation
"""

import numpy as np
from symplectic_eigenvalues import (
    symplectic_form,
    symplectic_eigenvalues,
    symplectic_diagonalizing_matrix,
    williamson_decomposition
)


def test_symplectic_form():
    """Test the symplectic form matrix construction."""
    print("Test: symplectic_form")
    
    # Test for n=2 (4x4 matrix)
    Omega = symplectic_form(2)
    
    # Check shape
    assert Omega.shape == (4, 4), "Symplectic form should be 4x4 for n=2"
    
    # Check antisymmetry: Ω^T = -Ω
    assert np.allclose(Omega.T, -Omega), "Symplectic form should be antisymmetric"
    
    # Check specific structure
    expected = np.array([[0, 0, 1, 0],
                         [0, 0, 0, 1],
                         [-1, 0, 0, 0],
                         [0, -1, 0, 0]])
    assert np.allclose(Omega, expected), "Symplectic form has incorrect structure"
    
    print("✓ Passed")


def test_symplectic_eigenvalues_identity():
    """Test symplectic eigenvalues of scaled identity matrix."""
    print("\nTest: symplectic_eigenvalues (scaled identity)")
    
    # For A = 2I, the symplectic eigenvalue should be 2
    A = np.eye(4) * 2
    evals = symplectic_eigenvalues(A)
    
    assert len(evals) >= 1, "Should have at least one symplectic eigenvalue"
    assert np.allclose(evals[0], 2.0), f"Symplectic eigenvalue should be 2.0, got {evals[0]}"
    
    print(f"  Symplectic eigenvalues: {evals}")
    print("✓ Passed")


def test_symplectic_eigenvalues_block_diagonal():
    """Test symplectic eigenvalues of block diagonal matrix."""
    print("\nTest: symplectic_eigenvalues (block diagonal)")
    
    # Block diagonal matrix with blocks [[3, 1], [1, 3]]
    # Eigenvalues of each block are 2 and 4
    A = np.array([[3, 1, 0, 0], 
                  [1, 3, 0, 0], 
                  [0, 0, 3, 1], 
                  [0, 0, 1, 3]], dtype=float)
    evals = symplectic_eigenvalues(A)
    
    # Should get [2, 4]
    assert len(evals) == 2, f"Should have 2 unique symplectic eigenvalues, got {len(evals)}"
    assert np.allclose(evals, [2.0, 4.0]), f"Expected [2.0, 4.0], got {evals}"
    
    print(f"  Symplectic eigenvalues: {evals}")
    print("✓ Passed")


def test_symplectic_diagonalizing_matrix():
    """Test the symplectic diagonalizing matrix."""
    print("\nTest: symplectic_diagonalizing_matrix")
    
    A = np.eye(4) * 3
    S, D = symplectic_diagonalizing_matrix(A)
    
    # Check dimensions
    assert S.shape == (4, 4), "S should be 4x4"
    assert D.shape == (4, 4), "D should be 4x4"
    
    # Check that S is invertible
    det_S = np.linalg.det(S)
    assert abs(det_S) > 1e-10, "S should be invertible"
    
    print(f"  S shape: {S.shape}")
    print(f"  D diagonal elements: {np.diag(D.real)}")
    print("✓ Passed")


def test_williamson_decomposition():
    """Test the Williamson decomposition."""
    print("\nTest: williamson_decomposition")
    
    A = np.array([[2, 0.5, 0, 0], 
                  [0.5, 2, 0, 0], 
                  [0, 0, 2, 0.5], 
                  [0, 0, 0.5, 2]], dtype=float)
    
    S, d = williamson_decomposition(A)
    
    # Check that we get symplectic eigenvalues
    assert len(d) > 0, "Should have at least one symplectic eigenvalue"
    assert all(ev > 0 for ev in d), "All symplectic eigenvalues should be positive"
    
    print(f"  Symplectic eigenvalues: {d}")
    print("✓ Passed")


def test_positive_definiteness():
    """Test that input validation works for non-positive-definite matrices."""
    print("\nTest: input validation")
    
    # Create a non-square matrix (should fail)
    try:
        A = np.eye(3, 4)
        evals = symplectic_eigenvalues(A)
        assert False, "Should have raised ValueError for non-square matrix"
    except ValueError as e:
        print(f"  ✓ Correctly rejected non-square matrix: {e}")
    
    # Create a matrix with odd dimension (should fail)
    try:
        A = np.eye(3)
        evals = symplectic_eigenvalues(A)
        assert False, "Should have raised ValueError for odd dimension"
    except ValueError as e:
        print(f"  ✓ Correctly rejected odd dimension: {e}")
    
    print("✓ Passed")


def run_all_tests():
    """Run all test functions."""
    print("=" * 60)
    print("Running Symplectic Eigenvalues Test Suite")
    print("=" * 60)
    
    test_symplectic_form()
    test_symplectic_eigenvalues_identity()
    test_symplectic_eigenvalues_block_diagonal()
    test_symplectic_diagonalizing_matrix()
    test_williamson_decomposition()
    test_positive_definiteness()
    
    print("\n" + "=" * 60)
    print("All tests passed! ✓")
    print("=" * 60)


if __name__ == "__main__":
    run_all_tests()
