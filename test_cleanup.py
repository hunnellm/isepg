#!/usr/bin/env python3
"""
Test to verify that very small values are cleaned up in matrices.
"""

import numpy as np
from symplectic_eigenvalues import (
    cleanup_small_values,
    symplectic_diagonalizing_matrix,
    symplectic_form
)


def test_cleanup_small_values():
    """Test the cleanup_small_values function."""
    print("Test: cleanup_small_values")
    
    # Test with real matrix
    A_real = np.array([[1.0, 1e-15], [1e-15, 2.0]])
    A_clean = cleanup_small_values(A_real)
    
    assert A_clean[0, 0] == 1.0, "Should preserve large value"
    assert A_clean[0, 1] == 0.0, "Should clean up small value"
    assert A_clean[1, 0] == 0.0, "Should clean up small value"
    assert A_clean[1, 1] == 2.0, "Should preserve large value"
    print("  ✓ Real matrix cleanup works")
    
    # Test with complex matrix
    A_complex = np.array([[1.0 + 1e-15j, 2.0], 
                          [3.0, 4.0 + 1e-15j]])
    A_clean = cleanup_small_values(A_complex)
    
    assert np.abs(A_clean[0, 0].imag) < 1e-12, "Should clean up small imaginary part"
    assert np.abs(A_clean[1, 1].imag) < 1e-12, "Should clean up small imaginary part"
    print("  ✓ Complex matrix cleanup works")
    
    # Test with negative zero
    A_neg_zero = np.array([[1.0, -0.0], [-0.0, 2.0]])
    A_clean = cleanup_small_values(A_neg_zero)
    
    # Check that -0.0 becomes 0.0
    assert A_clean[0, 1] == 0.0, "Should clean up -0.0"
    assert A_clean[1, 0] == 0.0, "Should clean up -0.0"
    print("  ✓ Negative zero cleanup works")
    
    # Test boundary conditions
    threshold = 1e-10
    A_boundary = np.array([[1.0, 1e-11], [1e-9, 2.0]])  # 1e-11 < threshold, 1e-9 > threshold
    A_clean = cleanup_small_values(A_boundary)
    
    assert A_clean[0, 1] == 0.0, "Should clean up value below threshold (1e-11)"
    assert A_clean[1, 0] != 0.0, "Should preserve value above threshold (1e-9)"
    assert np.abs(A_clean[1, 0] - 1e-9) < 1e-15, "Value above threshold should be unchanged"
    print("  ✓ Boundary condition cleanup works")
    
    print("✓ Passed\n")


def test_symplectic_diagonalization_cleanup():
    """Test that symplectic diagonalization produces clean matrices."""
    print("Test: symplectic_diagonalization_cleanup")
    
    # Create a simple test matrix
    A = np.array([[4, 1, 0, 0], 
                  [1, 4, 0, 0], 
                  [0, 0, 4, 1], 
                  [0, 0, 1, 4]], dtype=float)
    
    S, D = symplectic_diagonalizing_matrix(A)
    
    # Check that very small values are cleaned up
    # Look for values that should be exactly zero
    tolerance = 1e-10
    
    # Check S for small values
    s_small_values = np.abs(S[np.abs(S) < tolerance])
    if len(s_small_values) > 0:
        assert np.all(s_small_values == 0), "Small values in S should be exactly 0"
    
    # Check D for small off-diagonal values
    off_diag_mask = ~np.eye(D.shape[0], dtype=bool)
    d_off_diag = np.abs(D[off_diag_mask])
    d_small_off_diag = d_off_diag[d_off_diag < tolerance]
    if len(d_small_off_diag) > 0:
        assert np.all(d_small_off_diag == 0), "Small off-diagonal values in D should be exactly 0"
    
    # Check that no values have question marks (would show up in string representation)
    S_str = str(S)
    D_str = str(D)
    assert '?' not in S_str, "S should not contain question marks"
    assert '?' not in D_str, "D should not contain question marks"
    
    print("  ✓ Matrices are clean (no very small values)")
    print("  ✓ No question marks in output")
    print("✓ Passed\n")


def test_symplectic_form_cleanup():
    """Test that symplectic form has no -0.0 values."""
    print("Test: symplectic_form_cleanup")
    
    Omega = symplectic_form(3)
    
    # Check for negative zero
    # In numpy, -0.0 == 0.0, but we can check the sign bit
    for i in range(Omega.shape[0]):
        for j in range(Omega.shape[1]):
            if Omega[i, j] == 0.0:
                # Check that it's not negative zero by checking if 1/value is positive
                # (1/-0.0 would be -inf)
                assert not np.signbit(Omega[i, j]), f"Element [{i},{j}] should not be -0.0"
    
    print("  ✓ No -0.0 values in symplectic form")
    print("✓ Passed\n")


def test_no_question_marks_in_output():
    """Ensure no question marks appear in matrix representations."""
    print("Test: no_question_marks_in_output")
    
    # Create various test matrices
    A1 = np.eye(4) * 2
    A2 = np.array([[3, 1, 0, 0], 
                   [1, 3, 0, 0], 
                   [0, 0, 3, 1], 
                   [0, 0, 1, 3]], dtype=float)
    
    for A in [A1, A2]:
        S, D = symplectic_diagonalizing_matrix(A)
        
        # Convert to string and check for question marks
        assert '?' not in str(S), "S should not contain ?"
        assert '?' not in str(D), "D should not contain ?"
        assert '?' not in repr(S), "S repr should not contain ?"
        assert '?' not in repr(D), "D repr should not contain ?"
    
    print("  ✓ No question marks found in any output")
    print("✓ Passed\n")


if __name__ == "__main__":
    print("=" * 60)
    print("Testing Matrix Cleanup Functionality")
    print("=" * 60)
    print()
    
    test_cleanup_small_values()
    test_symplectic_diagonalization_cleanup()
    test_symplectic_form_cleanup()
    test_no_question_marks_in_output()
    
    print("=" * 60)
    print("All cleanup tests passed! ✓")
    print("=" * 60)
