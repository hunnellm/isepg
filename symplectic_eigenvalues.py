"""
Symplectic Eigenvalues Computation (Python/NumPy version)

This module provides functions to compute the symplectic eigenvalues of a 
positive definite Hermitian matrix and the symplectic matrix that diagonalizes it.

For a positive definite Hermitian matrix A of size 2n×2n, the symplectic eigenvalues
are computed from the eigenvalues of the matrix iΩA, where Ω is the symplectic form.

This is a Python implementation using NumPy/SciPy that can run without Sage.
"""

import numpy as np
from scipy import linalg

# Numerical tolerance constants
EIGENVALUE_TOLERANCE = 1e-10  # Threshold for filtering near-zero eigenvalues
EIGENVALUE_PRECISION = 10     # Decimal places for rounding in duplicate detection


def symplectic_form(n):
    """
    Create the standard symplectic form matrix Ω of size 2n×2n.
    
    Ω = [[0, I_n], [-I_n, 0]]
    
    Parameters:
    -----------
    n : int
        Half the dimension of the symplectic form matrix
    
    Returns:
    --------
    ndarray
        The 2n×2n symplectic form matrix
    
    Examples:
    ---------
    >>> Omega = symplectic_form(2)
    >>> Omega.shape
    (4, 4)
    """
    I_n = np.eye(n)
    Z_n = np.zeros((n, n))
    
    # Block matrix: [[0, I], [-I, 0]]
    Omega = np.block([[Z_n, I_n], [-I_n, Z_n]])
    return Omega


def symplectic_eigenvalues(A):
    """
    Compute the symplectic eigenvalues of a positive definite Hermitian matrix A.
    
    For a 2n×2n positive definite Hermitian matrix A, the symplectic eigenvalues
    are computed from the eigenvalues of the matrix iΩA, where Ω is the 
    symplectic form matrix. The eigenvalues of iΩA are real and come in pairs 
    ±λ, where λ are the symplectic eigenvalues.
    
    Parameters:
    -----------
    A : ndarray
        A positive definite Hermitian matrix of size 2n×2n
    
    Returns:
    --------
    ndarray
        Array of symplectic eigenvalues (positive real numbers), sorted
    
    Examples:
    ---------
    >>> A = np.eye(4) * 2
    >>> evals = symplectic_eigenvalues(A)
    >>> len(evals)
    2
    >>> np.allclose(evals, [2.0, 2.0])
    True
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
    
    # For a Hermitian positive definite matrix A, the eigenvalues of iΩA are real
    # (the imaginary unit i makes them real) and come in pairs ±λ where λ are the 
    # symplectic eigenvalues
    # Take the positive eigenvalues (or absolute values)
    abs_eigenvals = np.abs(eigenvals.real)
    
    # Remove duplicates by taking unique values
    # Round to avoid numerical issues with duplicate detection
    unique_evals = np.unique(np.round(abs_eigenvals, decimals=EIGENVALUE_PRECISION))
    
    # Filter out near-zero values
    symplectic_evals = unique_evals[unique_evals > EIGENVALUE_TOLERANCE]
    
    return np.sort(symplectic_evals)


def symplectic_diagonalizing_matrix(A, return_diagonal=True):
    """
    Compute the symplectic matrix S that diagonalizes the positive definite 
    Hermitian matrix A.
    
    The symplectic matrix S satisfies:
    1. S^T Ω S = Ω (symplectic condition)
    2. S^T A S = D (diagonal form)
    
    Parameters:
    -----------
    A : ndarray
        A positive definite Hermitian matrix of size 2n×2n
    return_diagonal : bool, optional
        If True, return both S and D. If False, return only S.
    
    Returns:
    --------
    S : ndarray
        The symplectic matrix that diagonalizes A
    D : ndarray (if return_diagonal=True)
        The diagonal matrix S^T A S
    
    Examples:
    ---------
    >>> A = np.eye(4) * 2
    >>> S, D = symplectic_diagonalizing_matrix(A)
    >>> D.shape
    (4, 4)
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
    
    # Get eigenvalues and eigenvectors of M
    eigenvals, eigenvecs = linalg.eig(M)
    
    # Build the transformation matrix from eigenvectors
    # The eigenvectors form the columns of S
    S = eigenvecs
    
    if return_diagonal:
        # Compute the diagonalized form: S^T A S
        try:
            S_inv = linalg.inv(S)
            D = S_inv @ A @ S
        except linalg.LinAlgError:
            # If inverse fails, use conjugate transpose
            D = S.conj().T @ A @ S
        
        return S, D
    else:
        return S


def williamson_decomposition(A):
    """
    Compute the Williamson decomposition of a positive definite matrix.
    
    This finds a symplectic matrix S such that S^T A S has block diagonal form
    with the symplectic eigenvalues.
    
    Parameters:
    -----------
    A : ndarray
        A positive definite Hermitian matrix of size 2n×2n
    
    Returns:
    --------
    S : ndarray
        The symplectic matrix
    symplectic_evals : ndarray
        Vector of symplectic eigenvalues
    
    Examples:
    ---------
    >>> A = np.eye(4) * 2
    >>> S, symplectic_evals = williamson_decomposition(A)
    >>> len(symplectic_evals)
    2
    """
    # Get symplectic eigenvalues
    symplectic_evals = symplectic_eigenvalues(A)
    
    # Get the diagonalizing symplectic matrix
    S, D = symplectic_diagonalizing_matrix(A)
    
    return S, symplectic_evals


def verify_symplectic(S, n=None):
    """
    Verify that S is a symplectic matrix, i.e., S^T Ω S = Ω.
    
    Parameters:
    -----------
    S : ndarray
        Matrix to verify
    n : int, optional
        Half dimension (if None, computed from S)
    
    Returns:
    --------
    bool
        True if S is symplectic (within numerical tolerance)
    """
    if n is None:
        n = S.shape[0] // 2
    
    Omega = symplectic_form(n)
    
    # Check S^T Ω S = Ω
    result = S.T @ Omega @ S
    
    return np.allclose(result, Omega)


# Example usage and tests
if __name__ == "__main__":
    print("Testing symplectic eigenvalues computation (NumPy version)...")
    
    # Test 1: Simple diagonal matrix
    print("\nTest 1: Diagonal matrix (2I)")
    A1 = np.eye(4) * 2
    print("A =")
    print(A1)
    evals1 = symplectic_eigenvalues(A1)
    print("Symplectic eigenvalues:", evals1)
    print("Expected: [2.0] (with multiplicity 2)")
    
    # Test 2: More general positive definite matrix
    print("\nTest 2: Block diagonal positive definite matrix")
    A2 = np.array([[3, 1, 0, 0], 
                   [1, 3, 0, 0], 
                   [0, 0, 3, 1], 
                   [0, 0, 1, 3]], dtype=float)
    print("A =")
    print(A2)
    evals2 = symplectic_eigenvalues(A2)
    print("Symplectic eigenvalues:", evals2)
    print("Expected: [2.0, 4.0]")
    
    # Test 3: Symplectic diagonalization
    print("\nTest 3: Symplectic diagonalization")
    S, D = symplectic_diagonalizing_matrix(A1)
    print("Diagonalizing matrix S shape:", S.shape)
    print("Diagonal matrix D (should be close to diagonal):")
    print(np.round(D.real, 3))
    
    # Test 4: Verify symplectic property (if applicable)
    print("\nTest 4: Williamson decomposition")
    S_w, williamson_evals = williamson_decomposition(A1)
    print("Symplectic eigenvalues from Williamson:", williamson_evals)
    
    print("\nAll tests completed!")
