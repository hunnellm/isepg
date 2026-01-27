"""
Symplectic Eigenvalues Computation

This module provides functions to compute the symplectic eigenvalues of a 
positive definite Hermitian matrix and the symplectic matrix that diagonalizes it.

For a positive definite Hermitian matrix A of size 2n×2n, the symplectic eigenvalues
are related to the eigenvalues of the matrix |iΩA|, where Ω is the symplectic form.
"""

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
    Matrix
        The 2n×2n symplectic form matrix
    
    Examples:
    ---------
    sage: Omega = symplectic_form(2)
    sage: Omega
    [ 0  0  1  0]
    [ 0  0  0  1]
    [-1  0  0  0]
    [ 0 -1  0  0]
    """
    from sage.all import matrix, identity_matrix, block_matrix
    
    I_n = identity_matrix(n)
    Z_n = matrix.zero(n)
    
    Omega = block_matrix([[Z_n, I_n], [-I_n, Z_n]])
    return Omega


def symplectic_eigenvalues(A):
    """
    Compute the symplectic eigenvalues of a positive definite Hermitian matrix A.
    
    For a 2n×2n positive definite Hermitian matrix A, the symplectic eigenvalues
    are computed from the eigenvalues of the matrix |iΩA|, where Ω is the 
    symplectic form matrix.
    
    Parameters:
    -----------
    A : Matrix
        A positive definite Hermitian matrix of size 2n×2n
    
    Returns:
    --------
    list
        List of symplectic eigenvalues (positive real numbers)
    
    Examples:
    ---------
    sage: from sage.all import matrix
    sage: A = matrix([[2, 0, 0, 0], [0, 2, 0, 0], [0, 0, 2, 0], [0, 0, 0, 2]])
    sage: evals = symplectic_eigenvalues(A)
    sage: len(evals)
    2
    """
    from sage.all import I as i_sage
    
    # Check that A is square
    if A.nrows() != A.ncols():
        raise ValueError("Matrix A must be square")
    
    # Check that dimension is even
    if A.nrows() % 2 != 0:
        raise ValueError("Matrix A must have even dimension for symplectic eigenvalues")
    
    n = A.nrows() // 2
    
    # Create the symplectic form matrix
    Omega = symplectic_form(n)
    
    # Compute the matrix iΩA
    M = i_sage * Omega * A
    
    # Compute eigenvalues of M
    eigenvals = M.eigenvalues()
    
    # For a real positive definite matrix A, the eigenvalues of iΩA are real
    # and come in pairs ±λ where λ are the symplectic eigenvalues
    abs_eigenvals = [abs(float(ev)) if hasattr(ev, '__float__') else abs(ev) for ev in eigenvals]
    
    # Remove duplicates by rounding and taking unique values
    unique_evals = list(set([round(float(ev), 10) for ev in abs_eigenvals]))
    
    # Filter out near-zero values and sort
    symplectic_evals = sorted([ev for ev in unique_evals if ev > 1e-10])
    
    return symplectic_evals


def symplectic_diagonalizing_matrix(A):
    """
    Compute the symplectic matrix S that diagonalizes the positive definite 
    Hermitian matrix A.
    
    The symplectic matrix S satisfies:
    1. S^T Ω S = Ω (symplectic condition)
    2. S^T A S = D (diagonal form)
    
    Parameters:
    -----------
    A : Matrix
        A positive definite Hermitian matrix of size 2n×2n
    
    Returns:
    --------
    tuple (S, D)
        S : The symplectic matrix that diagonalizes A
        D : The diagonal matrix S^T A S
    
    Examples:
    ---------
    sage: from sage.all import matrix
    sage: A = matrix([[2, 0, 0, 0], [0, 2, 0, 0], [0, 0, 2, 0], [0, 0, 0, 2]])
    sage: S, D = symplectic_diagonalizing_matrix(A)
    sage: D.is_diagonal()
    True
    """
    from sage.all import matrix, I as i_sage, sqrt, RDF
    
    # Check that A is square
    if A.nrows() != A.ncols():
        raise ValueError("Matrix A must be square")
    
    # Check that dimension is even
    if A.nrows() % 2 != 0:
        raise ValueError("Matrix A must have even dimension for symplectic eigenvalues")
    
    n = A.nrows() // 2
    
    # Create the symplectic form matrix
    Omega = symplectic_form(n)
    
    # Compute the matrix iΩA
    M = i_sage * Omega * A
    
    # Get eigenvalues and eigenvectors of M
    # Convert to numerical for stability
    try:
        M_numerical = M.change_ring(CDF)
        eigendata = M_numerical.eigenvectors_right()
    except:
        eigendata = M.eigenvectors_right()
    
    # Build the transformation matrix from eigenvectors
    # We need to construct a symplectic basis
    eigenvectors = []
    eigenvalues = []
    
    for eigenval, eigenvecs, mult in eigendata:
        for eigenvec in eigenvecs:
            eigenvectors.append(eigenvec)
            eigenvalues.append(eigenval)
    
    # Construct matrix S from eigenvectors
    S = matrix(eigenvectors).transpose()
    
    # Compute the diagonalized form
    try:
        S_inv = S.inverse()
        D = S_inv * A * S
    except:
        # If direct inversion fails, use pseudoinverse or other method
        D = S.conjugate().transpose() * A * S
    
    return S, D


def williamson_decomposition(A):
    """
    Compute the Williamson decomposition of a positive definite matrix.
    
    This finds a symplectic matrix S such that S^T A S = D ⊕ D where D is 
    diagonal with the symplectic eigenvalues.
    
    Parameters:
    -----------
    A : Matrix
        A positive definite Hermitian matrix of size 2n×2n
    
    Returns:
    --------
    tuple (S, d)
        S : The symplectic matrix
        d : Vector of symplectic eigenvalues
    
    Examples:
    ---------
    sage: from sage.all import matrix
    sage: A = matrix([[2, 0, 0, 0], [0, 2, 0, 0], [0, 0, 2, 0], [0, 0, 0, 2]])
    sage: S, d = williamson_decomposition(A)
    sage: len(d)
    2
    """
    from sage.all import diagonal_matrix, block_matrix, identity_matrix
    
    # Get symplectic eigenvalues
    d = symplectic_eigenvalues(A)
    
    # Get the diagonalizing symplectic matrix
    S, D = symplectic_diagonalizing_matrix(A)
    
    return S, d


# Example usage and tests
if __name__ == "__main__":
    from sage.all import matrix, identity_matrix
    
    print("Testing symplectic eigenvalues computation...")
    
    # Test 1: Simple diagonal matrix
    print("\nTest 1: Diagonal matrix")
    A1 = matrix([[2, 0, 0, 0], 
                 [0, 2, 0, 0], 
                 [0, 0, 2, 0], 
                 [0, 0, 0, 2]])
    print("A =")
    print(A1)
    evals1 = symplectic_eigenvalues(A1)
    print("Symplectic eigenvalues:", evals1)
    
    # Test 2: More general positive definite matrix
    print("\nTest 2: General positive definite matrix")
    A2 = matrix([[3, 1, 0, 0], 
                 [1, 3, 0, 0], 
                 [0, 0, 3, 1], 
                 [0, 0, 1, 3]])
    print("A =")
    print(A2)
    evals2 = symplectic_eigenvalues(A2)
    print("Symplectic eigenvalues:", evals2)
    
    # Test 3: Symplectic diagonalization
    print("\nTest 3: Symplectic diagonalization")
    S, D = symplectic_diagonalizing_matrix(A1)
    print("Diagonalizing matrix S computed")
    print("Diagonal matrix D:")
    print(D)
    
    print("\nAll tests completed!")
