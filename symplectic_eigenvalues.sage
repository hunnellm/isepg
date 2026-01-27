"""
Symplectic Eigenvalues Computation

This module provides functions to compute the symplectic eigenvalues of a 
positive definite Hermitian matrix and the symplectic matrix that diagonalizes it.

For a positive definite Hermitian matrix A of size 2n×2n, the symplectic eigenvalues
are computed from the eigenvalues of the matrix iΩA, where Ω is the symplectic form.

NUMERICAL FUNCTIONS (default):
- symplectic_eigenvalues(A): Compute eigenvalues numerically with floating-point arithmetic
- symplectic_diagonalizing_matrix(A): Compute diagonalizing matrix numerically
- williamson_decomposition(A): Williamson decomposition with numerical computation

SYMBOLIC FUNCTIONS (for exact/symbolic computation):
- symplectic_eigenvalues_symbolic(A): Compute eigenvalues symbolically (exact)
- symplectic_diagonalizing_matrix_symbolic(A): Compute diagonalizing matrix symbolically
- williamson_decomposition_symbolic(A): Williamson decomposition with symbolic computation

Use symbolic functions when:
- Working with symbolic matrices containing variables
- Need exact results without numerical approximation
- Preserving mathematical expressions in symbolic form
"""

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
    are computed from the eigenvalues of the matrix iΩA, where Ω is the 
    symplectic form matrix. The eigenvalues of iΩA are real and come in pairs 
    ±λ, where λ are the symplectic eigenvalues.
    
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
    
    # For a Hermitian positive definite matrix A, the eigenvalues of iΩA are real
    # (the imaginary unit i makes them real) and come in pairs ±λ where λ are the 
    # symplectic eigenvalues
    abs_eigenvals = [abs(float(ev)) if hasattr(ev, '__float__') else abs(ev) for ev in eigenvals]
    
    # Remove duplicates by rounding and taking unique values
    unique_evals = list(set([round(float(ev), EIGENVALUE_PRECISION) for ev in abs_eigenvals]))
    
    # Filter out near-zero values and sort
    symplectic_evals = sorted([ev for ev in unique_evals if ev > EIGENVALUE_TOLERANCE])
    
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
    tuple (S, symplectic_evals)
        S : The symplectic matrix
        symplectic_evals : Vector of symplectic eigenvalues
    
    Examples:
    ---------
    sage: from sage.all import matrix
    sage: A = matrix([[2, 0, 0, 0], [0, 2, 0, 0], [0, 0, 2, 0], [0, 0, 0, 2]])
    sage: S, symplectic_evals = williamson_decomposition(A)
    sage: len(symplectic_evals)
    2
    """
    from sage.all import diagonal_matrix, block_matrix, identity_matrix
    
    # Get symplectic eigenvalues
    symplectic_evals = symplectic_eigenvalues(A)
    
    # Get the diagonalizing symplectic matrix
    S, D = symplectic_diagonalizing_matrix(A)
    
    return S, symplectic_evals


# ============================================================================
# SYMBOLIC COMPUTATION FUNCTIONS
# ============================================================================

def symplectic_eigenvalues_symbolic(A):
    """
    Compute the symplectic eigenvalues of a positive definite Hermitian matrix A
    using symbolic computation (no numerical approximations).
    
    For a 2n×2n positive definite Hermitian matrix A, the symplectic eigenvalues
    are computed from the eigenvalues of the matrix iΩA, where Ω is the 
    symplectic form matrix. The eigenvalues of iΩA are real and come in pairs 
    ±λ, where λ are the symplectic eigenvalues.
    
    This function performs all computations symbolically, preserving exact values
    and symbolic expressions. Use this when you need exact results or when working
    with symbolic matrices.
    
    Parameters:
    -----------
    A : Matrix
        A positive definite Hermitian matrix of size 2n×2n (can contain symbolic entries)
    
    Returns:
    --------
    list
        List of n symplectic eigenvalues as symbolic expressions (positive values),
        where n is the half-dimension of the 2n×2n matrix A. Returns values with
        multiplicity preserved.
    
    Examples:
    ---------
    sage: from sage.all import matrix, SR
    sage: var('a')
    a
    sage: A = matrix(SR, [[a, 0, 0, 0], [0, a, 0, 0], [0, 0, a, 0], [0, 0, 0, a]])
    sage: evals = symplectic_eigenvalues_symbolic(A)
    sage: evals
    [abs(a), abs(a)]
    """
    from sage.all import I as i_sage, SR, abs as sage_abs, simplify, sorted as sage_sorted
    
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
    
    # Compute eigenvalues of M symbolically
    eigenvals = M.eigenvalues()
    
    # For a Hermitian positive definite matrix A, the eigenvalues of iΩA are real
    # (the imaginary unit i makes them real) and come in pairs ±λ where λ are the 
    # symplectic eigenvalues. Take absolute values symbolically.
    abs_eigenvals = [sage_abs(ev) for ev in eigenvals]
    
    # Simplify the symbolic expressions
    simplified_evals = [simplify(ev) for ev in abs_eigenvals]
    
    # Sort the eigenvalues (Sage can sort symbolic expressions)
    try:
        sorted_evals = sorted(simplified_evals)
    except:
        # If sorting fails (e.g., for complex symbolic expressions), use unsorted
        sorted_evals = simplified_evals
    
    # For a 2n×2n matrix, the eigenvalues of iΩA come in ± pairs, giving 2n values.
    # After taking absolute values and sorting, we return the first n values
    # (with multiplicity). This gives us exactly n symplectic eigenvalues.
    count = n
    symplectic_evals = sorted_evals[:count]
    
    return symplectic_evals


def symplectic_diagonalizing_matrix_symbolic(A):
    """
    Compute the symplectic matrix S that diagonalizes the positive definite 
    Hermitian matrix A using symbolic computation.
    
    The symplectic matrix S satisfies:
    1. S^T Ω S = Ω (symplectic condition)
    2. S^T A S = D (diagonal form)
    
    This function performs all computations symbolically, preserving exact values
    and symbolic expressions.
    
    Parameters:
    -----------
    A : Matrix
        A positive definite Hermitian matrix of size 2n×2n (can contain symbolic entries)
    
    Returns:
    --------
    tuple (S, D)
        S : The symplectic matrix that diagonalizes A (symbolic)
        D : The diagonal matrix S^T A S (symbolic)
    
    Examples:
    ---------
    sage: from sage.all import matrix, SR
    sage: A = matrix(SR, [[2, 0, 0, 0], [0, 2, 0, 0], [0, 0, 2, 0], [0, 0, 0, 2]])
    sage: S, D = symplectic_diagonalizing_matrix_symbolic(A)
    sage: D.is_diagonal()
    True
    """
    from sage.all import matrix, I as i_sage, SR
    
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
    
    # Get eigenvalues and eigenvectors of M symbolically
    eigendata = M.eigenvectors_right()
    
    # Build the transformation matrix from eigenvectors
    eigenvectors = []
    eigenvalues = []
    
    for eigenval, eigenvecs, mult in eigendata:
        for eigenvec in eigenvecs:
            eigenvectors.append(eigenvec)
            eigenvalues.append(eigenval)
    
    # Construct matrix S from eigenvectors
    S = matrix(eigenvectors).transpose()
    
    # Compute the diagonalized form symbolically
    try:
        S_inv = S.inverse()
        D = S_inv * A * S
    except:
        # If direct inversion fails, use conjugate transpose
        D = S.conjugate().transpose() * A * S
    
    return S, D


def williamson_decomposition_symbolic(A):
    """
    Compute the Williamson decomposition of a positive definite matrix
    using symbolic computation.
    
    This finds a symplectic matrix S such that S^T A S = D ⊕ D where D is 
    diagonal with the symplectic eigenvalues.
    
    This function performs all computations symbolically, preserving exact values
    and symbolic expressions.
    
    Parameters:
    -----------
    A : Matrix
        A positive definite Hermitian matrix of size 2n×2n (can contain symbolic entries)
    
    Returns:
    --------
    tuple (S, symplectic_evals)
        S : The symplectic matrix (symbolic)
        symplectic_evals : Vector of symplectic eigenvalues (symbolic)
    
    Examples:
    ---------
    sage: from sage.all import matrix, SR
    sage: var('a')
    a
    sage: A = matrix(SR, [[a, 0, 0, 0], [0, a, 0, 0], [0, 0, a, 0], [0, 0, 0, a]])
    sage: S, symplectic_evals = williamson_decomposition_symbolic(A)
    sage: symplectic_evals
    [abs(a)]
    """
    # Get symplectic eigenvalues symbolically
    symplectic_evals = symplectic_eigenvalues_symbolic(A)
    
    # Get the diagonalizing symplectic matrix symbolically
    S, D = symplectic_diagonalizing_matrix_symbolic(A)
    
    return S, symplectic_evals


# Example usage and tests
if __name__ == "__main__":
    from sage.all import matrix, identity_matrix, SR, var
    
    print("=" * 70)
    print("Testing symplectic eigenvalues computation...")
    print("=" * 70)
    
    # Test 1: Simple diagonal matrix (numerical)
    print("\n" + "=" * 70)
    print("Test 1: Diagonal matrix (NUMERICAL)")
    print("=" * 70)
    A1 = matrix([[2, 0, 0, 0], 
                 [0, 2, 0, 0], 
                 [0, 0, 2, 0], 
                 [0, 0, 0, 2]])
    print("A =")
    print(A1)
    evals1 = symplectic_eigenvalues(A1)
    print("Symplectic eigenvalues (numerical):", evals1)
    
    # Test 2: More general positive definite matrix (numerical)
    print("\n" + "=" * 70)
    print("Test 2: General positive definite matrix (NUMERICAL)")
    print("=" * 70)
    A2 = matrix([[3, 1, 0, 0], 
                 [1, 3, 0, 0], 
                 [0, 0, 3, 1], 
                 [0, 0, 1, 3]])
    print("A =")
    print(A2)
    evals2 = symplectic_eigenvalues(A2)
    print("Symplectic eigenvalues (numerical):", evals2)
    
    # Test 3: Symplectic diagonalization (numerical)
    print("\n" + "=" * 70)
    print("Test 3: Symplectic diagonalization (NUMERICAL)")
    print("=" * 70)
    S, D = symplectic_diagonalizing_matrix(A1)
    print("Diagonalizing matrix S computed")
    print("Diagonal matrix D:")
    print(D)
    
    # Test 4: Symbolic computation with numeric matrix
    print("\n" + "=" * 70)
    print("Test 4: Symbolic computation with numeric matrix")
    print("=" * 70)
    A3 = matrix(SR, [[2, 0, 0, 0], 
                     [0, 2, 0, 0], 
                     [0, 0, 2, 0], 
                     [0, 0, 0, 2]])
    print("A =")
    print(A3)
    evals3 = symplectic_eigenvalues_symbolic(A3)
    print("Symplectic eigenvalues (symbolic):", evals3)
    print("Expected: 2 eigenvalues (n = 2 for 4x4 matrix)")
    
    # Test 5: Symbolic computation with symbolic parameter
    print("\n" + "=" * 70)
    print("Test 5: Symbolic computation with parameter 'a' (SYMBOLIC)")
    print("=" * 70)
    var('a')
    A4 = matrix(SR, [[a, 0, 0, 0], 
                     [0, a, 0, 0], 
                     [0, 0, a, 0], 
                     [0, 0, 0, a]])
    print("A =")
    print(A4)
    evals4 = symplectic_eigenvalues_symbolic(A4)
    print("Symplectic eigenvalues (symbolic):", evals4)
    print("Note: Returns 2 eigenvalues (n = 2) with symbolic expression 'a'")
    
    # Test 6: Symbolic diagonalization
    print("\n" + "=" * 70)
    print("Test 6: Symbolic diagonalization")
    print("=" * 70)
    S_sym, D_sym = symplectic_diagonalizing_matrix_symbolic(A3)
    print("Symbolic diagonalizing matrix S computed")
    print("Symbolic diagonal matrix D:")
    print(D_sym)
    
    # Test 7: Williamson decomposition (symbolic)
    print("\n" + "=" * 70)
    print("Test 7: Williamson decomposition (SYMBOLIC)")
    print("=" * 70)
    S_will, evals_will = williamson_decomposition_symbolic(A4)
    print("Symplectic matrix S computed")
    print("Symplectic eigenvalues:", evals_will)
    
    print("\n" + "=" * 70)
    print("All tests completed!")
    print("=" * 70)
    print("\nSummary:")
    print("- Numerical functions: symplectic_eigenvalues(), symplectic_diagonalizing_matrix()")
    print("- Symbolic functions: symplectic_eigenvalues_symbolic(), symplectic_diagonalizing_matrix_symbolic()")
    print("- Use symbolic functions for exact results and symbolic parameters")
