# Symplectic Eigenvalues Computation

This repository provides functions to compute the symplectic eigenvalues of positive definite Hermitian matrices and the symplectic matrices that diagonalize them.

## Overview

For a positive definite Hermitian matrix **A** of size 2n×2n, the **symplectic eigenvalues** are computed from the eigenvalues of the matrix **iΩA**, where **Ω** is the standard symplectic form matrix:

```
Ω = [[0, I_n], [-I_n, 0]]
```

The eigenvalues of **iΩA** come in pairs **±λ**, where **λ** are the symplectic eigenvalues.

## Files

- `symplectic_eigenvalues.py` - Python/NumPy implementation (works without Sage)
- `symplectic_eigenvalues.sage` - SageMath implementation (requires Sage)
- `symplectic_eigenvalues.mac` - Maxima implementation (symbolic computation)
- `examples_maxima.mac` - Comprehensive Maxima examples

## Installation

### For Python version (NumPy/SciPy):

```bash
pip install numpy scipy
```

### For Sage version:

Install SageMath from https://www.sagemath.org/

### For Maxima version:

Install Maxima from https://maxima.sourceforge.io/ or using your package manager:

```bash
# Ubuntu/Debian
sudo apt-get install maxima

# macOS
brew install maxima

# Fedora
sudo dnf install maxima
```

## Usage

### Python/NumPy Version

```python
import numpy as np
from symplectic_eigenvalues import (
    symplectic_eigenvalues,
    symplectic_diagonalizing_matrix,
    williamson_decomposition,
    symplectic_form
)

# Example 1: Compute symplectic eigenvalues of a diagonal matrix
A = np.eye(4) * 2  # 4x4 identity matrix scaled by 2
eigenvalues = symplectic_eigenvalues(A)
print("Symplectic eigenvalues:", eigenvalues)
# Output: [2.0]

# Example 2: Compute symplectic eigenvalues of a block diagonal matrix
A = np.array([[3, 1, 0, 0], 
              [1, 3, 0, 0], 
              [0, 0, 3, 1], 
              [0, 0, 1, 3]], dtype=float)
eigenvalues = symplectic_eigenvalues(A)
print("Symplectic eigenvalues:", eigenvalues)
# Output: [2.0, 4.0]

# Example 3: Get the symplectic matrix that diagonalizes A
S, D = symplectic_diagonalizing_matrix(A)
print("Diagonalizing matrix S shape:", S.shape)
print("Diagonal matrix D:\n", D)

# Example 4: Williamson decomposition
S, d = williamson_decomposition(A)
print("Symplectic eigenvalues:", d)
```

### SageMath Version

```python
# Load the Sage file
load("symplectic_eigenvalues.sage")

# Example usage
from sage.all import matrix

# Create a positive definite matrix
A = matrix([[3, 1, 0, 0], 
            [1, 3, 0, 0], 
            [0, 0, 3, 1], 
            [0, 0, 1, 3]])

# Compute symplectic eigenvalues
eigenvalues = symplectic_eigenvalues(A)
print("Symplectic eigenvalues:", eigenvalues)

# Get the symplectic diagonalizing matrix
S, D = symplectic_diagonalizing_matrix(A)
print("Diagonal matrix D:")
print(D)
```

### Maxima Version (Symbolic Computation)

```maxima
/* Load the Maxima file */
load("symplectic_eigenvalues.mac");

/* Example 1: Compute symplectic eigenvalues */
A: matrix([3,1,0,0], [1,3,0,0], [0,0,3,1], [0,0,1,3]);
eigenvalues: symplectic_eigenvalues(A);
/* Output: [2, 4] */

/* Example 2: Symbolic computation with parameters */
/* Define a matrix with symbolic parameter a */
A_symbolic: matrix([a,0,0,0], [0,a,0,0], [0,0,a,0], [0,0,0,a]);
evals: symplectic_eigenvalues(A_symbolic);
/* Output: [abs(a)] - symbolic result! */

/* Example 3: Get the symplectic form matrix */
Omega: symplectic_form(2);
/* Output: 4×4 antisymmetric matrix */

/* Example 4: Williamson decomposition */
result: williamson_decomposition(A);
S: first(result);
d: second(result);
/* S is the symplectic matrix, d is the list of eigenvalues */

/* Example 5: Run comprehensive examples */
/* maxima --batch=examples_maxima.mac */
```

## Functions

### `symplectic_form(n)`
Creates the standard 2n×2n symplectic form matrix Ω.

**Parameters:**
- `n`: Half the dimension of the symplectic form matrix

**Returns:**
- The 2n×2n symplectic form matrix

### `symplectic_eigenvalues(A)`
Computes the symplectic eigenvalues of a positive definite Hermitian matrix A.

**Parameters:**
- `A`: A positive definite Hermitian matrix of size 2n×2n

**Returns:**
- Array/list of unique symplectic eigenvalues (positive real numbers), sorted

### `symplectic_diagonalizing_matrix(A)`
Computes the symplectic matrix S that diagonalizes the matrix A.

**Parameters:**
- `A`: A positive definite Hermitian matrix of size 2n×2n

**Returns:**
- `S`: The symplectic matrix that diagonalizes A
- `D`: The diagonal matrix S^T A S (or S^(-1) A S)

### `williamson_decomposition(A)`
Computes the Williamson decomposition of a positive definite matrix.

**Parameters:**
- `A`: A positive definite Hermitian matrix of size 2n×2n

**Returns:**
- `S`: The symplectic matrix
- `d`: Vector of symplectic eigenvalues

## Mathematical Background

The symplectic eigenvalues have applications in:
- Quantum information theory (entanglement measures)
- Quantum optics (symplectic transformations of Gaussian states)
- Control theory
- Hamiltonian mechanics

For a positive definite matrix **A**, the symplectic eigenvalues satisfy:
- They are always positive real numbers
- They appear with specific multiplicity determined by the structure of **A**
- They are invariant under symplectic transformations

## Testing

Run the test suite:

```bash
# Python version
python3 symplectic_eigenvalues.py

# Sage version (if Sage is installed)
sage symplectic_eigenvalues.sage
```

## References

- Williamson, J. (1936). "On the Algebraic Problem Concerning the Normal Forms of Linear Dynamical Systems". American Journal of Mathematics.
- Simon, R., et al. (1994). "Quantum-noise matrix for multimode systems: U(n) invariance, squeezing, and normal forms". Physical Review A.

## License

This code is provided as-is for educational and research purposes.

## Contributing

This repository is for computations related to the inverse symplectic eigenvalue for graphs. Contributions are welcome!
