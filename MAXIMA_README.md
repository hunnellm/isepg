# Symplectic Eigenvalues - Maxima Implementation

This file contains the Maxima implementation for computing symplectic eigenvalues symbolically.

## What is Maxima?

Maxima is a computer algebra system (CAS) for symbolic mathematical computation. Unlike numerical implementations (Python/NumPy), Maxima can work with exact symbolic expressions, parameters, and algebraic expressions.

## Features

The Maxima implementation provides:

1. **Symbolic computation** - Work with exact symbolic expressions instead of floating-point approximations
2. **Parametric matrices** - Compute eigenvalues for matrices with symbolic parameters
3. **Exact arithmetic** - Results like `3/2` instead of `1.5`, or `abs(a)` instead of numerical approximation

## Installation

### Ubuntu/Debian
```bash
sudo apt-get install maxima
```

### macOS
```bash
brew install maxima
```

### Fedora
```bash
sudo dnf install maxima
```

### Windows
Download from: https://maxima.sourceforge.io/

## Quick Start

### Interactive Session

```bash
maxima
```

Then in the Maxima prompt:

```maxima
load("symplectic_eigenvalues.mac");

/* Define a matrix */
A: matrix([2,0,0,0], [0,2,0,0], [0,0,2,0], [0,0,0,2]);

/* Compute symplectic eigenvalues */
evals: symplectic_eigenvalues(A);
```

### Batch Mode

```bash
maxima --batch=examples_maxima.mac
```

### Single Command

```bash
maxima --very-quiet --batch-string="load(\"symplectic_eigenvalues.mac\"); A: matrix([3,1,0,0],[1,3,0,0],[0,0,3,1],[0,0,1,3]); symplectic_eigenvalues(A);"
```

## Available Functions

### `symplectic_form(n)`
Creates the standard 2n×2n symplectic form matrix Ω.

```maxima
Omega: symplectic_form(2);
/* Returns:
   [ 0   0   1  0 ]
   [ 0   0   0  1 ]
   [-1   0   0  0 ]
   [ 0  -1   0  0 ]
*/
```

### `symplectic_eigenvalues(A)`
Computes the symplectic eigenvalues of a positive definite Hermitian matrix A.

```maxima
A: matrix([3,1,0,0], [1,3,0,0], [0,0,3,1], [0,0,1,3]);
evals: symplectic_eigenvalues(A);
/* Returns: [2, 4] */
```

### `symplectic_diagonalizing_matrix(A)`
Computes the symplectic matrix S that diagonalizes the matrix A.

```maxima
A: matrix([2,0,0,0], [0,2,0,0], [0,0,2,0], [0,0,0,2]);
result: symplectic_diagonalizing_matrix(A);
S: first(result);   /* The symplectic matrix */
D: second(result);  /* The diagonalized form */
```

### `williamson_decomposition(A)`
Computes the Williamson decomposition of a positive definite matrix.

```maxima
A: matrix([3,1,0,0], [1,3,0,0], [0,0,3,1], [0,0,1,3]);
result: williamson_decomposition(A);
S: first(result);              /* The symplectic matrix */
eigenvals: second(result);     /* The symplectic eigenvalues */
```

### `verify_symplectic(S, n)`
Verifies that S is a symplectic matrix (S^T Ω S = Ω).

```maxima
S: matrix([1,0,0,0], [0,1,0,0], [0,0,1,0], [0,0,0,1]);
is_symplectic: verify_symplectic(S, 2);
```

## Symbolic Computation Examples

### Example 1: Exact Arithmetic

```maxima
/* Matrix with rational entries */
A: matrix([2, 1/2, 0, 0], 
          [1/2, 2, 0, 0], 
          [0, 0, 2, 1/2], 
          [0, 0, 1/2, 2]);
evals: symplectic_eigenvalues(A);
/* Result: [3/2, 5/2] - exact fractions! */
```

### Example 2: Parametric Matrix

```maxima
/* Matrix with symbolic parameter a */
A: matrix([a,0,0,0], [0,a,0,0], [0,0,a,0], [0,0,0,a]);
evals: symplectic_eigenvalues(A);
/* Result: [abs(a)] - symbolic result depending on parameter a */
```

### Example 3: Algebraic Expressions

```maxima
/* Matrix with algebraic expressions */
A: matrix([sqrt(2), 0, 0, 0],
          [0, sqrt(2), 0, 0],
          [0, 0, sqrt(2), 0],
          [0, 0, 0, sqrt(2)]);
evals: symplectic_eigenvalues(A);
/* Result contains sqrt(2) symbolically */
```

## Running the Examples

The file `examples_maxima.mac` contains comprehensive examples. Run it with:

```bash
cd /path/to/isepg
maxima --batch=examples_maxima.mac
```

This will demonstrate:
- Computing symplectic eigenvalues for various matrices
- Symbolic computation with parameters
- Creating symplectic form matrices
- Williamson decomposition
- Diagonalizing matrices

## Comparison with Other Implementations

| Feature | Python/NumPy | SageMath | Maxima |
|---------|--------------|----------|---------|
| Numerical computation | ✓ | ✓ | ✓ |
| Symbolic computation | ✗ | ✓ | ✓ |
| Exact arithmetic | ✗ | ✓ | ✓ |
| Parametric matrices | ✗ | ✓ | ✓ |
| Easy to install | ✓ | ✗ | ✓ |
| No dependencies | ✓ | ✗ | ✓ |

The Maxima implementation is particularly useful when you need:
- Exact symbolic results
- To work with parametric matrices
- To analyze how eigenvalues depend on parameters
- To avoid floating-point rounding errors

## Limitations

- For very large matrices (>10×10), symbolic computation can be slow
- Complex symbolic expressions may not simplify automatically
- Numerical evaluation may be needed for some complex symbolic results

## Tips

1. Use `ratsimp()` to simplify symbolic results:
   ```maxima
   evals: ratsimp(symplectic_eigenvalues(A));
   ```

2. Convert to floating-point if needed:
   ```maxima
   evals_float: float(evals);
   ```

3. Substitute specific values for parameters:
   ```maxima
   /* For matrix with parameter a */
   evals: symplectic_eigenvalues(A);
   /* Substitute a=3 */
   evals_a3: subst(3, a, evals);
   ```

## Further Reading

- Maxima documentation: https://maxima.sourceforge.io/docs/manual/maxima.html
- Eigenvalues in Maxima: https://maxima.sourceforge.io/docs/manual/maxima_25.html#SEC137
- Matrix operations: https://maxima.sourceforge.io/docs/manual/maxima_23.html

## License

This code is provided as-is for educational and research purposes.
