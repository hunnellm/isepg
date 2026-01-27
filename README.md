# isepg
Functions for computations related to the inverse symplectic eigenvalue for a graph

## Symplectic Eigenvalues Computation

This repository provides implementations for computing:
- Symplectic eigenvalues of positive definite Hermitian matrices
- Symplectic matrices that diagonalize these matrices
- Williamson decomposition

### Quick Start

```bash
# Install dependencies
pip install -r requirements.txt

# Run tests
python3 test_symplectic_eigenvalues.py

# Run examples
python3 examples.py
```

### Files

- `symplectic_eigenvalues.py` - Python/NumPy implementation
- `symplectic_eigenvalues.sage` - SageMath implementation
- `symplectic_eigenvalues.mac` - Maxima implementation (symbolic computation)
- `test_symplectic_eigenvalues.py` - Test suite
- `examples_maxima.mac` - Maxima examples
- `SYMPLECTIC_EIGENVALUES.md` - Detailed documentation

### Example Usage

```python
import numpy as np
from symplectic_eigenvalues import symplectic_eigenvalues

# Compute symplectic eigenvalues of a matrix
A = np.array([[3, 1, 0, 0], 
              [1, 3, 0, 0], 
              [0, 0, 3, 1], 
              [0, 0, 1, 3]], dtype=float)

eigenvalues = symplectic_eigenvalues(A)
print("Symplectic eigenvalues:", eigenvalues)
# Output: [2.0, 4.0]
```

For more details, see [SYMPLECTIC_EIGENVALUES.md](SYMPLECTIC_EIGENVALUES.md).
