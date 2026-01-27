#!/usr/bin/env python3
"""
Example script demonstrating symplectic eigenvalues computation

This script shows how to use the symplectic eigenvalues functions
with various types of positive definite matrices.
"""

import numpy as np
from symplectic_eigenvalues import (
    symplectic_eigenvalues,
    symplectic_diagonalizing_matrix,
    williamson_decomposition,
    symplectic_form,
    verify_symplectic
)


def example1_diagonal_matrix():
    """Example 1: Diagonal matrix"""
    print("\n" + "="*70)
    print("Example 1: Diagonal Matrix")
    print("="*70)
    
    # Create a diagonal positive definite matrix
    A = np.diag([2.0, 3.0, 2.0, 3.0])
    print("\nMatrix A:")
    print(A)
    
    # Compute symplectic eigenvalues
    evals = symplectic_eigenvalues(A)
    print(f"\nSymplectic eigenvalues: {evals}")
    print(f"Number of unique eigenvalues: {len(evals)}")


def example2_block_diagonal():
    """Example 2: Block diagonal matrix"""
    print("\n" + "="*70)
    print("Example 2: Block Diagonal Matrix")
    print("="*70)
    
    # Create a block diagonal matrix
    # Two blocks: [[3, 1], [1, 3]] and [[5, 2], [2, 5]]
    A = np.array([[3, 1, 0, 0], 
                  [1, 3, 0, 0], 
                  [0, 0, 5, 2], 
                  [0, 0, 2, 5]], dtype=float)
    print("\nMatrix A (block diagonal):")
    print(A)
    
    # Compute eigenvalues of the blocks
    block1 = A[:2, :2]
    block2 = A[2:, 2:]
    evals_block1 = np.linalg.eigvalsh(block1)
    evals_block2 = np.linalg.eigvalsh(block2)
    
    print(f"\nEigenvalues of first block: {evals_block1}")
    print(f"Eigenvalues of second block: {evals_block2}")
    
    # Compute symplectic eigenvalues
    evals = symplectic_eigenvalues(A)
    print(f"\nSymplectic eigenvalues: {evals}")


def example3_diagonalization():
    """Example 3: Symplectic diagonalization"""
    print("\n" + "="*70)
    print("Example 3: Symplectic Diagonalization")
    print("="*70)
    
    # Create a simple positive definite matrix
    A = np.array([[4, 1, 0, 0], 
                  [1, 4, 0, 0], 
                  [0, 0, 4, 1], 
                  [0, 0, 1, 4]], dtype=float)
    print("\nOriginal matrix A:")
    print(A)
    
    # Get the symplectic diagonalizing matrix
    S, D = symplectic_diagonalizing_matrix(A)
    
    print("\nDiagonalizing matrix S (shape):", S.shape)
    print("Determinant of S:", np.linalg.det(S))
    
    print("\nDiagonal matrix D (real part):")
    print(np.round(D.real, 3))
    
    # Verify that D is approximately diagonal
    off_diagonal = np.sum(np.abs(D - np.diag(np.diag(D))))
    print(f"\nSum of off-diagonal elements in D: {off_diagonal:.2e}")
    if off_diagonal < 1e-10:
        print("✓ D is effectively diagonal")


def example4_williamson():
    """Example 4: Williamson decomposition"""
    print("\n" + "="*70)
    print("Example 4: Williamson Decomposition")
    print("="*70)
    
    # Create a covariance matrix (positive definite)
    A = np.array([[2.5, 0.3, 0.2, 0.1], 
                  [0.3, 2.5, 0.1, 0.2], 
                  [0.2, 0.1, 2.5, 0.3], 
                  [0.1, 0.2, 0.3, 2.5]], dtype=float)
    print("\nCovariance matrix A:")
    print(A)
    
    # Check that A is positive definite
    eigenvalues_A = np.linalg.eigvalsh(A)
    print(f"\nEigenvalues of A: {eigenvalues_A}")
    print(f"Minimum eigenvalue: {np.min(eigenvalues_A):.4f}")
    if np.min(eigenvalues_A) > 0:
        print("✓ A is positive definite")
    
    # Williamson decomposition
    S, d = williamson_decomposition(A)
    print(f"\nSymplectic eigenvalues from Williamson: {d}")
    print(f"Product of symplectic eigenvalues: {np.prod(d):.4f}")


def example5_symplectic_form():
    """Example 5: Symplectic form properties"""
    print("\n" + "="*70)
    print("Example 5: Symplectic Form Properties")
    print("="*70)
    
    n = 3
    Omega = symplectic_form(n)
    
    print(f"\nSymplectic form Ω for n={n} (size {2*n}×{2*n}):")
    print(Omega)
    
    # Check antisymmetry
    print("\nChecking antisymmetry: Ω^T = -Ω")
    print(f"||Ω^T + Ω|| = {np.linalg.norm(Omega.T + Omega):.2e}")
    if np.allclose(Omega.T, -Omega):
        print("✓ Ω is antisymmetric")
    
    # Check Ω² = -I
    Omega2 = Omega @ Omega
    neg_I = -np.eye(2*n)
    print(f"\nChecking Ω² = -I:")
    print(f"||Ω² + I|| = {np.linalg.norm(Omega2 - neg_I):.2e}")
    if np.allclose(Omega2, neg_I):
        print("✓ Ω² = -I")


def example6_random_matrix():
    """Example 6: Random positive definite matrix"""
    print("\n" + "="*70)
    print("Example 6: Random Positive Definite Matrix")
    print("="*70)
    
    # Generate a random positive definite matrix
    np.random.seed(42)  # For reproducibility
    n = 2
    # Create a random matrix and make it positive definite
    M = np.random.randn(2*n, 2*n)
    A = M.T @ M + np.eye(2*n)  # Ensures positive definiteness
    
    print("\nRandom positive definite matrix A:")
    print(np.round(A, 3))
    
    # Check eigenvalues
    evals_A = np.linalg.eigvalsh(A)
    print(f"\nEigenvalues of A: {np.round(evals_A, 3)}")
    
    # Compute symplectic eigenvalues
    sym_evals = symplectic_eigenvalues(A)
    print(f"\nSymplectic eigenvalues: {np.round(sym_evals, 3)}")


def main():
    """Run all examples"""
    print("\n" + "="*70)
    print("SYMPLECTIC EIGENVALUES - EXAMPLES")
    print("="*70)
    print("\nThis script demonstrates various uses of symplectic eigenvalue")
    print("computation for positive definite Hermitian matrices.")
    
    example1_diagonal_matrix()
    example2_block_diagonal()
    example3_diagonalization()
    example4_williamson()
    example5_symplectic_form()
    example6_random_matrix()
    
    print("\n" + "="*70)
    print("All examples completed successfully!")
    print("="*70)


if __name__ == "__main__":
    main()
