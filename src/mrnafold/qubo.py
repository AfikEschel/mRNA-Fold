"""
QUBO formulation for mRNA secondary structure prediction.

This module implements the conversion of quartet-based preprocessing
into a Quadratic Unconstrained Binary Optimization (QUBO) problem,
following the IBM–Moderna paper approach.
"""

from typing import List, Dict, Tuple, Set
import numpy as np
from mrnafold.quartets import Quartet, find_conflicting_quartets, find_stackable_quartets


def build_quartet_qubo(quartets: List[Quartet],
                       conflict_penalty: float = 1.0,
                       stacking_reward: float = 0.5) -> np.ndarray:
    """
    Build QUBO matrix from quartet preprocessing.

    The QUBO formulation uses binary variables x_q where x_q = 1 if quartet q
    is selected in the structure, 0 otherwise.

    Objective function:
    H = sum_{q in Q} h_q * x_q + sum_{(q,r) in conflicts} J_{qr} * x_q * x_r

    Where:
    - h_q: linear coefficient (can be used for quartet energy)
    - J_{qr}: quadratic coefficient (penalty for conflicts, reward for stacking)

    Args:
        quartets: List of Quartet objects
        conflict_penalty: Penalty for selecting conflicting quartets
        stacking_reward: Reward for selecting stacking quartets

    Returns:
        QUBO matrix as numpy array (size n_quartets x n_quartets)
    """
    n = len(quartets)
    qubo = np.zeros((n, n))

    # Get conflict and stacking relationships
    conflicts = find_conflicting_quartets(quartets)
    stacking = find_stackable_quartets(quartets)

    # Build quadratic terms
    for i in range(n):
        for j in range(i + 1, n):
            qi, qj = quartets[i], quartets[j]

            # Check if these quartets conflict
            if (qi, qj) in conflicts or (qj, qi) in conflicts:
                qubo[i, j] = qubo[j, i] = conflict_penalty

            # Check if these quartets can stack
            if qj in stacking.get(qi, []) or qi in stacking.get(qj, []):
                # Stacking reward (negative coefficient since we minimize)
                qubo[i, j] = qubo[j, i] = -stacking_reward

    return qubo


def qubo_to_ising(qubo: np.ndarray) -> Tuple[np.ndarray, float]:
    """
    Convert QUBO to Ising model form.

    Ising form: H = sum_i h_i * s_i + sum_{i<j} J_{ij} * s_i * s_j
    where s_i = {-1, +1}

    QUBO form: H = sum_i h_i * x_i + sum_{i<j} J_{ij} * x_i * x_j
    where x_i = {0, 1}

    Conversion: s_i = 2*x_i - 1
    Then: h_i = h_i_qubo + sum_j J_{ij}_qubo
         J_{ij} = J_{ij}_qubo / 4

    Args:
        qubo: QUBO matrix

    Returns:
        Tuple of (Ising matrix, constant offset)
    """
    n = qubo.shape[0]
    ising = np.zeros((n, n))

    # Convert quadratic terms for upper-triangular QUBO
    for i in range(n):
        for j in range(i + 1, n):
            ising[i, j] = ising[j, i] = qubo[i, j] / 2

    # Linear terms become diagonal in Ising
    for i in range(n):
        linear_term = qubo[i, i] / 2 + sum(qubo[i, j] / 2 for j in range(n) if j != i)
        ising[i, i] = linear_term

    # Constant offset uses upper-triangular entries only
    constant = sum(qubo[i, i] for i in range(n)) / 2 + sum(qubo[i, j] / 2 for i in range(n) for j in range(i + 1, n))

    return ising, constant


def solve_qubo_brute_force(qubo: np.ndarray) -> Tuple[np.ndarray, float]:
    """
    Solve small QUBO instances using brute force.

    This is for validation and testing on tiny instances only.
    Complexity is O(2^n), so only works for n <= 20 or so.

    Args:
        qubo: QUBO matrix

    Returns:
        Tuple of (best_solution vector, best_energy)
    """
    n = qubo.shape[0]
    if n > 20:
        raise ValueError(f"Brute force only works for n <= 20, got n={n}")

    if n == 0:
        return np.array([], dtype=int), 0.0

    best_energy = float('inf')
    best_solution = None

    # Try all possible binary assignments
    for assignment in range(2**n):
        # Convert to binary vector
        x = np.array([int(bit) for bit in format(assignment, f'0{n}b')])

        # Compute energy: x^T @ Q @ x
        energy = x @ qubo @ x

        if energy < best_energy:
            best_energy = energy
            best_solution = x.copy()

    return best_solution, best_energy


def quartets_to_structure(selected_quartets: List[Quartet],
                         sequence_length: int) -> List[Tuple[int, int]]:
    """
    Convert selected quartets to base pair structure.

    Args:
        selected_quartets: List of selected Quartet objects
        sequence_length: Length of the RNA sequence

    Returns:
        List of (i,j) base pairs
    """
    base_pairs = set()

    for quartet in selected_quartets:
        # Add both pairs from the quartet
        base_pairs.add((quartet.i, quartet.j))      # outer pair
        base_pairs.add((quartet.inner_i, quartet.inner_j))  # inner pair

    return sorted(list(base_pairs))


def structure_to_dot_bracket(base_pairs: List[Tuple[int, int]],
                           sequence_length: int) -> str:
    """
    Convert base pairs to dot-bracket notation.

    Args:
        base_pairs: List of (i,j) base pairs (0-based indexing)
        sequence_length: Length of sequence

    Returns:
        Dot-bracket string
    """
    structure = ['.'] * sequence_length

    for i, j in base_pairs:
        if structure[i] != '.' or structure[j] != '.':
            # Conflict - this shouldn't happen with valid quartet selection
            continue
        structure[i] = '('
        structure[j] = ')'

    return ''.join(structure)