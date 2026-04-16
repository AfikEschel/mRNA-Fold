# src/mrnafold/quantum/base_solver.py
"""
Abstract base class for QUBO solvers.

Defines the interface that all solvers (classical, quantum, etc.) should implement.
This enables pluggable solver architecture for future extensibility.
"""

from abc import ABC, abstractmethod
from typing import Tuple
import numpy as np


class BaseSolver(ABC):
    """
    Abstract base class for QUBO solvers.

    All solvers should implement solve() to take a QUBO matrix and return
    the best bitstring found and its energy.
    """

    @abstractmethod
    def solve(self, qubo: np.ndarray) -> Tuple[np.ndarray, float]:
        """
        Solve a QUBO problem.

        Args:
            qubo: QUBO matrix (n x n) representing the optimization problem.

        Returns:
            Tuple of (best_bitstring, best_energy)
                - best_bitstring: np.ndarray of shape (n,) with binary values {0, 1}
                - best_energy: float, energy of the best solution found
        """
        raise NotImplementedError("Subclasses must implement solve().")
