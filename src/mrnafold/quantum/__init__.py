"""
Quantum module for mRNA secondary structure prediction.

Phase 10: Full CVaR-VQE implementation

This module provides:
- BaseSolver: Abstract interface for pluggable solvers
- TwoLocalAnsatz: Hardware-efficient ansatz (Y rotations + CZ gates, p=2 layers)
- CVaRVQE: Conditional Value at Risk Variational Quantum Eigensolver

Runs on Qiskit's statevector simulator for small instances.
"""

from .base_solver import BaseSolver
from .cvar_vqe import CVaRVQE
from .ansatz import TwoLocalAnsatz

__all__ = ["BaseSolver", "CVaRVQE", "TwoLocalAnsatz"]
__version__ = "0.2.0"  # Phase 10 complete
