# src/mrnafold/quantum/ansatz.py
"""
Placeholder for ansatz implementation.

Inspired by the IBM-Moderna paper: Uses hardware-efficient "two-local" ansatz
with single-qubit Y rotations and two-qubit CZ gates.

Compatible with NFT optimizer.

This is a scaffold - NOT implemented yet. Future work.
"""

from typing import List
import numpy as np


class TwoLocalAnsatz:
    """
    Placeholder class for two-local ansatz.

    In the paper: Single-qubit Pauli-Y rotations and CZ gates.
    Layers: Even qubits entangled with next, then odd.
    """

    def __init__(self, num_qubits: int, depth: int = 1):
        """
        Initialize two-local ansatz.

        Args:
            num_qubits: Number of qubits (from QUBO variables).
            depth: Circuit depth (number of layers).
        """
        self.num_qubits = num_qubits
        self.depth = depth
        self.parameters = np.random.rand(num_qubits * depth)  # Placeholder

    def build_circuit(self, parameters: np.ndarray) -> object:
        """
        Placeholder for building the parameterized circuit.

        Args:
            parameters: Rotation angles.

        Returns:
            Placeholder circuit object.
        """
        # TODO: Implement with Qiskit or similar
        raise NotImplementedError("Ansatz circuit not implemented yet. This is scaffolding.")

    def get_parameter_count(self) -> int:
        """
        Get number of parameters.

        Returns:
            Number of parameters.
        """
        return len(self.parameters)