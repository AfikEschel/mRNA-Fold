# src/mrnafold/quantum/ansatz.py
"""
Two-local ansatz implementation for CVaR-VQE.

Inspired by the IBM-Moderna paper: Hardware-efficient "two-local" ansatz
with single-qubit Pauli-Y rotations and two-qubit CZ (control-Z) gates.

Architecture:
- Layer 1: Single-qubit Y rotations on all qubits
- Layers 1 to p (p=2 in paper): Pairwise CZ entanglement + Y rotations
  - Even i in layer 1: CZ(i, i+1) for even i
  - Odd i in layer 2: CZ(i, i+1) for odd i
  - Repeat for p layers total
"""

from typing import Optional
import numpy as np
from qiskit import QuantumCircuit, QuantumRegister
from qiskit.circuit import ParameterVector


class TwoLocalAnsatz:
    """
    Two-local ansatz with single-qubit Y rotations and CZ gates.

    Structure from IBM-Moderna paper:
    - p layers (typically p=2 for simulations)
    - Each layer: RY rotations on all qubits, then pairwise CZ gates
    - CZ pattern alternates between even and odd pairs
    """

    def __init__(self, num_qubits: int, num_layers: int = 2, seed: Optional[int] = None):
        """
        Initialize two-local ansatz.

        Args:
            num_qubits: Number of qubits in the circuit.
            num_layers: Number of two-local layers (p in paper). Default: 2.
            seed: Random seed for parameter initialization. Default: None.
        """
        if num_qubits < 1:
            raise ValueError(f"num_qubits must be >= 1, got {num_qubits}")
        if num_layers < 1:
            raise ValueError(f"num_layers must be >= 1, got {num_layers}")

        self.num_qubits = num_qubits
        self.num_layers = num_layers
        self.seed = seed

        # Calculate number of parameters
        # Initial RY layer: num_qubits parameters
        # Each subsequent layer: num_qubits (RY) + variable CZ gates
        # For simplicity, we use num_qubits parameters per layer
        self.num_parameters = num_qubits * (num_layers + 1)

        # Initialize random parameters
        if seed is not None:
            np.random.seed(seed)
        self.initial_parameters = np.random.uniform(0, 2 * np.pi, self.num_parameters)

    def build_circuit(self, parameters: np.ndarray) -> QuantumCircuit:
        """
        Build the parameterized quantum circuit.

        Args:
            parameters: Array of rotation angles (must have length num_parameters).

        Returns:
            QuantumCircuit representing the ansatz.
        """
        if len(parameters) != self.num_parameters:
            raise ValueError(
                f"Expected {self.num_parameters} parameters, got {len(parameters)}"
            )

        qr = QuantumRegister(self.num_qubits, "q")
        circuit = QuantumCircuit(qr, name="TwoLocal")

        param_idx = 0

        # Initial RY layer
        for i in range(self.num_qubits):
            circuit.ry(parameters[param_idx], qr[i])
            param_idx += 1

        # p layers of entanglement and rotation
        for layer in range(self.num_layers):
            # CZ gates in pairwise fashion
            # Layer 0 (even pairs): CZ(0,1), CZ(2,3), ...
            # Layer 1 (odd pairs): CZ(1,2), CZ(3,4), ...
            if layer % 2 == 0:
                # Even layer: entangle even qubits with next
                for i in range(0, self.num_qubits - 1, 2):
                    circuit.cz(qr[i], qr[i + 1])
            else:
                # Odd layer: entangle odd qubits with next
                for i in range(1, self.num_qubits - 1, 2):
                    circuit.cz(qr[i], qr[i + 1])

            # RY rotations after entanglement
            for i in range(self.num_qubits):
                circuit.ry(parameters[param_idx], qr[i])
                param_idx += 1

        return circuit

    def get_parameter_count(self) -> int:
        """
        Get the number of parameters in the ansatz.

        Returns:
            Number of parameters.
        """
        return self.num_parameters