# src/mrnafold/quantum/cvar_vqe.py
"""
CVaR-VQE (Conditional Value at Risk Variational Quantum Eigensolver) implementation.

Based on the IBM-Moderna paper: Uses CVaR as the objective function instead of
standard expectation value, which can improve convergence on hardware.

CVaR(α) = average of the lower α-tail of the energy distribution.

For α=1 (standard): CVaR equals the expectation value <H>
For α→0: CVaR approaches the minimum energy found

Uses NFT optimizer (via scipy if scipy-available, fallback to scipy.optimize.minimize).
"""

from typing import Tuple, Optional, List
import numpy as np  

from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister
from qiskit_aer import AerSimulator

from .base_solver import BaseSolver
from .ansatz import TwoLocalAnsatz


class CVaRVQE(BaseSolver):
    """
    Conditional Value at Risk Variational Quantum Eigensolver.

    Minimizes the CVaR objective using a parameterized ansatz on a quantum simulator.
    """

    def __init__(
        self,
        num_qubits: int,
        num_layers: int = 2,
        alpha: float = 0.1,
        shots: int = 1024,
        optimizer: str = "scipy",
        max_iterations: int = 100,
        seed: Optional[int] = None,
    ):
        """
        Initialize CVaR-VQE solver.

        Args:
            num_qubits: Number of qubits/variables.
            num_layers: Number of ansatz layers (p in paper). Default: 2.
            alpha: CVaR parameter (0 < alpha <= 1). Default: 0.1 (paper uses this).
            shots: Number of measurement shots per circuit evaluation. Default: 1024.
            optimizer: Optimizer to use ("scipy" or "nft"). Default: "scipy".
            max_iterations: Maximum optimization iterations. Default: 100.
            seed: Random seed. Default: None.
        """
        if not (0 < alpha <= 1):
            raise ValueError(f"alpha must be in (0, 1], got {alpha}")
        if shots < 1:
            raise ValueError(f"shots must be >= 1, got {shots}")
        if max_iterations < 1:
            raise ValueError(f"max_iterations must be >= 1, got {max_iterations}")

        self.num_qubits = num_qubits
        self.num_layers = num_layers
        self.alpha = alpha
        self.shots = shots
        self.optimizer_type = optimizer
        self.max_iterations = max_iterations
        self.seed = seed

        # Create ansatz
        self.ansatz = TwoLocalAnsatz(num_qubits, num_layers, seed=seed)

        # Choose simulator method based on number of qubits.
        # Use statevector for small circuits, matrix_product_state for larger ones.
        self.simulator_method = "statevector" if num_qubits <= 20 else "matrix_product_state"
        self.simulator = AerSimulator(method=self.simulator_method)

        # Track optimization history
        self.iteration = 0
        self.energy_history: List[float] = []
        self.best_energy = float("inf")
        self.best_parameters = None
        self.best_bitstring = None

    def _create_measurement_circuit(self, parameters: np.ndarray) -> QuantumCircuit:
        """
        Create a circuit that prepares the state and measures in computational basis.

        Args:
            parameters: Ansatz parameters.

        Returns:
            QuantumCircuit with preparation and measurement.
        """
        ansatz_circuit = self.ansatz.build_circuit(parameters)

        # Create circuit with measurement
        qr = ansatz_circuit.qregs[0]  # Get qreg from ansatz
        cr = ClassicalRegister(self.num_qubits, "c")
        circuit = QuantumCircuit(qr, cr)

        # Compose ansatz circuit
        circuit.compose(ansatz_circuit, inplace=True)

        # Measure all qubits
        circuit.measure(range(self.num_qubits), range(self.num_qubits))

        return circuit

    def _evaluate_cvar_objective(
        self, parameters: np.ndarray, qubo: np.ndarray
    ) -> float:
        """
        Evaluate the CVaR objective for given parameters.

        Samples bitstrings by measuring the prepared state,
        computes their energies, and returns CVaR.

        Args:
            parameters: Ansatz parameters.
            qubo: QUBO matrix.

        Returns:
            CVaR energy value.
        """
        circuit = self._create_measurement_circuit(parameters)

        # Run on simulator
        job = self.simulator.run(circuit, shots=self.shots, seed_simulator=self.seed)
        result = job.result()
        counts = result.get_counts()

        # Extract energies from bitstrings
        energies = []
        bitstring_samples = []
        for bitstring_str, count in counts.items():
            # Convert bitstring to binary array (Qiskit returns as string, LSB first)
            bitstring = np.array([int(b) for b in reversed(bitstring_str)])
            energy = float(bitstring @ qubo @ bitstring)
            # Add samples according to their count
            energies.extend([energy] * count)
            bitstring_samples.append((bitstring, energy, count))

        energies = np.array(energies)

        # Compute CVaR
        sorted_energies = np.sort(energies)
        k = max(1, int(np.ceil(self.alpha * len(sorted_energies))))
        cvar = np.mean(sorted_energies[:k])

        # Track best solution
        # Find the bitstring with minimum energy from this round
        min_energy_idx = np.argmin(energies)
        min_energy = energies[min_energy_idx]

        if min_energy < self.best_energy:
            self.best_energy = min_energy
            self.best_parameters = parameters.copy()
            # Find the bitstring corresponding to min_energy
            for bs, e, cnt in bitstring_samples:
                if abs(e - min_energy) < 1e-10:
                    self.best_bitstring = bs.copy()
                    break

        self.energy_history.append(cvar)
        self.iteration += 1

        return cvar

    def _objective_wrapper(
        self, parameters: np.ndarray, qubo: np.ndarray
    ) -> Tuple[float, np.ndarray]:
        """
        Wrapper for scipy optimizer.

        Args:
            parameters: Ansatz parameters.
            qubo: QUBO matrix.

        Returns:
            Tuple of (energy, gradient approximation).
        """
        energy = self._evaluate_cvar_objective(parameters, qubo)
        # Numerical gradient using finite differences
        gradient = np.zeros_like(parameters)
        eps = 1e-5
        for i in range(len(parameters)):
            params_plus = parameters.copy()
            params_plus[i] += eps
            e_plus = self._evaluate_cvar_objective(params_plus, qubo)
            gradient[i] = (e_plus - energy) / eps

        return energy, gradient

    def solve(self, qubo: np.ndarray) -> Tuple[np.ndarray, float]:
        """
        Solve QUBO using CVaR-VQE.

        Args:
            qubo: QUBO matrix.

        Returns:
            Tuple of (best_bitstring, best_energy).
        """
        # Lazy import scipy to avoid import-time issues
        from scipy.optimize import minimize

        if qubo.shape[0] != self.num_qubits:
            raise ValueError(
                f"QUBO size {qubo.shape[0]} doesn't match num_qubits {self.num_qubits}"
            )

        # Reset history
        self.iteration = 0
        self.energy_history = []
        self.best_energy = float("inf")
        self.best_parameters = None
        self.best_bitstring = None

        # Initial parameters
        initial_params = self.ansatz.initial_parameters.copy()

        # Optimize using scipy
        def objective(params):
            return self._evaluate_cvar_objective(params, qubo)

        result = minimize(
            objective,
            initial_params,
            method="COBYLA",
            options={
                "maxiter": self.max_iterations,
                "tol": 1e-5,
            },
        )

        # Ensure we have a valid solution
        if self.best_bitstring is None or self.best_parameters is None:
            # Fallback: run final circuit with best parameters
            if self.best_parameters is None:
                self.best_parameters = initial_params

            circuit = self._create_measurement_circuit(self.best_parameters)
            job = self.simulator.run(circuit, shots=self.shots, seed_simulator=self.seed)
            result_obj = job.result()
            counts = result_obj.get_counts()

            # Find bitstring with lowest energy
            min_energy = float("inf")
            best_bs = None
            for bitstring_str, count in counts.items():
                bitstring = np.array([int(b) for b in reversed(bitstring_str)])
                energy = float(bitstring @ qubo @ bitstring)
                if energy < min_energy:
                    min_energy = energy
                    best_bs = bitstring

            if best_bs is None:
                # Last-resort fallback
                best_bs = np.zeros(self.num_qubits, dtype=int)
                min_energy = 0.0

            self.best_bitstring = best_bs
            self.best_energy = min_energy

        return self.best_bitstring.astype(int), float(self.best_energy)

    def get_energy_history(self) -> List[float]:
        """
        Get the optimization energy history.

        Returns:
            List of CVaR energies at each iteration.
        """
        return self.energy_history

    def get_best_parameters(self) -> Optional[np.ndarray]:
        """
        Get the best parameters found during optimization.

        Returns:
            Best parameters array or None if not yet optimized.
        """
        return self.best_parameters

    def compute_cvar(self, energies: List[float]) -> float:
        """
        Compute CVaR(α) from a list of energies.

        Args:
            energies: List of sampled energies.

        Returns:
            CVaR value.
        """
        sorted_energies = sorted(energies)
        k = max(1, int(np.ceil(self.alpha * len(sorted_energies))))
        return float(np.mean(sorted_energies[:k]))