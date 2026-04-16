"""
Tests for QUBO formulation module.
"""

import pytest
import numpy as np
from mrnafold.quartets import Quartet, generate_quartets
from mrnafold.qubo import build_quartet_qubo, qubo_to_ising, solve_qubo_brute_force, quartets_to_structure, structure_to_dot_bracket


class TestQUBO:
    """Test QUBO matrix construction."""

    def test_build_quartet_qubo_simple(self):
        """Test QUBO construction with simple quartet set."""
        # Create a simple set of quartets
        quartets = [
            Quartet(i=0, j=3, pair_type="canonical"),
            Quartet(i=1, j=4, pair_type="canonical"),
            Quartet(i=2, j=5, pair_type="canonical")
        ]

        qubo = build_quartet_qubo(quartets)

        # Check matrix shape
        assert qubo.shape == (3, 3)

        # Check that conflicting quartets have positive coefficients
        # Quartet(0,3) and Quartet(1,4) conflict (share positions 1,3)
        assert qubo[0, 1] > 0
        assert qubo[1, 0] > 0

        # Quartet(1,4) and Quartet(2,5) conflict (share positions 2,4)
        assert qubo[1, 2] > 0
        assert qubo[2, 1] > 0

    def test_build_quartet_qubo_from_sequence(self):
        """Test QUBO construction from real sequence."""
        seq = 'AUCUGCAUGGCCAAGAGGGUUA'
        quartets = generate_quartets(seq)
        qubo = build_quartet_qubo(quartets)

        assert qubo.shape == (len(quartets), len(quartets))

        # Matrix should be symmetric
        assert np.allclose(qubo, qubo.T)

        # Diagonal should be zero (no self-interactions)
        assert np.all(np.diag(qubo) == 0)

    def test_qubo_parameters(self):
        """Test QUBO with custom parameters."""
        quartets = [
            Quartet(i=0, j=3, pair_type="canonical"),
            Quartet(i=1, j=4, pair_type="canonical")
        ]

        # Test with different penalties
        qubo1 = build_quartet_qubo(quartets, conflict_penalty=2.0)
        qubo2 = build_quartet_qubo(quartets, conflict_penalty=1.0)

        # Higher penalty should give higher coefficients
        assert qubo1[0, 1] == 2.0 * qubo2[0, 1]


class TestQUBOToIsing:
    """Test QUBO to Ising conversion."""

    def test_qubo_to_ising_conversion(self):
        """Test basic QUBO to Ising conversion."""
        # Simple 2x2 QUBO
        qubo = np.array([
            [0, 1],
            [1, 0]
        ])

        ising, constant = qubo_to_ising(qubo)

        # Expected Ising matrix
        expected_ising = np.array([
            [0.5, 0.5],
            [0.5, 0.5]
        ])

        assert np.allclose(ising, expected_ising)
        assert constant == 0.5

    def test_ising_matrix_properties(self):
        """Test that Ising matrix has correct properties."""
        seq = 'AUGC'
        quartets = generate_quartets(seq)
        qubo = build_quartet_qubo(quartets)
        ising, constant = qubo_to_ising(qubo)

        # Ising matrix should be symmetric
        assert np.allclose(ising, ising.T)

        # Constant should be finite
        assert np.isfinite(constant)


class TestBruteForceSolver:
    """Test brute force QUBO solver."""

    def test_solve_small_qubo(self):
        """Test brute force solver on small instance."""
        # Simple symmetric QUBO: minimize -2*x1*x2 when both variables are 1
        qubo = np.array([
            [0, -1],
            [-1, 0]
        ])

        solution, energy = solve_qubo_brute_force(qubo)

        # Should find the optimal solution with the symmetric negative off-diagonal coefficient
        assert energy == -2
        assert np.array_equal(solution, [1, 1])

    def test_solve_empty_qubo(self):
        """Test solver with no quartets."""
        qubo = np.zeros((0, 0))
        solution, energy = solve_qubo_brute_force(qubo)

        assert len(solution) == 0
        assert energy == 0

    def test_brute_force_size_limit(self):
        """Test that large QUBOs are rejected."""
        large_qubo = np.zeros((25, 25))

        with pytest.raises(ValueError, match="Brute force only works for n <= 20"):
            solve_qubo_brute_force(large_qubo)


class TestStructureConversion:
    """Test conversion from quartets to structure."""

    def test_quartets_to_structure(self):
        """Test conversion of selected quartets to base pairs."""
        quartets = [
            Quartet(i=0, j=5, pair_type="canonical"),  # spans 0-5, inner 1-4
        ]

        base_pairs = quartets_to_structure(quartets, 6)

        expected = [(0, 5), (1, 4)]
        assert base_pairs == expected

    def test_multiple_quartets_to_structure(self):
        """Test conversion with multiple non-conflicting quartets."""
        quartets = [
            Quartet(i=0, j=3, pair_type="canonical"),  # (0,3) and (1,2)
            Quartet(i=4, j=7, pair_type="canonical"),  # (4,7) and (5,6)
        ]

        base_pairs = quartets_to_structure(quartets, 8)

        expected = [(0, 3), (1, 2), (4, 7), (5, 6)]
        assert base_pairs == expected

    def test_structure_to_dot_bracket(self):
        """Test conversion to dot-bracket notation."""
        base_pairs = [(0, 5), (1, 4)]
        dot_bracket = structure_to_dot_bracket(base_pairs, 6)

        assert dot_bracket == "((..))"

    def test_dot_bracket_unpaired(self):
        """Test dot-bracket with unpaired bases."""
        base_pairs = [(1, 4)]
        dot_bracket = structure_to_dot_bracket(base_pairs, 6)

        assert dot_bracket == ".(..)."

    def test_empty_structure(self):
        """Test conversion with no base pairs."""
        base_pairs = []
        dot_bracket = structure_to_dot_bracket(base_pairs, 4)

        assert dot_bracket == "...."


class TestIntegration:
    """Integration tests for QUBO pipeline."""

    def test_full_pipeline_small(self):
        """Test complete pipeline on small sequence."""
        seq = 'AUGC'  # Length 4
        quartets = generate_quartets(seq)

        # Build QUBO
        qubo = build_quartet_qubo(quartets)

        # Solve (brute force for small instance)
        solution, energy = solve_qubo_brute_force(qubo)

        # Convert to structure
        selected_quartets = [q for q, selected in zip(quartets, solution) if selected]
        base_pairs = quartets_to_structure(selected_quartets, len(seq))
        dot_bracket = structure_to_dot_bracket(base_pairs, len(seq))

        # Should produce valid dot-bracket string
        assert len(dot_bracket) == len(seq)
        assert all(c in '.()' for c in dot_bracket)

    def test_full_pipeline_real_sequence(self):
        """Test pipeline on real sequence (without solving large QUBO)."""
        seq = 'AUCUGCAUGGCCAAGAGGGUUA'
        quartets = generate_quartets(seq)

        # Build QUBO
        qubo = build_quartet_qubo(quartets)

        # Just check that QUBO has reasonable properties
        assert qubo.shape[0] == len(quartets)
        assert np.allclose(qubo, qubo.T)  # Symmetric
        assert np.all(np.diag(qubo) == 0)  # No self-terms