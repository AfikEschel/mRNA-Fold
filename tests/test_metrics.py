"""
Unit tests for structure evaluation metrics.
"""

import pytest
from mrnafold.metrics import (structure_to_base_pairs, base_pairs_to_structure,
                             calculate_f1_score, calculate_sensitivity, calculate_ppv,
                             evaluate_structure_prediction)


class TestStructureConversion:
    """Test conversion between dot-bracket and base pairs."""

    def test_structure_to_base_pairs_simple(self):
        """Test conversion of simple structure."""
        structure = "((..))"
        pairs = structure_to_base_pairs(structure)
        expected = [(0, 5), (1, 4)]  # length 6
        assert pairs == expected

    def test_structure_to_base_pairs_no_pairs(self):
        """Test structure with no base pairs."""
        structure = "....."
        pairs = structure_to_base_pairs(structure)
        assert pairs == []

    def test_base_pairs_to_structure(self):
        """Test conversion from base pairs to dot-bracket."""
        pairs = [(0, 5), (1, 4)]
        structure = base_pairs_to_structure(pairs, 6)
        assert structure == "((..))"

    def test_round_trip_conversion(self):
        """Test that conversion is reversible."""
        original = "((..))"
        pairs = structure_to_base_pairs(original)
        reconstructed = base_pairs_to_structure(pairs, len(original))
        assert reconstructed == original


class TestMetrics:
    """Test evaluation metrics."""

    def test_perfect_prediction(self):
        """Test metrics for perfect prediction."""
        ref = "((..))"
        pred = "((..))"
        metrics = evaluate_structure_prediction(pred, ref)

        assert metrics['f1_score'] == 1.0
        assert metrics['sensitivity'] == 1.0
        assert metrics['ppv'] == 1.0
        assert metrics['predicted_pairs'] == 2
        assert metrics['reference_pairs'] == 2

    def test_no_prediction(self):
        """Test metrics when no pairs are predicted."""
        ref = "((..))"  # length 6
        pred = "......"  # length 6, no pairs
        metrics = evaluate_structure_prediction(pred, ref)

        assert metrics['f1_score'] == 0.0
        assert metrics['sensitivity'] == 0.0
        assert metrics['ppv'] == 1.0  # No false positives
        assert metrics['predicted_pairs'] == 0
        assert metrics['reference_pairs'] == 2

    def test_no_reference_pairs(self):
        """Test metrics when reference has no pairs."""
        ref = "....."
        pred = "....."
        metrics = evaluate_structure_prediction(pred, ref)

        assert metrics['f1_score'] == 1.0  # Perfect match
        assert metrics['sensitivity'] == 1.0
        assert metrics['ppv'] == 1.0
        assert metrics['predicted_pairs'] == 0
        assert metrics['reference_pairs'] == 0

    def test_partial_overlap(self):
        """Test metrics with partial overlap."""
        ref = "((..))"  # length 6, pairs (0,5),(1,4)
        pred = "(....)"  # length 6, pairs (0,5) only
        metrics = evaluate_structure_prediction(pred, ref)

        # pred has (0,5), ref has (0,5) and (1,4)
        # TP=1, FP=0, FN=1
        expected_f1 = 2 * (1/1) * (1/2) / (1/1 + 1/2)  # 2/3
        assert abs(metrics['f1_score'] - expected_f1) < 1e-6
        assert metrics['sensitivity'] == 0.5  # 1/2
        assert metrics['ppv'] == 1.0  # 1/1

    def test_calculate_f1_score(self):
        """Test F1 score calculation directly."""
        pred = [(0, 4)]
        ref = [(0, 4), (1, 3)]
        f1 = calculate_f1_score(pred, ref)
        expected = 2/3
        assert abs(f1 - expected) < 1e-6

    def test_calculate_sensitivity(self):
        """Test sensitivity calculation."""
        pred = [(0, 4)]
        ref = [(0, 4), (1, 3)]
        sens = calculate_sensitivity(pred, ref)
        assert sens == 0.5

    def test_calculate_ppv(self):
        """Test PPV calculation."""
        pred = [(0, 4), (2, 5)]  # One correct, one incorrect
        ref = [(0, 4), (1, 3)]
        ppv = calculate_ppv(pred, ref)
        assert ppv == 0.5