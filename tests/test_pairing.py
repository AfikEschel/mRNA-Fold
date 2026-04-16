"""
Unit tests for base pair validation logic.
"""

import pytest
from mrnafold.pairing import (
    is_canonical_pair,
    is_wobble_pair,
    is_valid_pair,
    can_pair,
    get_pair_type,
    CANONICAL_PAIRS,
    WOBBLE_PAIRS,
    VALID_PAIRS
)


class TestCanonicalPairs:
    """Test canonical Watson-Crick base pair detection."""
    
    def test_canonical_pairs_true(self):
        """Test that canonical pairs are correctly identified."""
        assert is_canonical_pair("A", "U") is True
        assert is_canonical_pair("U", "A") is True
        assert is_canonical_pair("C", "G") is True
        assert is_canonical_pair("G", "C") is True
    
    def test_canonical_pairs_false(self):
        """Test that non-canonical pairs are correctly rejected."""
        assert is_canonical_pair("A", "A") is False
        assert is_canonical_pair("G", "U") is False  # This is wobble, not canonical
        assert is_canonical_pair("C", "U") is False
        assert is_canonical_pair("A", "G") is False


class TestWobblePairs:
    """Test wobble (G-U) base pair detection."""
    
    def test_wobble_pairs_true(self):
        """Test that wobble pairs are correctly identified."""
        assert is_wobble_pair("G", "U") is True
        assert is_wobble_pair("U", "G") is True
    
    def test_wobble_pairs_false(self):
        """Test that non-wobble pairs are correctly rejected."""
        assert is_wobble_pair("A", "U") is False  # This is canonical
        assert is_wobble_pair("C", "G") is False  # This is canonical
        assert is_wobble_pair("A", "A") is False
        assert is_wobble_pair("G", "G") is False


class TestValidPairs:
    """Test general valid pair detection (canonical + wobble)."""
    
    def test_valid_pairs_canonical(self):
        """Test that canonical pairs are valid."""
        assert is_valid_pair("A", "U") is True
        assert is_valid_pair("U", "A") is True
        assert is_valid_pair("C", "G") is True
        assert is_valid_pair("G", "C") is True
    
    def test_valid_pairs_wobble(self):
        """Test that wobble pairs are valid."""
        assert is_valid_pair("G", "U") is True
        assert is_valid_pair("U", "G") is True
    
    def test_invalid_pairs(self):
        """Test that invalid pairs are correctly rejected."""
        assert is_valid_pair("A", "A") is False
        assert is_valid_pair("U", "U") is False
        assert is_valid_pair("C", "C") is False
        assert is_valid_pair("G", "G") is False
        assert is_valid_pair("A", "C") is False
        assert is_valid_pair("A", "G") is False
        assert is_valid_pair("C", "U") is False


class TestPairType:
    """Test pair type classification."""
    
    def test_canonical_type(self):
        """Test canonical pair type detection."""
        assert get_pair_type("A", "U") == "canonical"
        assert get_pair_type("U", "A") == "canonical"
        assert get_pair_type("C", "G") == "canonical"
        assert get_pair_type("G", "C") == "canonical"
    
    def test_wobble_type(self):
        """Test wobble pair type detection."""
        assert get_pair_type("G", "U") == "wobble"
        assert get_pair_type("U", "G") == "wobble"
    
    def test_invalid_type(self):
        """Test invalid pair type detection."""
        assert get_pair_type("A", "A") == "invalid"
        assert get_pair_type("A", "C") == "invalid"
        assert get_pair_type("G", "A") == "invalid"


class TestSequencePairing:
    """Test pairing validation within sequences."""
    
    def test_can_pair_valid(self):
        """Test valid pairing in sequences."""
        seq = "AUCUGCAUGGCCAAGAGGGUUA"
        
        # Test actual pairs from the structure .......(((((.....)))))
        # Pairs are: 7↔21, 8↔20, 9↔19, 10↔18, 11↔17
        assert can_pair(7, 21, seq) is True   # U-A (canonical)
        assert can_pair(8, 20, seq) is True   # G-U (wobble)
        assert can_pair(9, 19, seq) is True   # G-U (wobble)
        assert can_pair(10, 18, seq) is True  # C-G (canonical)
        assert can_pair(11, 17, seq) is True  # C-G (canonical)
    
    def test_can_pair_invalid(self):
        """Test invalid pairing in sequences."""
        seq = "AUCUGCAUGGCCAAGAGGGUUA"
        
        # Test some invalid pairs (nucleotides that cannot form base pairs)
        assert can_pair(0, 2, seq) is False   # A-C
        assert can_pair(1, 3, seq) is False   # U-U
        assert can_pair(0, 4, seq) is False   # A-G
        assert can_pair(2, 3, seq) is False   # C-U
        assert can_pair(4, 6, seq) is False   # G-A
    
    def test_can_pair_bounds(self):
        """Test boundary conditions."""
        seq = "AUCG"
        
        # Mix of valid and invalid pairs
        assert can_pair(0, 3, seq) is False  # A-G (invalid pair)
        assert can_pair(1, 3, seq) is True   # U-G (valid wobble pair)
        
        # Invalid indices should raise IndexError
        with pytest.raises(IndexError):
            can_pair(-1, 0, seq)
        with pytest.raises(IndexError):
            can_pair(0, 4, seq)
        with pytest.raises(IndexError):
            can_pair(4, 0, seq)
        
        # Invalid indices
        with pytest.raises(IndexError):
            can_pair(-1, 0, seq)
        with pytest.raises(IndexError):
            can_pair(0, 4, seq)
        with pytest.raises(IndexError):
            can_pair(4, 0, seq)
    
    def test_invalid_nucleotides(self):
        """Test error handling for invalid nucleotides."""
        seq = "AUCX"  # X is invalid
        
        with pytest.raises(ValueError):
            can_pair(0, 3, seq)  # A-X
        with pytest.raises(ValueError):
            can_pair(3, 0, seq)  # X-A


class TestConstants:
    """Test that constants are properly defined."""
    
    def test_canonical_pairs_constant(self):
        """Test CANONICAL_PAIRS constant."""
        expected = {("A", "U"), ("U", "A"), ("C", "G"), ("G", "C")}
        assert CANONICAL_PAIRS == expected
    
    def test_wobble_pairs_constant(self):
        """Test WOBBLE_PAIRS constant."""
        expected = {("G", "U"), ("U", "G")}
        assert WOBBLE_PAIRS == expected
    
    def test_valid_pairs_constant(self):
        """Test VALID_PAIRS constant (union of canonical and wobble)."""
        expected = CANONICAL_PAIRS | WOBBLE_PAIRS
        assert VALID_PAIRS == expected
        assert len(VALID_PAIRS) == 6  # 4 canonical + 2 wobble