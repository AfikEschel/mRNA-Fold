"""
Unit tests for quartet generation and logic.
"""

import pytest
from mrnafold.quartets import Quartet, generate_quartets, find_conflicting_quartets, find_stackable_quartets


class TestQuartet:
    """Test Quartet dataclass."""
    
    def test_valid_quartet(self):
        """Test creating a valid quartet."""
        q = Quartet(i=0, j=5, pair_type="canonical")
        assert q.i == 0
        assert q.j == 5
        assert q.inner_i == 1
        assert q.inner_j == 4
        assert q.pair_type == "canonical"
    
    def test_invalid_positions(self):
        """Test that invalid positions raise ValueError."""
        with pytest.raises(ValueError):
            Quartet(i=0, j=2, pair_type="canonical")  # j - i < 3
        
        with pytest.raises(ValueError):
            Quartet(i=5, j=0, pair_type="canonical")  # i > j
    
    def test_invalid_pair_type(self):
        """Test that invalid pair type raises ValueError."""
        with pytest.raises(ValueError):
            Quartet(i=0, j=5, pair_type="invalid")
    
    def test_quartet_string(self):
        """Test quartet string representation."""
        q = Quartet(i=1, j=6, pair_type="wobble")
        assert str(q) == "Quartet(1,6)[wobble]"
    
    def test_conflicts_with(self):
        """Test quartet conflict detection."""
        q1 = Quartet(i=0, j=5, pair_type="canonical")
        q2 = Quartet(i=1, j=4, pair_type="canonical")  # Shares positions 1,4
        q3 = Quartet(i=6, j=9, pair_type="canonical")  # No shared positions
        
        assert q1.conflicts_with(q2) is True
        assert q1.conflicts_with(q3) is False
        assert q2.conflicts_with(q3) is False


class TestQuartetGeneration:
    """Test quartet generation from sequences."""
    
    def test_generate_quartets_simple(self):
        """Test quartet generation on a simple sequence."""
        # Sequence: A U C G A U
        # Positions: 0 1 2 3 4 5
        # Valid pairs: (0,5):A-U, (1,4):U-A, (2,3):C-G
        # Valid quartets: (0,5) with inner (1,4), and (1,4) with inner (2,3)
        seq = "AUCGAU"
        quartets = generate_quartets(seq)
        
        # Should find two quartets: i=0,j=5 and i=1,j=4
        assert len(quartets) == 2
        quartets_sorted = sorted(quartets, key=lambda q: (q.i, q.j))
        assert quartets_sorted[0].i == 0 and quartets_sorted[0].j == 5
        assert quartets_sorted[1].i == 1 and quartets_sorted[1].j == 4
    
    def test_generate_quartets_no_quartets(self):
        """Test sequence with no valid quartets."""
        seq = "AAAA"  # No valid pairs
        quartets = generate_quartets(seq)
        assert len(quartets) == 0
    
    def test_generate_quartets_from_dataset(self):
        """Test quartet generation on actual dataset sequence."""
        # From the dataset: AUCUGCAUGGCCAAGAGGGUUA
        seq = "AUCUGCAUGGCCAAGAGGGUUA"
        quartets = generate_quartets(seq)
        
        # Should find several quartets
        assert len(quartets) > 0
        
        # Check that all quartets are valid
        for q in quartets:
            assert q.i < q.inner_i < q.inner_j < q.j
            assert q.j - q.i >= 3  # Minimum quartet span
    
    def test_quartet_properties(self):
        """Test quartet properties on generated quartets."""
        seq = "AUCUGCAUGGCCAAGAGGGUUA"
        quartets = generate_quartets(seq)
        
        for q in quartets:
            # Check inner positions
            assert q.inner_i == q.i + 1
            assert q.inner_j == q.j - 1
            
            # Check pair type is valid
            assert q.pair_type in ["canonical", "wobble"]


class TestQuartetConflicts:
    """Test quartet conflict detection."""
    
    def test_no_conflicts(self):
        """Test quartets with no conflicts."""
        q1 = Quartet(i=0, j=5, pair_type="canonical")
        q2 = Quartet(i=6, j=9, pair_type="canonical")
        
        conflicts = find_conflicting_quartets([q1, q2])
        assert len(conflicts) == 0
    
    def test_with_conflicts(self):
        """Test quartets with conflicts."""
        q1 = Quartet(i=0, j=7, pair_type="canonical")
        q2 = Quartet(i=1, j=6, pair_type="canonical")  # Shares positions with q1
        q3 = Quartet(i=8, j=11, pair_type="canonical")  # No conflict
        
        conflicts = find_conflicting_quartets([q1, q2, q3])
        assert len(conflicts) == 1
        
        # Check the conflict pair
        conflict_pair = list(conflicts)[0]
        assert q1 in conflict_pair and q2 in conflict_pair
    
    def test_multiple_conflicts(self):
        """Test multiple conflicting quartets."""
        q1 = Quartet(i=0, j=5, pair_type="canonical")
        q2 = Quartet(i=1, j=4, pair_type="canonical")  # Conflicts with q1
        q3 = Quartet(i=2, j=5, pair_type="canonical")  # Conflicts with q1 and q2
        
        conflicts = find_conflicting_quartets([q1, q2, q3])
        assert len(conflicts) == 3  # All pairs conflict


class TestQuartetStacking:
    """Test quartet stacking detection."""
    
    def test_stackable_quartets(self):
        """Test finding stackable quartets."""
        q1 = Quartet(i=0, j=7, pair_type="canonical")
        q2 = Quartet(i=2, j=5, pair_type="canonical")  # Can stack inside q1
        q3 = Quartet(i=8, j=11, pair_type="canonical")  # Cannot stack
        
        stackable = find_stackable_quartets([q1, q2, q3])
        
        assert q2 in stackable[q1]
        assert len(stackable[q1]) == 1
        assert len(stackable[q2]) == 0  # q2 cannot stack with anything
        assert len(stackable[q3]) == 0
    
    def test_no_stacking(self):
        """Test quartets that cannot stack."""
        q1 = Quartet(i=0, j=5, pair_type="canonical")
        q2 = Quartet(i=6, j=9, pair_type="canonical")  # Too far apart
        
        stackable = find_stackable_quartets([q1, q2])
        
        assert len(stackable[q1]) == 0
        assert len(stackable[q2]) == 0


class TestIntegration:
    """Integration tests with real sequences."""
    
    def test_full_pipeline(self):
        """Test the full quartet generation pipeline."""
        seq = "AUCUGCAUGGCCAAGAGGGUUA"
        
        # Generate quartets
        quartets = generate_quartets(seq)
        assert len(quartets) > 0
        
        # Find conflicts
        conflicts = find_conflicting_quartets(quartets)
        # Some conflicts expected in real sequences
        
        # Find stackable
        stackable = find_stackable_quartets(quartets)
        # Some stacking expected
        
        # Verify all data structures are consistent
        assert isinstance(conflicts, set)
        assert isinstance(stackable, dict)
        assert all(isinstance(q, Quartet) for q in quartets)
        assert all(isinstance(stack_list, list) for stack_list in stackable.values())
