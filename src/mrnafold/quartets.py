"""
Quartet representation and generation for mRNA secondary structure prediction.

A quartet is two nested consecutive base pairs: (i, j) and (i+1, j-1).
This encoding is central to the IBM–Moderna paper's QUBO formulation.
"""

from dataclasses import dataclass
from typing import List, Set, Tuple
from mrnafold.pairing import is_valid_pair, get_pair_type


@dataclass(frozen=True)
class Quartet:
    """
    Represents a quartet: two nested consecutive base pairs.
    
    A quartet consists of positions (i, j) and (i+1, j-1) where:
    - i < i+1 < j-1 < j
    - Both (i,j) and (i+1,j-1) are valid base pairs
    
    Attributes:
        i: First position of the outer pair
        j: Second position of the outer pair
        pair_type: Type of the outer pair ("canonical" or "wobble")
    """
    i: int
    j: int
    pair_type: str
    
    def __post_init__(self):
        """Validate quartet properties."""
        if not (self.i < self.i + 1 < self.j - 1 < self.j):
            raise ValueError(f"Invalid quartet positions: i={self.i}, j={self.j}")
        if self.j - self.i < 3:
            raise ValueError(f"Quartet too small: j-i={self.j - self.i} < 3")
        if self.pair_type not in ["canonical", "wobble"]:
            raise ValueError(f"Invalid pair type: {self.pair_type}")
    
    @property
    def inner_i(self) -> int:
        """Position i of the inner pair."""
        return self.i + 1
    
    @property
    def inner_j(self) -> int:
        """Position j of the inner pair."""
        return self.j - 1
    
    def conflicts_with(self, other: 'Quartet') -> bool:
        """
        Check if this quartet conflicts with another quartet.
        
        Quartets conflict if they share any base pair positions.
        """
        self_positions = {self.i, self.j, self.inner_i, self.inner_j}
        other_positions = {other.i, other.j, other.inner_i, other.inner_j}
        return bool(self_positions & other_positions)
    
    def __str__(self) -> str:
        return f"Quartet({self.i},{self.j})[{self.pair_type}]"


def generate_quartets(sequence: str) -> List[Quartet]:
    """
    Generate all valid quartets from a sequence.
    
    A quartet is valid if:
    1. Positions i < i+1 < j-1 < j
    2. Both (i,j) and (i+1,j-1) are valid base pairs
    3. Minimum quartet size (j - i >= 3)
    
    Args:
        sequence: RNA sequence string
    
    Returns:
        List of all valid Quartet objects
    """
    quartets = []
    n = len(sequence)
    
    # Minimum quartet requires j - i >= 3
    for i in range(n - 3):  # i from 0 to n-4
        for j in range(i + 3, n):  # j from i+3 to n-1
            # Check if both pairs are valid
            if is_valid_pair(sequence[i], sequence[j]) and is_valid_pair(sequence[i+1], sequence[j-1]):
                # Determine pair type for the outer pair
                pair_type = get_pair_type(sequence[i], sequence[j])
                
                quartet = Quartet(i=i, j=j, pair_type=pair_type)
                quartets.append(quartet)
    
    return quartets


def find_conflicting_quartets(quartets: List[Quartet]) -> Set[Tuple[Quartet, Quartet]]:
    """
    Find all pairs of conflicting quartets.
    
    Args:
        quartets: List of Quartet objects
    
    Returns:
        Set of (quartet1, quartet2) tuples where quartet1 < quartet2
        to avoid duplicates
    """
    conflicts = set()
    n = len(quartets)
    
    for i in range(n):
        for j in range(i + 1, n):
            if quartets[i].conflicts_with(quartets[j]):
                # Store in sorted order to avoid duplicates
                conflicts.add((quartets[i], quartets[j]))
    
    return conflicts


def find_stackable_quartets(quartets: List[Quartet]) -> dict[Quartet, List[Quartet]]:
    """
    Find quartets that can stack with each other.
    
    Two quartets can stack if they are adjacent and compatible.
    For quartet Q1(i,j) and Q2(k,l), they can stack if:
    - Q2.i == Q1.i + 2 and Q2.j == Q1.j - 2 (nested stacking)
    
    Args:
        quartets: List of Quartet objects
    
    Returns:
        Dictionary mapping each quartet to list of quartets it can stack with
    """
    stackable = {q: [] for q in quartets}
    
    for q1 in quartets:
        for q2 in quartets:
            if q1 == q2:
                continue
            
            # Check if q2 can stack inside q1
            if q2.i == q1.i + 2 and q2.j == q1.j - 2:
                stackable[q1].append(q2)
    
    return stackable
