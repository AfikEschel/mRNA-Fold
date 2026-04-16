"""
Base pair validation for mRNA secondary structure prediction.

This module implements canonical and wobble base pair rules
as described in the IBM–Moderna paper and standard RNA folding literature.
"""

from typing import Set, Tuple


# Canonical Watson-Crick base pairs
CANONICAL_PAIRS: Set[Tuple[str, str]] = {
    ("A", "U"), ("U", "A"),
    ("C", "G"), ("G", "C")
}

# Wobble pairs (G-U base pairs)
WOBBLE_PAIRS: Set[Tuple[str, str]] = {
    ("G", "U"), ("U", "G")
}

# All valid base pairs (canonical + wobble)
VALID_PAIRS: Set[Tuple[str, str]] = CANONICAL_PAIRS | WOBBLE_PAIRS


def is_canonical_pair(nuc1: str, nuc2: str) -> bool:
    """
    Check if two nucleotides form a canonical Watson-Crick base pair.
    
    Canonical pairs: A-U, U-A, C-G, G-C
    
    Args:
        nuc1: First nucleotide (A, U, C, G)
        nuc2: Second nucleotide (A, U, C, G)
    
    Returns:
        True if nucleotides form a canonical pair, False otherwise
    """
    return (nuc1, nuc2) in CANONICAL_PAIRS


def is_wobble_pair(nuc1: str, nuc2: str) -> bool:
    """
    Check if two nucleotides form a wobble (G-U) base pair.
    
    Wobble pairs: G-U, U-G
    
    Args:
        nuc1: First nucleotide (A, U, C, G)
        nuc2: Second nucleotide (A, U, C, G)
    
    Returns:
        True if nucleotides form a wobble pair, False otherwise
    """
    return (nuc1, nuc2) in WOBBLE_PAIRS


def is_valid_pair(nuc1: str, nuc2: str) -> bool:
    """
    Check if two nucleotides can form any valid base pair.
    
    Valid pairs include both canonical and wobble pairs.
    
    Args:
        nuc1: First nucleotide (A, U, C, G)
        nuc2: Second nucleotide (A, U, C, G)
    
    Returns:
        True if nucleotides can form a valid pair, False otherwise
    """
    return (nuc1, nuc2) in VALID_PAIRS


def can_pair(i: int, j: int, sequence: str) -> bool:
    """
    Check if positions i and j in a sequence can form a base pair.
    
    This is a convenience function that extracts nucleotides from a sequence
    and checks if they can form a valid base pair.
    
    Args:
        i: First position (0-based)
        j: Second position (0-based)
        sequence: RNA sequence string
    
    Returns:
        True if positions i,j can form a valid base pair
    
    Raises:
        IndexError: if i or j are out of bounds
        ValueError: if sequence contains invalid nucleotides
    """
    if i < 0 or j < 0 or i >= len(sequence) or j >= len(sequence):
        raise IndexError(f"Position out of bounds: i={i}, j={j}, len={len(sequence)}")
    
    nuc1 = sequence[i]
    nuc2 = sequence[j]
    
    # Validate nucleotides
    valid_nucs = set("AUGC")
    if nuc1 not in valid_nucs or nuc2 not in valid_nucs:
        raise ValueError(f"Invalid nucleotides: '{nuc1}' or '{nuc2}' (must be A, U, C, G)")
    
    return is_valid_pair(nuc1, nuc2)


def get_pair_type(nuc1: str, nuc2: str) -> str:
    """
    Get the type of base pair formed by two nucleotides.
    
    Args:
        nuc1: First nucleotide
        nuc2: Second nucleotide
    
    Returns:
        "canonical" if Watson-Crick pair
        "wobble" if G-U pair
        "invalid" if no valid pair
    """
    if is_canonical_pair(nuc1, nuc2):
        return "canonical"
    elif is_wobble_pair(nuc1, nuc2):
        return "wobble"
    else:
        return "invalid"
