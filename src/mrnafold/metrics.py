"""
Metrics and evaluation for mRNA secondary structure prediction.

This module provides functions to compare predicted structures against
reference structures and evaluate prediction quality.
"""

from typing import List, Tuple
import numpy as np


def structure_to_base_pairs(dot_bracket: str) -> List[Tuple[int, int]]:
    """
    Convert dot-bracket notation to list of base pairs.

    Args:
        dot_bracket: Structure in dot-bracket notation

    Returns:
        List of (i,j) base pairs (0-based indexing)
    """
    base_pairs = []
    stack = []

    for i, char in enumerate(dot_bracket):
        if char == '(':
            stack.append(i)
        elif char == ')':
            if stack:
                j = stack.pop()
                base_pairs.append((j, i))

    return sorted(base_pairs)


def base_pairs_to_structure(base_pairs: List[Tuple[int, int]],
                           sequence_length: int) -> str:
    """
    Convert base pairs back to dot-bracket notation.

    Args:
        base_pairs: List of (i,j) base pairs
        sequence_length: Length of sequence

    Returns:
        Dot-bracket string
    """
    structure = ['.'] * sequence_length

    for i, j in base_pairs:
        if structure[i] == '.' and structure[j] == '.':
            structure[i] = '('
            structure[j] = ')'

    return ''.join(structure)


def calculate_f1_score(predicted_pairs: List[Tuple[int, int]],
                      reference_pairs: List[Tuple[int, int]]) -> float:
    """
    Calculate F1 score for base pair prediction.

    F1 = 2 * precision * recall / (precision + recall)

    Args:
        predicted_pairs: Predicted base pairs
        reference_pairs: Reference base pairs

    Returns:
        F1 score (0.0 to 1.0)
    """
    pred_set = set(predicted_pairs)
    ref_set = set(reference_pairs)

    true_positives = len(pred_set & ref_set)
    false_positives = len(pred_set - ref_set)
    false_negatives = len(ref_set - pred_set)

    # Special case: both empty sets (perfect match)
    if len(pred_set) == 0 and len(ref_set) == 0:
        return 1.0

    if true_positives == 0:
        return 0.0

    precision = true_positives / (true_positives + false_positives)
    recall = true_positives / (true_positives + false_negatives)

    f1 = 2 * precision * recall / (precision + recall)
    return f1


def calculate_sensitivity(predicted_pairs: List[Tuple[int, int]],
                         reference_pairs: List[Tuple[int, int]]) -> float:
    """
    Calculate sensitivity (recall) for base pair prediction.

    Sensitivity = TP / (TP + FN)

    Args:
        predicted_pairs: Predicted base pairs
        reference_pairs: Reference base pairs

    Returns:
        Sensitivity score (0.0 to 1.0)
    """
    pred_set = set(predicted_pairs)
    ref_set = set(reference_pairs)

    true_positives = len(pred_set & ref_set)
    false_negatives = len(ref_set - pred_set)

    if len(ref_set) == 0:
        return 1.0 if len(pred_set) == 0 else 0.0

    return true_positives / (true_positives + false_negatives)


def calculate_ppv(predicted_pairs: List[Tuple[int, int]],
                 reference_pairs: List[Tuple[int, int]]) -> float:
    """
    Calculate positive predictive value (precision) for base pair prediction.

    PPV = TP / (TP + FP)

    Args:
        predicted_pairs: Predicted base pairs
        reference_pairs: Reference base pairs

    Returns:
        PPV score (0.0 to 1.0)
    """
    pred_set = set(predicted_pairs)
    ref_set = set(reference_pairs)

    true_positives = len(pred_set & ref_set)
    false_positives = len(pred_set - ref_set)

    if len(pred_set) == 0:
        return 1.0

    return true_positives / (true_positives + false_positives)


def evaluate_structure_prediction(predicted_structure: str,
                                reference_structure: str) -> dict:
    """
    Evaluate a predicted structure against a reference structure.

    Args:
        predicted_structure: Predicted dot-bracket string
        reference_structure: Reference dot-bracket string

    Returns:
        Dictionary with evaluation metrics
    """
    if len(predicted_structure) != len(reference_structure):
        raise ValueError("Predicted and reference structures must have same length")

    pred_pairs = structure_to_base_pairs(predicted_structure)
    ref_pairs = structure_to_base_pairs(reference_structure)

    metrics = {
        'predicted_pairs': len(pred_pairs),
        'reference_pairs': len(ref_pairs),
        'f1_score': calculate_f1_score(pred_pairs, ref_pairs),
        'sensitivity': calculate_sensitivity(pred_pairs, ref_pairs),
        'ppv': calculate_ppv(pred_pairs, ref_pairs),
    }

    return metrics


def print_evaluation_report(predicted_structure: str,
                          reference_structure: str,
                          sequence: str = None) -> None:
    """
    Print a formatted evaluation report.

    Args:
        predicted_structure: Predicted dot-bracket string
        reference_structure: Reference dot-bracket string
        sequence: Optional RNA sequence for context
    """
    metrics = evaluate_structure_prediction(predicted_structure, reference_structure)

    print("Structure Evaluation Report")
    print("=" * 40)
    if sequence:
        print(f"Sequence: {sequence}")
    print(f"Reference:  {reference_structure}")
    print(f"Predicted:  {predicted_structure}")
    print()
    print("Base Pairs:")
    print(f"  Reference: {metrics['reference_pairs']}")
    print(f"  Predicted: {metrics['predicted_pairs']}")
    print()
    print("Metrics:")
    print(".3f")
    print(".3f")
    print(".3f")
    print()