"""
Data loading and validation for mRNA sequences.

This module handles reading the semicolon-separated dataset,
normalizing column names (with fuzzy matching for typos),
and validating data consistency.
"""

from dataclasses import dataclass
from pathlib import Path
from typing import Optional
import pandas as pd


@dataclass
class MRNASequence:
    """
    Dataclass representing a single mRNA sequence record.
    
    Attributes:
        sequence: nucleotide sequence (AUGC)
        length: length column from dataset
        reference_structure: dot-bracket notation from dataset
        reference_energy: MFE in kcal/mol from dataset
    """
    sequence: str
    length: int
    reference_structure: str
    reference_energy: float


def load_dataset(path: Path | str) -> tuple[pd.DataFrame, list[str]]:
    """
    Load mRNA dataset from semicolon-separated CSV.
    
    Returns:
        (DataFrame, list of validation warnings)
    
    Raises:
        FileNotFoundError: if path does not exist
        ValueError: if critical columns are missing
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Dataset not found: {path}")
    
    # Read with semicolon separator
    df = pd.read_csv(path, sep=";")
    
    warnings = []
    
    # Normalize column names: handle typo "Dot-braket" (missing 'c')
    # Strategy: fuzzy match column names
    col_map = {}
    expected_cols = {
        "sequence": ["Sequence"],
        "length": ["Length"],
        "structure": ["Dot-bracket min free energy conformation", 
                     "Dot-braket min free energy conformation"],  # typo version
        "energy": ["Energy (kcal/mol)", "Energy(kcal/mol)"]
    }
    
    # Match each expected column
    actual_cols = {col.lower() for col in df.columns}
    
    for key, candidates in expected_cols.items():
        found = False
        for candidate in candidates:
            if candidate.lower() in actual_cols:
                col_map[key] = candidate
                found = True
                break
        if not found:
            # Try fuzzy matching: look for partial matches
            for col in df.columns:
                col_lower = col.lower()
                if any(x.lower() in col_lower for x in key.split()):
                    col_map[key] = col
                    found = True
                    break
        if not found:
            raise ValueError(
                f"Column '{key}' not found. Available: {list(df.columns)}"
            )
    
    # Rename to canonical names
    df = df.rename(columns=col_map)
    df = df.rename(columns={v: k for k, v in col_map.items()})
    
    # Keep only required columns
    df = df[["sequence", "length", "structure", "energy"]]
    
    return df, warnings


def validate_dataset(df: pd.DataFrame) -> tuple[int, int, list[str]]:
    """
    Validate dataset consistency.
    
    Returns:
        (num_valid, num_invalid, list of validation issues)
    """
    issues = []
    invalid_count = 0
    
    for idx, row in df.iterrows():
        seq = row["sequence"]
        declared_len = row["length"]
        actual_len = len(seq)
        
        # Check length mismatch
        if actual_len != declared_len:
            issues.append(
                f"Row {idx}: sequence length {actual_len} != declared {declared_len}"
            )
            invalid_count += 1
        
        # Check for invalid nucleotides
        valid_nucleotides = set("AUGC")
        if not all(n in valid_nucleotides for n in seq):
            issues.append(f"Row {idx}: invalid nucleotides in sequence")
            invalid_count += 1
        
        # Check for missing structure or energy
        if pd.isna(row["structure"]) or pd.isna(row["energy"]):
            issues.append(f"Row {idx}: missing structure or energy")
            invalid_count += 1
    
    valid_count = len(df) - invalid_count
    return valid_count, invalid_count, issues


def load_sequences(path: Path | str) -> list[MRNASequence]:
    """
    Load dataset and return list of MRNASequence objects.
    """
    df, warnings = load_dataset(path)
    sequences = []
    for _, row in df.iterrows():
        seq = MRNASequence(
            sequence=row["sequence"],
            length=row["length"],
            reference_structure=row["structure"],
            reference_energy=row["energy"]
        )
        sequences.append(seq)
    return sequences
