"""
Inspect and validate the mRNA dataset.

Usage:
    python scripts/inspect_dataset.py
"""

from pathlib import Path
import sys

# Add src to path so we can import mrnafold
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from mrnafold.data_loader import load_dataset, validate_dataset


def main():
    """Load dataset and print inspection summary."""
    dataset_path = Path(__file__).parent.parent / "data" / "raw" / "mrna_sequence_dataset.csv"
    
    print("=" * 70)
    print("mRNA Dataset Inspection")
    print("=" * 70)
    print(f"\nDataset path: {dataset_path}")
    print(f"Exists: {dataset_path.exists()}\n")
    
    if not dataset_path.exists():
        print(f"ERROR: Dataset not found at {dataset_path}")
        sys.exit(1)
    
    # Load dataset
    try:
        df, load_warnings = load_dataset(dataset_path)
        print(f"✓ Successfully loaded {len(df)} rows")
    except Exception as e:
        print(f"✗ Failed to load dataset: {e}")
        sys.exit(1)
    
    if load_warnings:
        print(f"⚠ Load warnings: {len(load_warnings)}")
        for w in load_warnings:
            print(f"  - {w}")
    
    print(f"\nColumns: {list(df.columns)}\n")
    
    # Validate dataset
    valid_count, invalid_count, issues = validate_dataset(df)
    
    print(f"Validation Results:")
    print(f"  Valid rows:   {valid_count}")
    print(f"  Invalid rows: {invalid_count}")
    
    if issues:
        print(f"\n⚠ Validation issues: {len(issues)}")
        for issue in issues[:10]:  # Show first 10
            print(f"  - {issue}")
        if len(issues) > 10:
            print(f"  ... and {len(issues) - 10} more")
    else:
        print(f"✓ No validation issues found")
    
    # Print statistics
    print(f"\n" + "=" * 70)
    print("Dataset Statistics")
    print("=" * 70)
    
    lengths = df["length"].describe()
    print(f"\nSequence Length Statistics:")
    print(f"  Count:   {int(lengths['count'])}")
    print(f"  Min:     {int(lengths['min'])}")
    print(f"  Max:     {int(lengths['max'])}")
    print(f"  Mean:    {lengths['mean']:.1f}")
    print(f"  Median:  {df['length'].median():.0f}")
    
    energies = df["energy"].describe()
    print(f"\nReference Energy Statistics (kcal/mol):")
    print(f"  Count:   {int(energies['count'])}")
    print(f"  Min:     {energies['min']:.2f}")
    print(f"  Max:     {energies['max']:.2f}")
    print(f"  Mean:    {energies['mean']:.2f}")
    print(f"  Median:  {df['energy'].median():.2f}")
    
    # Show sample rows
    print(f"\n" + "=" * 70)
    print("Sample Rows (first 3)")
    print("=" * 70)
    for i, (idx, row) in enumerate(df.head(3).iterrows()):
        print(f"\nRow {idx}:")
        print(f"  Sequence: {row['sequence'][:50]}{'...' if len(row['sequence']) > 50 else ''}")
        print(f"  Length: {row['length']}")
        print(f"  Structure: {row['structure'][:50]}{'...' if len(row['structure']) > 50 else ''}")
        print(f"  Energy: {row['energy']} kcal/mol")
    
    print("\n" + "=" * 70)
    print("✓ Inspection complete")
    print("=" * 70)


if __name__ == "__main__":
    main()
