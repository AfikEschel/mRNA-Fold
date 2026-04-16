#!/usr/bin/env python3
"""Demo script running the complete pipeline on dataset sequences."""

import sys
sys.path.insert(0, 'src')

from mrnafold.data_loader import load_dataset
from mrnafold.quartets import generate_quartets
from mrnafold.qubo import build_quartet_qubo, solve_qubo_brute_force, quartets_to_structure, structure_to_dot_bracket
from mrnafold.metrics import evaluate_structure_prediction, print_evaluation_report


def run_pipeline_on_sequence(sequence: str, reference_structure: str, sequence_name: str = "Test"):
    """Run the complete pipeline on a single sequence."""
    print(f"\n{'='*60}")
    print(f"Running pipeline on {sequence_name}")
    print(f"Sequence: {sequence}")
    print(f"Reference: {reference_structure}")
    print(f"Length: {len(sequence)}")
    print('='*60)

    # Phase 1-4: Generate quartets
    print("\n1. Quartet Generation:")
    quartets = generate_quartets(sequence)
    print(f"   Generated {len(quartets)} valid quartets")

    if len(quartets) == 0:
        print("   No quartets found - cannot proceed with QUBO")
        return

    # Phase 5: Build QUBO
    print("\n2. QUBO Formulation:")
    qubo = build_quartet_qubo(quartets, conflict_penalty=1.0, stacking_reward=0.5)
    print(f"   QUBO matrix: {qubo.shape[0]}x{qubo.shape[1]}")

    # Solve QUBO
    print("\n3. QUBO Optimization:")
    if len(quartets) > 20:
        print(f"   Skipping brute force optimization (n={len(quartets)} > 20)")
        print("   QUBO formulation complete - ready for quantum/approximate solvers")
        return
    
    solution, energy = solve_qubo_brute_force(qubo)
    print(f"   Optimal solution found with energy: {energy}")

    # Convert to structure
    print("\n4. Structure Prediction:")
    selected_quartets = [q for q, selected in zip(quartets, solution) if selected]
    print(f"   Selected {len(selected_quartets)} quartets")

    base_pairs = quartets_to_structure(selected_quartets, len(sequence))
    predicted_structure = structure_to_dot_bracket(base_pairs, len(sequence))
    print(f"   Predicted structure: {predicted_structure}")

    # Evaluate
    print("\n5. Evaluation:")
    metrics = evaluate_structure_prediction(predicted_structure, reference_structure)
    print_evaluation_report(predicted_structure, reference_structure, sequence)


def main():
    """Run pipeline on dataset sequences."""
    print("mRNA Folding Pipeline Demo")
    print("Loading dataset...")

    # Load dataset
    df, warnings = load_dataset('data/raw/mrna_sequence_dataset.csv')
    if warnings:
        print(f"Warnings: {warnings}")
    print(f"Loaded {len(df)} sequences")

    # Run on shortest sequence
    shortest_row = df.loc[df['length'].idxmin()]
    sequence = shortest_row['sequence']
    reference = shortest_row['structure']

    run_pipeline_on_sequence(sequence, reference, f"Shortest sequence (length {len(sequence)})")

    # Try to find a sequence with manageable quartets for full pipeline
    print("\n" + "="*60)
    print("Looking for sequence with <=20 quartets for full pipeline demo...")
    
    from mrnafold.quartets import generate_quartets
    for idx, row in df.iterrows():
        quartets = generate_quartets(row['sequence'])
        if len(quartets) <= 20 and len(quartets) > 0:
            print(f"Found suitable sequence: length {row['length']}, quartets: {len(quartets)}")
            run_pipeline_on_sequence(row['sequence'], row['structure'], 
                                   f"Demo sequence (length {row['length']}, {len(quartets)} quartets)")
            break
    else:
        print("No sequence found with <=20 quartets. All have too many for brute force demo.")


if __name__ == '__main__':
    main()