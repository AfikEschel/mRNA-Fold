#!/usr/bin/env python3
"""Demo script showing the complete QUBO pipeline for mRNA folding."""

import sys
sys.path.insert(0, 'src')

from mrnafold.quartets import generate_quartets
from mrnafold.qubo import build_quartet_qubo, solve_qubo_brute_force, quartets_to_structure, structure_to_dot_bracket


def main():
    """Demonstrate the complete QUBO pipeline."""
    # Use a small sequence for demonstration
    seq = 'AUGC'  # Length 4
    print(f"Sequence: {seq} (length {len(seq)})")
    print()

    # Phase 1-4: Generate quartets
    print("Phase 1-4: Quartet Generation")
    quartets = generate_quartets(seq)
    print(f"Generated {len(quartets)} valid quartets:")
    for i, q in enumerate(quartets):
        print(f"  Q{i}: {q}")
    print()

    # Phase 5: Build QUBO
    print("Phase 5: QUBO Formulation")
    qubo = build_quartet_qubo(quartets, conflict_penalty=1.0, stacking_reward=0.5)
    print(f"QUBO matrix ({qubo.shape[0]}x{qubo.shape[1]}):")
    print(qubo)
    print()

    if qubo.shape[0] == 0:
        print("No quartets were generated, so there is no QUBO to solve.")
        print("Selected quartets: []")
        print(f"Base pairs: []")
        print(f"Dot-bracket: {'.' * len(seq)}")
        print()
        print("Pipeline complete! 🎉")
        return

    # Solve QUBO (brute force for small instance)
    print("Solving QUBO (brute force):")
    solution, energy = solve_qubo_brute_force(qubo)
    print(f"Optimal solution: {solution}")
    print(f"Minimum energy: {energy}")
    print()

    # Convert to structure
    print("Converting to RNA structure:")
    selected_quartets = [q for q, selected in zip(quartets, solution) if selected]
    print(f"Selected quartets: {[str(q) for q in selected_quartets]}")

    base_pairs = quartets_to_structure(selected_quartets, len(seq))
    print(f"Base pairs: {base_pairs}")

    dot_bracket = structure_to_dot_bracket(base_pairs, len(seq))
    print(f"Dot-bracket: {dot_bracket}")
    print()

    print("Pipeline complete! 🎉")


if __name__ == '__main__':
    main()