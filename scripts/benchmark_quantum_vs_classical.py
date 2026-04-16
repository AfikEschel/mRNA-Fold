#!/usr/bin/env python3
# scripts/benchmark_quantum_vs_classical.py
"""
Benchmark CVaR-VQE quantum solver against classical brute-force baseline.

This script:
1. Loads the shortest sequences from the mRNA dataset
2. Generates quartets and QUBO formulation
3. Solves with both classical brute-force and quantum CVaR-VQE
4. Compares energies and solution quality
5. Produces a report

This is the correctness check for Phase 10: quantum and classical should
agree on tiny instances when run on the simulator (no hardware noise).
"""

import sys
import time
import numpy as np
from typing import List, Tuple

try:
    from mrnafold.data_loader import load_dataset
    from mrnafold.quartets import generate_quartets
    from mrnafold.qubo import build_quartet_qubo, solve_qubo_brute_force, quartets_to_structure, structure_to_dot_bracket
    from mrnafold.quantum import CVaRVQE
except ImportError as e:
    print(f"Error importing mrnafold modules: {e}")
    sys.exit(1)


def benchmark_sequence(
    sequence: str,
    reference_structure: str,
    max_quartets: int = 20,
) -> dict:
    """
    Benchmark a single sequence.

    Args:
        sequence: RNA sequence
        reference_structure: Reference dot-bracket structure
        max_quartets: Max quartets for brute-force (skip if exceeded)

    Returns:
        Dictionary with results
    """
    result = {
        "sequence": sequence,
        "length": len(sequence),
        "quartets": 0,
        "classical_energy": None,
        "classical_time": 0.0,
        "quantum_energy": None,
        "quantum_time": 0.0,
        "energy_diff": None,
        "classical_solved": False,
        "quantum_solved": False,
        "reference_structure": reference_structure,
    }

    # Generate quartets
    try:
        quartets = generate_quartets(sequence)
    except Exception as e:
        print(f"  ❌ Error generating quartets: {e}")
        return result

    if len(quartets) == 0:
        print(f"  ℹ️  No quartets found for this sequence")
        return result

    result["quartets"] = len(quartets)
    print(f"  Quartets generated: {len(quartets)}")

    # Build QUBO
    try:
        qubo = build_quartet_qubo(quartets, conflict_penalty=1.0, stacking_reward=0.5)
    except Exception as e:
        print(f"  ❌ Error building QUBO: {e}")
        return result

    # Classical solver (if small enough)
    if len(quartets) <= max_quartets:
        print(f"  Running classical brute-force solver...")
        try:
            start = time.time()
            classical_bitstring, classical_energy = solve_qubo_brute_force(qubo)
            result["classical_time"] = time.time() - start
            result["classical_energy"] = classical_energy
            result["classical_solved"] = True
            print(f"    ✓ Classical energy: {classical_energy:.6f} (time: {result['classical_time']:.3f}s)")
        except Exception as e:
            print(f"    ❌ Classical solver error: {e}")
    else:
        print(f"  ⊘ Skipping classical (too many quartets: {len(quartets)} > {max_quartets})")

    # Quantum solver
    print(f"  Running CVaR-VQE quantum solver...")
    try:
        start = time.time()
        solver = CVaRVQE(
            num_qubits=len(quartets),
            num_layers=2,
            alpha=0.1,
            shots=1024,
            max_iterations=100,
            seed=42,
        )
        quantum_bitstring, quantum_energy = solver.solve(qubo)
        result["quantum_time"] = time.time() - start
        result["quantum_energy"] = quantum_energy
        result["quantum_solved"] = True
        print(f"    ✓ Quantum energy: {quantum_energy:.6f} (time: {result['quantum_time']:.3f}s)")

        # Compute difference
        if result["classical_solved"]:
            result["energy_diff"] = abs(quantum_energy - result["classical_energy"])
            print(f"    Energy diff (quantum - classical): {result['energy_diff']:.6f}")

    except Exception as e:
        print(f"    ❌ Quantum solver error: {e}")

    return result


def main():
    """Main benchmark routine."""
    print("=" * 70)
    print("CVaR-VQE Quantum vs Classical Brute-Force Benchmark")
    print("=" * 70)
    print()

    # Load dataset
    print("Loading mRNA dataset...")
    try:
        df, warnings = load_dataset("data/raw/mrna_sequence_dataset.csv")
        print(f"✓ Loaded {len(df)} sequences")
        if warnings:
            print(f"  Warnings: {warnings}")
    except Exception as e:
        print(f"❌ Error loading dataset: {e}")
        return

    # Convert DataFrame to list of dicts
    data = df.to_dict('records')

    # Sort by sequence length
    data_sorted = sorted(data, key=lambda x: x["length"])

    print()
    print(f"Benchmarking shortest sequences (target: <= 20 quartets)...")
    print()

    results = []
    for seq_idx, row in enumerate(data_sorted):
        sequence = row["sequence"]
        reference = row["structure"]
        print(f"[{seq_idx + 1}/{len(data)}] Sequence length {len(sequence)}: {sequence[:30]}...")

        result = benchmark_sequence(sequence, reference)
        results.append(result)

        print()

        # Stop after a few successful benchmarks
        if len([r for r in results if r["quantum_solved"]]) >= 2:
            print("✓ Sufficient benchmarks completed")
            break

    # Summary report
    print()
    print("=" * 70)
    print("BENCHMARK SUMMARY")
    print("=" * 70)
    print()

    completed = [r for r in results if r["quantum_solved"]]
    if not completed:
        print("❌ No successful quantum solves")
        return

    print(f"Completed: {len(completed)} sequences")
    print()

    for r in completed:
        print(f"  Sequence: {r['sequence'][:40]}... (len={r['length']})")
        print(f"  Quartets: {r['quartets']}")
        if r["classical_solved"]:
            print(f"    Classical: {r['classical_energy']:8.4f} (time: {r['classical_time']:6.3f}s)")
            print(f"    Quantum:   {r['quantum_energy']:8.4f} (time: {r['quantum_time']:6.3f}s)")
            print(f"    Diff:      {r['energy_diff']:8.4f}")
        else:
            print(f"    Quantum:   {r['quantum_energy']:8.4f} (time: {r['quantum_time']:6.3f}s)")
        print()

    # Statistics
    classical_energies = [r["classical_energy"] for r in completed if r["classical_solved"]]
    quantum_energies = [r["quantum_energy"] for r in completed if r["quantum_solved"]]
    diffs = [r["energy_diff"] for r in completed if r["energy_diff"] is not None]

    print("Statistics:")
    print(f"  Quantum solves:     {len(completed)}")
    print(f"  Classical solves:   {len(classical_energies)}")
    if diffs:
        print(f"  Avg energy diff:    {np.mean(diffs):.6f}")
        print(f"  Max energy diff:    {np.max(diffs):.6f}")
        print(f"  Min energy diff:    {np.min(diffs):.6f}")

    print()
    print(f"✓ Benchmark complete")
    print("=" * 70)


if __name__ == "__main__":
    main()