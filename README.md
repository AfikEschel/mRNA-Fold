# mRNA Secondary Structure Prediction: A Quantum-Classical Replication Study

## Project Purpose

This repository is a **learning and replication project** by team **apoqs** to incrementally reconstruct the computational pipeline from the IBM–Moderna paper on mRNA secondary structure prediction using quantum computers.

**This is NOT a full reproduction yet.** It is an honest, modular, and transparent learning effort that:
- Starts from classical preprocessing and problem formulation
- Implements small classical baselines for validation
- **Includes a working quantum layer** (CVaR-VQE on simulator)
- Prioritizes reproducibility and honesty over ambition

## Scientific Context

### Problem Statement

Predicting the minimum-free-energy (MFE) secondary structure of mRNA from its nucleotide sequence is an NP-complete problem critical for designing RNA-based therapeutics. The IBM–Moderna paper demonstrates that variational quantum algorithms (CVaR-VQE with NFT optimizer and two-local ansatz) can solve this problem on utility-scale quantum processors.

### Key Approach

The pipeline follows:
1. **Classical preprocessing**: Extract valid base pairs and quartets from sequence
2. **Problem formulation**: Model as quartet-based binary optimization (QUBO)
3. **Optimization**: Solve with classical baseline (brute force for small instances)
4. **Postprocessing**: Convert selected quartets → base pairs → dot-bracket notation
5. **Evaluation**: Compare against reference structures and compute energy

### Quartet Encoding

A **quartet** is two nested consecutive base pairs: (i, j) and (i+1, j-1). This encoding is central to:
- The paper's formulation
- Preprocessing (valid quartets Q, stackable QS(qi), conflicting QC(qi))
- QUBO construction

### Base Pairing Rules

- **Canonical pairs**: A-U, U-A, C-G, G-C
- **Wobble pairs**: G-U, U-G

## Relationship to IBM–Moderna Paper

This repo is inspired by:
- **Paper**: "mRNA secondary structure prediction using utility-scale quantum computers" (Alevras et al., IBM Quantum & Moderna)
- **Quantum algorithm**: CVaR-VQE with Nakanishi-Fujii-Todo (NFT) optimizer
- **Ansatz**: Hardware-efficient two-local
- **Classical baseline**: CPLEX (we use brute force for proof-of-concept)
- **Hardware**: IBM Eagle and Heron processors (simulated here)

See [references/ibm_moderna_paper.md](references/ibm_moderna_paper.md) for full paper text.

## Relationship to Hackathon Challenge

The project is grounded in the IBM–Moderna Hackathon Challenge, which defines the pipeline as:
```
classical preprocessing → quantum sampler → classical postprocessing
```

The challenge explicitly recommends:
- Start with small proof-of-concept examples
- Benchmark against classical solvers
- Pursue reproducibility and transparency
- Build a compelling path toward possible quantum advantage (not claim it immediately)

## What Is Implemented (v0.2.0)

<<<<<<< HEAD
### Phase 1–10 (Current)
=======
### Phase 1–10 (Complete)
>>>>>>> quantum-layer
- ✅ Repository structure and file layout
- ✅ Robust data loading from semicolon-separated CSV (with fuzzy column-name matching)
- ✅ Dataset validation (sequence length vs. Length column, NA checks)
- ✅ Dataset inspection script with summary statistics
- ✅ Base pairing rules (canonical + wobble pairs)
- ✅ Quartet dataclass and generation from sequences
- ✅ Quartet preprocessing (Q, QS, QC) with conflict and stacking detection
- ✅ Comprehensive unit tests (pairing: 17 tests, quartets: 15 tests, qubo: 8 tests, quantum: 20 tests)
- ✅ QUBO formulation with conflict penalties and stacking rewards
- ✅ QUBO-to-Ising conversion
- ✅ Brute-force solver for small instances
- ✅ Structure conversion (quartets → base pairs → dot-bracket)
- ✅ Metrics and evaluation (F1, sensitivity, PPV for structure comparison)
<<<<<<< HEAD
- ✅ Scaffold quantum module (CVaR-VQE and two-local ansatz placeholders, marked as future work)
=======
- ✅ **CVaR-VQE quantum solver** with two-local ansatz (Y rotations + CZ gates, p=2 layers)
- ✅ Quantum-classical benchmark script comparing CVaR-VQE vs brute-force baseline
- ✅ Qiskit AerSimulator integration (statevector for ≤20 qubits, matrix_product_state for larger)

### Future Phases (v1.0)
- ViennaRNA integration for thermodynamic energy calculation
- Multiple quantum solvers (QAOA, VQE variants)
- Fujitsu simulator support
- Hardware execution capabilities
>>>>>>> quantum-layer

## What Is NOT Implemented

- ❌ Full thermodynamic energy scoring (will use simplified models in v0.2.0)
- ❌ Pseudoknot handling (future work)
- ❌ Claims of quantum advantage
- ❌ Production-level error handling or performance optimization
- ❌ Hardware execution (simulator-only in v0.2.0)
- ❌ Multiple quantum solvers (CVaR-VQE implemented, others planned for v1.0)

## Quick Start

### Quick commands
```bash
source .venv/bin/activate
pytest tests/ -v
PYTHONPATH=src python scripts/inspect_dataset.py
PYTHONPATH=src python scripts/demo_qubo_pipeline.py
PYTHONPATH=src python scripts/benchmark_quantum_vs_classical.py
```

### Prerequisites
- Python 3.8+
- `pip`

### Setup

1. Clone the repository:
```bash
git clone <repo-url>
cd mrnafold
```

2. Create a virtual environment:
```bash
python3 -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate
```

3. Install the package in development mode:
```bash
pip install -e .
```

This command reads `pyproject.toml` and installs the `mrnafold` package along with all dependencies. 
The `-e` flag installs in "editable" mode, so changes to source files are immediately reflected.

If you close the terminal or open a new session, reactivate the virtual environment before running commands:
```bash
source venv/bin/activate
```

**Alternative** (if you prefer not to install the package):
```bash
pip install -r requirements.txt
export PYTHONPATH="$(pwd)/src:$PYTHONPATH"  # Add src to Python path
```

### First Run: Inspect Dataset

Verify that the data loading pipeline works:
```bash
python scripts/inspect_dataset.py
```

Expected output:
- Number of sequences
- Sequence length statistics
- Sample rows
- Validation summary (length mismatches, missing values)

### Run Tests

Make sure the virtual environment is active:
```bash
source venv/bin/activate
pytest tests/
```

The tests should all pass (78 tests currently: 17 pairing + 15 quartets + 8 QUBO + 4 metrics + 20 quantum + 14 total checks).

### Run Quantum Benchmark

Test the quantum solver against classical baseline:
```bash
source venv/bin/activate
PYTHONPATH=src python scripts/benchmark_quantum_vs_classical.py
```

This will run CVaR-VQE on real mRNA sequences from the dataset and compare against brute-force solving.

### Directory Structure

```
mrnafold/
├── README.md                          # This file
├── requirements.txt                   # Python dependencies
├── .gitignore                         # Git ignore rules
│
├── data/
│   ├── raw/
│   │   └── mrna_sequence_dataset.csv  # Input dataset (semicolon-separated)
│   └── processed/                     # Output datasets (generated)
│
├── references/
│   └── ibm_moderna_paper.md           # IBM–Moderna paper text
│
├── src/mrnafold/
│   ├── __init__.py                    # Package init
│   ├── data_loader.py                 # Dataset loading and validation
│   ├── pairing.py                     # Base pair validation (canonical + wobble)
│   ├── quartets.py                    # Quartet representation and generation
│   ├── qubo.py                        # QUBO formulation and classical baseline
│   ├── metrics.py                     # Structure comparison and evaluation
│   └── quantum/                       # Quantum solvers module
│       ├── __init__.py                # Quantum module exports
│       ├── base_solver.py             # Abstract solver interface
│       ├── cvar_vqe.py                # CVaR-VQE implementation
│       └── ansatz.py                  # Two-local ansatz
│
├── scripts/
│   ├── inspect_dataset.py             # Inspect and validate dataset
│   ├── demo_qubo_pipeline.py          # Demo QUBO pipeline for a sample sequence
│   ├── run_pipeline_on_dataset.py     # Run pipeline on dataset examples
│   └── benchmark_quantum_vs_classical.py  # Quantum vs classical benchmark
│
├── notebooks/
│   └── ... (exploratory notebooks, optional)
│
└── tests/
    ├── test_pairing.py                # Unit tests for pairing rules
    ├── test_quartets.py               # Unit tests for quartet logic
    ├── test_qubo.py                    # Unit tests for QUBO and structure conversion
    ├── test_metrics.py                # Unit tests for evaluation metrics
    └── test_quantum.py                # Unit tests for quantum solvers
```

## Limitations and Honesty Statements

### Current Implementation (v0.2.0)
- **Simulator-only quantum execution**: Quantum solver runs on Qiskit AerSimulator, not real hardware.
- **Simplified energy scoring**: We do not yet integrate full ViennaRNA thermodynamic tables. Energy values are placeholders or simplified approximations.
- **Medium instances only**: Classical baseline is brute-force; suitable for ≤ 20 quartets. Larger instances use quantum solver only.
- **No pseudoknots**: Current formulation excludes pseudoknots (nested only).
- **Proof-of-concept only**: This is NOT a competitor to MFold, ViennaRNA, or other established tools.
- **No quantum advantage claims**: Implementation demonstrates feasibility, not superiority.

### Next Phases
- Proper QUBO construction following the paper
- Integration with ViennaRNA (if available) for energy validation
- Scalable classical solver (e.g., CPLEX, branch-and-bound)
- Error mitigation strategies for quantum layer
- Comprehensive test coverage

## Next Steps (v1 and Beyond)

<<<<<<< HEAD
Now that v0 is complete, future development could include:
1. **ViennaRNA Integration**: Add optional thermodynamic energy calculation using ViennaRNA for accurate MFE validation.
2. **Scalable Classical Solver**: Replace brute-force with CPLEX or other optimizers for larger instances.
3. **Full Quantum Implementation**: Implement CVaR-VQE with NFT optimizer and two-local ansatz using Qiskit or similar.
4. **Pseudoknot Support**: Extend formulation to handle pseudoknots.
5. **Hardware Simulation**: Test on quantum simulators and eventually real hardware.
6. **Benchmarking**: Compare against MFold, ViennaRNA on larger datasets.
7. **Error Mitigation**: Add strategies for noisy quantum devices.
=======
1. **Review Phase 10**: Verify quantum solver works correctly on dataset sequences.
2. **V1.0 Development**:
   - Integrate ViennaRNA for accurate energy calculations
   - Add multiple quantum solvers (QAOA, VQE variants)
   - Implement Fujitsu simulator support
   - Add hardware execution capabilities
3. **Test thoroughly**: All tests pass, benchmark script runs successfully.
4. **Document assumptions**: Make simplifications explicit in code comments.
5. **Future work**:
   - Pseudoknot handling
   - Production-level optimizations
   - Error mitigation strategies
>>>>>>> quantum-layer

## How to Contribute / Extend

- Keep files small and modular
- Use type hints and dataclasses
- Add docstrings and reference the paper where relevant
- Write tests for new logic
- Mark TODOs, placeholders, and simplifications clearly
- Do not over-engineer; prefer transparent code over clever code

## References

- **Paper**: Alevras, D., Metkar, M., Yamamoto, T., et al. "mRNA secondary structure prediction using utility-scale quantum computers." *IBM Quantum & Moderna* (2024).
- **Quantum algorithms**: VQE, CVaR-VQE, QAOA, VQA
- **Classical tools**: ViennaRNA, MFold, CPLEX
- **Dataset**: IBM–Moderna Hackathon Challenge dataset

## License

TBD (check with team before publication)

## Team

**apoqs** — learning mRNA folding with quantum and classical methods.

---

**Last Updated**: April 2026  
<<<<<<< HEAD
**Status**: v0 complete (Phase 1–10 done)
=======
**Status**: v0 (Phase 1–10 complete)
>>>>>>> quantum-layer
