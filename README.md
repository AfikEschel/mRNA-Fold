# mRNA Secondary Structure Prediction: A Quantum-Classical Replication Study

## Project Purpose

This repository is a **learning and replication project** by team **apoqs**. It reconstructs the IBM–Moderna mRNA folding pipeline in a transparent, modular way.

The code includes:
- classical sequence preprocessing
- quartet generation and QUBO formulation
- a classical brute-force baseline for small instances
- a working CVaR-VQE quantum solver on Qiskit simulator
- evaluation metrics and postprocessing to dot-bracket format

## Current Status (v0.2.0)

- ✅ Phase 1–10 complete
- ✅ Data loading and validation
- ✅ Quartet preprocessing and QUBO construction
- ✅ Classical brute-force solver for small cases
- ✅ Working CVaR-VQE quantum solver with two-local ansatz
- ✅ Benchmark script for quantum vs classical comparison
- ✅ Qiskit AerSimulator integration

## Quick commands

```bash
source .venv/bin/activate
pytest tests/ -v
PYTHONPATH=src python scripts/inspect_dataset.py
PYTHONPATH=src python scripts/demo_qubo_pipeline.py
PYTHONPATH=src python scripts/benchmark_quantum_vs_classical.py
```

## Setup

### Prerequisites
- Python 3.8+
- pip

### Install

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -e .
```

If you prefer not to install the package:

```bash
pip install -r requirements.txt
export PYTHONPATH="$(pwd)/src:$PYTHONPATH"
```

## Validation commands

Inspect dataset:

```bash
PYTHONPATH=src python scripts/inspect_dataset.py
```

Run all tests:

```bash
pytest tests/ -v
```

Run the quantum benchmark:

```bash
PYTHONPATH=src python scripts/benchmark_quantum_vs_classical.py
```

## Project structure

```
mrnafold_rep/
├── README.md
├── requirements.txt
├── pyproject.toml
├── data/
│   ├── raw/
│   │   └── mrna_sequence_dataset.csv
│   └── processed/
├── references/
│   └── ibm_moderna_paper.md
├── scripts/
│   ├── inspect_dataset.py
│   ├── demo_qubo_pipeline.py
│   ├── run_pipeline_on_dataset.py
│   └── benchmark_quantum_vs_classical.py
├── src/mrnafold/
│   ├── __init__.py
│   ├── data_loader.py
│   ├── pairing.py
│   ├── quartets.py
│   ├── qubo.py
│   ├── metrics.py
│   └── quantum/
│       ├── __init__.py
│       ├── base_solver.py
│       ├── cvar_vqe.py
│       └── ansatz.py
└── tests/
    ├── test_pairing.py
    ├── test_quartets.py
    ├── test_qubo.py
    ├── test_metrics.py
    └── test_quantum.py
```

## Limitations and honesty

- Simulator-only quantum execution: no real hardware runs yet
- Simplified energy scoring: ViennaRNA is not integrated
- No pseudoknots: nested structures only
- Proof-of-concept: not a production tool
- No quantum advantage claims

## References

- IBM–Moderna paper: "mRNA secondary structure prediction using utility-scale quantum computers"
- IBM–Moderna Hackathon Challenge dataset

## License

TBD

## Team

**apoqs** — replication and learning for mRNA folding with quantum and classical methods.
