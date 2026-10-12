# Quantum Computing CA2 — QFT, Period Finding & Toffoli Gate

## Project Overview

A university quantum computing project developed using Python and Qiskit to explore periodic quantum states, the Quantum Fourier Transform (QFT), and reversible quantum logic.

The project demonstrates key concepts used in quantum algorithms, including Shor-style period finding and three-qubit Toffoli gate operations.

## Key Features

- **Periodic State Preparation:** Construct quantum superpositions with regularly spaced basis states.
- **Quantum Fourier Transform:** Apply the QFT to reveal periodic structure through measurement peaks.
- **Periodicity Analysis:** Compare cases where the period divides the register size and where it does not.
- **Toffoli Gate Verification:** Test the controlled-controlled-NOT (CCX) gate against its classical truth table.
- **Superposition Testing:** Examine Toffoli gate behavior on a quantum superposition.
- **Result Visualization:** Generate quantum circuit diagrams and measurement histograms.

## Technologies Used

- Python
- Qiskit
- Qiskit Aer
- NumPy
- Matplotlib

## Simulation Results

The simulations demonstrated:

- Correct periodic state preparation for both tested periods.
- QFT measurement peaks consistent with the expected periodic structure.
- Correct Toffoli gate behavior across all eight computational basis states.
- Approximately equal measurement probabilities for the two expected superposition outcomes.

## Project Files

- `CA2-1.py` — Main Python implementation.
- `task1_*.png` — Periodic state circuits and histograms.
- `task2_*.png` — QFT circuits and measurement histograms.
- `toffoli_*.png` — Toffoli gate circuit and superposition results.

## Running the Project

Install the required Python packages:

```bash
pip install qiskit qiskit-aer numpy matplotlib
```

Run the main script:

```bash
python CA2-1.py
```

The program generates simulation results and saves circuit diagrams and histograms to the `outputs` directory.

## Project Context

Developed as part of university coursework in Quantum Computing and Security at the University of Liverpool.
