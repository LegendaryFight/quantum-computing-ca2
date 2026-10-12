# CA2 – Periodic Superposition, QFT Peaks, and Toffoli Verification
# Tasks:
#   1. Prepare periodic quantum states (Shor-style collapsed register)
#   2. Apply QFT and observe peaks that encode the period
#   3. Verify Toffoli (CCX) on basis states and a superposition
# All simulations use Qiskit's AerSimulator.

import os
from collections import Counter
import numpy as np

from qiskit import QuantumCircuit, transpile
from qiskit_aer import AerSimulator
from qiskit.circuit.library import StatePreparation, QFT
from qiskit.visualization import plot_histogram

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SEED=42
SHOTS=4096
OUTDIR="outputs"
os.makedirs(OUTDIR,exist_ok=True)

sim=AerSimulator(seed_simulator=SEED,method="automatic")

def run_counts(qc,shots=SHOTS,label=None,filename=None):
    tqc=transpile(qc,sim,optimization_level=1,seed_transpiler=SEED)
    result=sim.run(tqc,shots=shots).result()
    counts=result.get_counts()
    if filename:
        fig=plot_histogram(counts,title=label or "Histogram")
        fig.savefig(os.path.join(OUTDIR,filename),bbox_inches="tight",dpi=160)
        plt.close(fig)
    return counts

def save_circuit_draw(qc,name):
    fig=qc.draw("mpl")
    fig.savefig(os.path.join(OUTDIR,name),bbox_inches="tight",dpi=160)
    plt.close(fig)

def top_peaks(counts,k=6):
    total=sum(counts.values())
    items=[(b,v/total) for b,v in counts.items()]
    return sorted(items,key=lambda x:x[1],reverse=True)[:k]

# Task 1 — Periodic State Preparation
# Build |ψ> = (1/√K) Σ_k |x0 + k*r>, with x0+kr < 2^n
def periodic_state_amplitudes(n,r,x0):
    Q=2**n
    inds=[]
    x=x0
    while x<Q:
        inds.append(x)
        x+=r
    K=len(inds)
    amps=np.zeros(Q,dtype=complex)
    if K==0:
        raise ValueError("No indices fit within 2**n; choose smaller r or x0.")
    amps[inds]=1/np.sqrt(K)
    return amps,inds

def build_periodic_circuit(n,r,x0,measure=True):
    # Always allocate n classical bits so later measure() calls are valid
    amps,inds=periodic_state_amplitudes(n,r,x0)
    sp=StatePreparation(amps)
    qc=QuantumCircuit(n,n)
    qc.append(sp,range(n))
    if measure:
        qc.measure(range(n),range(n))
    return qc,inds

def verify_periodicity(counts,r,x0,n):
    Q=2**n
    support={int(b,2) for b in counts.keys()}
    expected=set()
    x=x0
    while x<Q:
        expected.add(x)
        x+=r
    only_expected=all((s in expected) for s in support)
    all_present=all(((format(e,f"0{n}b") in counts) for e in expected))
    ok=only_expected and all_present
    return ok,expected,support

# Task 2 — Apply QFT and analyze peaks
def apply_qft_to_state(n,periodic_qc_no_meas):
    qc=periodic_qc_no_meas.copy()
    qc.append(QFT(num_qubits=n,do_swaps=True),range(n))
    qc.measure(range(n),range(n))
    return qc

def rational_peaks(Q,r,kmax=6):
    return [int(round(m*Q/r)) for m in range(r)][:kmax]

# Task 3 — Toffoli (CCX) verification
# NOTE: Qiskit prints bitstrings as c_{n-1}...c_0, so we reverse to interpret as (a,b,c).
def toffoli_truth_table():
    truth={}
    for a in (0,1):
        for b in (0,1):
            for c in (0,1):
                q=QuantumCircuit(3,3)
                if a:q.x(0)
                if b:q.x(1)
                if c:q.x(2)
                q.ccx(0,1,2)
                q.measure(range(3),range(3))
                counts=run_counts(q,shots=1)
                raw=next(iter(counts.keys()))
                out_bits=raw[::-1]   # reverse so order is (a,b,c)
                truth[f"{a}{b}{c}"]=out_bits
    return truth

def toffoli_on_superposition():
    q=QuantumCircuit(3,3)
    q.h(0)   # |+> on first qubit
    q.x(1)   # force second qubit to |1>
    q.ccx(0,1,2)
    q.measure(range(3),range(3))
    raw_counts=run_counts(q,shots=SHOTS,label="CCX on |+>|1>|0>",
                          filename="toffoli_superposition_hist.png")
    save_circuit_draw(q,"toffoli_superposition_circuit.png")
    # Reverse keys so printed bitstrings are in (a,b,c) order
    counts={"".join(reversed(k)):v for k,v in raw_counts.items()}
    return counts

if __name__=="__main__":
    print("Task 1 — Periodic State Preparation")

    # Case A: r divides Q (clean QFT peaks)
    n_A,r_A,x0_A=8,16,3
    qcA_no_meas,indices_A=build_periodic_circuit(n_A,r_A,x0_A,measure=False)
    qcA=qcA_no_meas.copy()
    qcA.measure(range(n_A),range(n_A))
    save_circuit_draw(qcA_no_meas,"task1_periodic_state_circuit_divides.png")
    countsA=run_counts(qcA,label=f"Periodic state (n={n_A}, r={r_A}, x0={x0_A})",
                       filename="task1_periodic_state_hist_divides.png")
    okA,expectedA,supportA=verify_periodicity(countsA,r_A,x0_A,n_A)
    print("[Divides] Periodicity OK:",okA)
    print("Expected indices:",sorted(expectedA))
    print("Measured support:",sorted(supportA))
    print("Top outcomes:",top_peaks(countsA,8))

    # Case B: r does NOT divide Q (QFT gives smeared peaks)
    n_B,r_B,x0_B=8,12,1
    qcB_no_meas,indices_B=build_periodic_circuit(n_B,r_B,x0_B,measure=False)
    qcB=qcB_no_meas.copy()
    qcB.measure(range(n_B),range(n_B))
    save_circuit_draw(qcB_no_meas,"task1_periodic_state_circuit_notdivides.png")
    countsB=run_counts(qcB,label=f"Periodic state (n={n_B}, r={r_B}, x0={x0_B})",
                       filename="task1_periodic_state_hist_notdivides.png")
    okB,expectedB,supportB=verify_periodicity(countsB,r_B,x0_B,n_B)
    print("\n[Not divides] Periodicity OK:",okB)
    print("Expected indices:",sorted(expectedB))
    print("Measured support:",sorted(supportB))
    print("Top outcomes:",top_peaks(countsB,8))

    print("\nTask 2 — Apply QFT")
    qftA=apply_qft_to_state(n_A,qcA_no_meas)
    save_circuit_draw(qftA.remove_final_measurements(inplace=False),"task2_qft_circuit_divides.png")
    qft_countsA=run_counts(qftA,label="QFT(|ψ_periodic|) — r divides Q",
                           filename="task2_qft_hist_divides.png")
    print("QFT peaks (divides):",top_peaks(qft_countsA,10))
    print("Expected Q/r multiples:",rational_peaks(2**n_A,r_A,r_A))

    qftB=apply_qft_to_state(n_B,qcB_no_meas)
    save_circuit_draw(qftB.remove_final_measurements(inplace=False),"task2_qft_circuit_notdivides.png")
    qft_countsB=run_counts(qftB,label="QFT(|ψ_periodic|) — r not dividing Q",
                           filename="task2_qft_hist_notdivides.png")
    print("QFT peaks (not divides):",top_peaks(qft_countsB,10))
    print("Expected approximate peaks:",rational_peaks(2**n_B,r_B,r_B))

    print("\nTask 3 — Toffoli Verification")
    truth=toffoli_truth_table()
    print("Truth table (basis states, in->out, a b c order):")
    for k in sorted(truth.keys()):
        print(f"{k} -> {truth[k]}")
    sup_counts=toffoli_on_superposition()
    print("Superposition test counts (a b c order):",sup_counts)

    print("\nAll outputs saved under ./outputs")
