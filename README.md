# Recursive PCA for Adaptive Process Monitoring

Implementation of **Li et al. (2000)** — *Recursive PCA for Adaptive Process Monitoring*  
Journal of Process Control, Vol. 10, pp. 471–486

**Course:** Adaptive Learning, Estimation and Supervision of Dynamical Systems  
**A.Y.:** 2025/2026  
**Author:** Amin Borqal — 1073928 — University of Bergamo

---

## Repository Structure

├── DataPreparation.mlx         # Dataset generation and visualization
├── PCA_batch.mlx               # Batch PCA — fixed model monitoring
├── RPCA.mlx                    # RPCA — sample-wise and block-wise Lanczos
├──── lanczos_tridiag.m           # Lanczos tridiagonalization
├──── bisection_eig.m             # Eigenvalue computation via bisection + Sturm
├──── sturm_count.m               # Sturm sequence eigenvalue counter
├──── inverse_iteration.m         # Eigenvector recovery via inverse iteration
├── plots/                      # All exported figures (PNG)
└── LIVESCRIPTS-HTML&PDF/
    ├── DataPreparation.html
    ├── DataPreparation.pdf
    ├── PCA_batch.html
    ├── PCA_batch.pdf
    ├── RPCA.html
    └── RPCA.pdf

---

## How to View Results

### Option A — HTML (recommended)
Download and open in any browser — all plots rendered inline:

1. Go to `LIVESCRIPTS-HTML&PDF/`
2. Download `DataPreparation.html`, `PCA_batch.html`, `RPCA.html`
3. Open in browser

### Option B — PDF
Static version of all results — download from `LIVESCRIPTS-HTML&PDF/`

### Option C — Run in MATLAB (R2022b or later)
Run scripts in order — all `.m` files must be in the same folder:

  DataPreparation.mlx   →  generates monitoring_data.mat
  PCA_batch.mlx         →  requires monitoring_data.mat
  RPCA.mlx              →  requires monitoring_data.mat


---

## Dataset

| Phase | Samples | Description |
|---|---|---|
| Training | $N = 10\,000$ | Stationary healthy — used to fit initial model |
| Monitoring | $N_2 = 1\,000$ | Time-varying — 3 drifts + 3 fault types |

**3 simultaneous drifts** (always active):
- Mean of $t_1$ increases linearly
- Correlation structure: $A_{11}(k) = A_{11}(0) + 0.0002k$
- Sensor aging: $x_2(k) = x_2(k) + 0.002k$

**3 injected fault types:**

| Fault | Samples | Sensor | Type | Magnitude |
|---|---|---|---|---|
| Impulsive | 500–520 | $x_5$ | Additive spike | $+5\sigma_{x_5}$ |
| Incipient | 600–700 | $x_4$ | Sinusoidal growth | $+5\sin(\pi(k-600)/100)$ |
| Step | 800–1000 | $x_1$ | Permanent bias | $+5\sigma_{x_1}$ |

---

## Methods

| Method | Update | Eigenstructure | $q$ selection |
|---|---|---|---|
| Batch PCA | never | SVD once | fixed $q=3$ |
| RPCA Sample-wise | every sample | Lanczos + Bisection + Inverse iteration | VRE adaptive |
| RPCA Block-wise | every $B=10$ samples | Lanczos + Bisection + Inverse iteration | VRE adaptive |

**Forgetting factors tested:** $\lambda \in \{0.980,\, 0.995,\, 0.999\}$  
**Effective memory:** $W_{\text{eff}} = 1/(1-\lambda) \in \{50,\, 200,\, 1000\}$

---

## Metrics

- **FAR** — False Alarm Rate on healthy drift samples
- **DR** — Detection Rate on fault samples
- **Precision, Recall, F1** — per fault type
- **Confusion matrices** — SPE-based, $\lambda = 0.995$

---

## Reference

W. Li, H.H. Yue, S. Valle-Cervantes, S.J. Qin  
*Recursive PCA for Adaptive Process Monitoring*  
Journal of Process Control, Vol. 10, pp. 471–486, 2000  
https://doi.org/10.1016/S0959-1524(00)00022-6
