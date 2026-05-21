# Recursive PCA for Adaptive Process Monitoring

Implementation of **Li et al. (2000)** — *Recursive PCA for Adaptive Process Monitoring*  
Journal of Process Control, Vol. 10, pp. 471–486

**Course:** Adaptive Learning, Estimation and Supervision of Dynamical Systems  
**A.Y.:** 2025/2026  
**Author:** Amin Borqal — 1073928 — University of Bergamo

---

## Repository Structure

```
├── DataPreparation.mlx              # Dataset generation and visualization
├── PCA_batch.mlx                    # Batch PCA — fixed model monitoring
├── RPCA.mlx                         # RPCA — sample-wise and block-wise
│   ├── lanczos_tridiag.m            # Lanczos tridiagonalization
│   ├── bisection_eig.m              # Eigenvalue computation via bisection + Sturm
│   ├── sturm_count.m                # Sturm sequence eigenvalue counter
│   └── inverse_iteration.m          # Eigenvector recovery via inverse iteration
├── plots/  # All exported figures (PNG)
└── LIVESCRIPTS-HTML&PDF/
    ├── DataPreparation.html / .pdf
    ├── PCA_batch.html / .pdf
    └── RPCA.html / .pdf
```

---

## How to View Results

### Option A — HTML (recommended)

All plots rendered inline, no MATLAB needed:

1. Go to `LIVESCRIPTS-HTML&PDF/`
2. Download `DataPreparation.html`, `PCA_batch.html`, `RPCA.html`
3. Open in any browser

### Option B — PDF

Static version — download from `LIVESCRIPTS-HTML&PDF/`

### Option C — Run in MATLAB (R2022b or later)

All `.m` helper files must be in the same folder. Run in order:

```
DataPreparation.mlx   →  generates monitoring_data.mat
PCA_batch.mlx         →  requires monitoring_data.mat
RPCA.mlx              →  requires monitoring_data.mat
```

---

## Dataset

| Phase | Samples | Description |
|---|---|---|
| Training | N = 10,000 | Stationary healthy — used to fit batch PCA model |
| Monitoring | N₂ = 1,000 | Time-varying — 3 drifts + 3 fault types |

Data generated from the **exact synthetic model of Li et al. (2000), Equation 50**:

6 measured variables, 3 latent variables (t₁ ~ N(0,1.0), t₂ ~ N(0,0.8), t₃ ~ N(0,0.6)), noise variance 0.04.

**3 simultaneous drifts** (always active during monitoring):

- Mean of t₁ increases linearly
- Correlation structure: A₁₁(k) = A₁₁(0) + 0.0002k
- Sensor aging: x₂(k) = x₂(k) + 0.002k

**3 injected fault types:**

| Fault | Samples | Sensor | Type | Magnitude |
|---|---|---|---|---|
| Impulsive | 500–520 | x₅ | Additive spike | +5σ(x₅) |
| Incipient | 600–700 | x₄ | Sinusoidal growth | +5sin(π(k-600)/100) |
| Step | 800–1000 | x₁ | Permanent bias | +5σ(x₁) |

![Monitoring signals overview](plots/monitoring_signals.png)

---

## Part 1 — The Problem: Batch PCA on Time-Varying Processes

Standard Batch PCA learns a **fixed** correlation matrix from training data. When the process drifts, the model becomes obsolete.

### Scree plot — selecting q = 3

PC1: 49.2% · PC2: 26.0% · PC3: 16.8% · **Total: 92.0%**

![Scree plot](plots/scree_plot.png)


### Batch PCA: healthy stationary data ✅

Alarm rates at nominal 5% on stationary training distribution:

![Batch PCA stationary](plots/pca_4a_stationary.png)

### Batch PCA: drift only (no fault) ❌

With time-varying drift and **no real fault**, T² alarm rate reaches **32.6%** — 6.5× the expected rate:

![Batch PCA drift failure](plots/pca_4b_drift.png)

**Diagnosis:** Model becomes obsolete as process drifts → fixed loadings misrepresent current correlation structure.

### Batch PCA: fault detection

**Impulsive fault (x₅, k=500–540):**

![Batch PCA fault 1](plots/pca_4c_fault1.png)

**Incipient fault (x₄, k=600–700):**

![Batch PCA fault 2](plots/pca_4d_fault2.png)

**Step fault (x₁, k=800–1000):**

![Batch PCA fault 3](plots/pca_4e_fault3.png)

**Full monitoring run (drift + all 3 faults):**

![Batch PCA full monitoring](plots/pca_4f_full.png)

**SPE contribution plot — root cause localization:**

x₁ dominates with ~57% of SPE at post-fault samples. Correctly identifies the injected bias source.

![Batch PCA contributions](plots/pca_4g_contributions.png)

---

## Part 2 — Solution: Recursive PCA (RPCA)

RPCA updates the covariance estimate recursively using an **Exponentially Weighted Moving Average**:

**R**_{k+1} = λ **R**_k + **v**₁**v**₁ᵀ + **v**₂**v**₂ᵀ

where **v**₁ captures the effect of mean change and **v**₂ = √(1−λ) · **x̃** incorporates the new standardized sample.

**Forgetting factor** λ ∈ (0,1) controls effective memory:

| λ | W_eff = 1/(1−λ) | Behavior |
|---|---|---|
| 0.980 | 50 samples | Fast adaptation — risk of adapting to faults |
| 0.995 | 200 samples | Medium — balanced |
| 0.999 | 1000 samples | Long memory — reliable sustained detection |

### Eigenstructure update: Lanczos algorithm

Rather than full SVD (O(m³)), RPCA uses **Lanczos tridiagonalization**:

1. Tridiagonalize **R**_k → **Γ** in O(4m²m₁) operations
2. Compute q largest eigenvalues via **bisection + Sturm sequences**
3. Recover eigenvectors via **inverse iteration**

Termination: stop when trace(**Γ**) / trace(**R**_k) > 0.995 (99.5% variance captured).

### Adaptive q via VRE method

VRE(q) = Σᵢ uᵢ / var(xᵢ)

where uᵢ = reconstruction error variance for variable i. U-shaped curve: minimum at optimal q.

**VRE snapshot at k=250 (mid-monitoring, pre-fault):**

![VRE snapshot](plots/rpca_4b_vre_snapshot.png)

**Adaptive q evolution — sample-wise, all three λ:**

![Adaptive q evolution SW](plots/rpca_4a_q_evolution.png)

- λ=0.980: q oscillates 1–3, hyper-responsive to every perturbation
- λ=0.995: q transitions during drift, recovers post-fault  
- λ=0.999: q stable at 3, long memory stabilizes dimensionality

---

## Part 3 — RPCA Results: Sample-Wise Update

### Stationary baseline

![RPCA stationary](plots/rpca_3a_stationary.png)

### Drift only (no fault) — RPCA tracks successfully

![RPCA drift](plots/rpca_3b_drift.png)

RPCA adapts to process drift → alarm rates remain controlled, unlike Batch PCA.

### Fault detection

**Impulsive fault (x₅, k=500–520):**

![RPCA fault 1](plots/rpca_3c_fault1.png)

**Incipient fault (x₄, k=600–700):**

![RPCA fault 2](plots/rpca_3d_fault2.png)

**Step fault (x₁, k=800–1000):**

![RPCA fault 3](plots/rpca_3e_fault3.png)

**Full monitoring run:**

![RPCA full monitoring](plots/rpca_3f_full.png)

---

## Part 4 — RPCA Results: Block-Wise Update (B = 10)

Block-wise RPCA updates every **B = 10 samples**, reducing real-time compute overhead.

**Drift only:**

![Block-wise drift](plots/rpca_6a_bw_drift.png)

**Impulsive fault:**

![Block-wise fault 1](plots/rpca_6b_bw_fault1.png)

**Incipient fault:**

![Block-wise fault 2](PresentationRPCA/plots/rpca_6c_bw_fault2.png)

**Step fault:**

![Block-wise fault 3](plots/rpca_6d_bw_fault3.png)

**Full run:**

![Block-wise full monitoring](plots/rpca_6e_bw_full.png)

**Adaptive q evolution (block-wise):**

![Block-wise q evolution](plots/rpca_7_bw_q_evolution.png)

---

## Part 5 — Sample-Wise vs Block-Wise Comparison

Direct head-to-head comparison across λ values:

![SW vs BW — 1](plots/rpca_8_1_sw_vs_bw.png)

![SW vs BW — 2](plots/rpca_8_2_sw_vs_bw.png)

![SW vs BW — 3](plots/rpca_8_3_sw_vs_bw.png)

![SW vs BW — 4](plots/rpca_8_4_sw_vs_bw.png)

![SW vs BW — 5](plots/rpca_8_5_sw_vs_bw.png)

**SPE contribution decomposition (RPCA, post-fault sample k=350):**

![RPCA contributions](plots/rpca_8b_contributions.png)

x₁ dominates at ~57% of SPE — correctly identifies the root cause of the step fault.

---

## Part 6 — Quantitative Evaluation (λ = 0.995)

### Confusion matrices — SPE-based detection, per fault type

![Confusion matrices](plots/rpca_9a_confusion_matrices.png)

### F1 score: Batch PCA vs RPCA sample-wise vs RPCA block-wise

![F1 comparison](plots/rpca_9b_f1_comparison.png)

### FAR (False Alarm Rate) and DR (Detection Rate) summary

![FAR and DR](plots/rpca_9d_far_dr.png)

---

## Methods Summary

| Method | Update | Eigenstructure | q selection |
|---|---|---|---|
| Batch PCA | never | full SVD once | fixed q = 3 |
| RPCA sample-wise | every sample | Lanczos + Bisection + Inverse iteration | VRE adaptive |
| RPCA block-wise | every B = 10 samples | Lanczos + Bisection + Inverse iteration | VRE adaptive |

**Forgetting factors tested:** λ ∈ {0.980, 0.995, 0.999}  
**Effective memory:** W_eff ∈ {50, 200, 1000} samples

---

## Key Results

| Method | FAR (drift only) | DR (step fault) | Verdict |
|---|---|---|---|
| Batch PCA | **32.6% T²** (6.5× nominal) | high but unreliable | Model obsolete during drift ❌ |
| RPCA λ=0.980 | ~5% | adapts to fault → loses it | Dangerous for sustained faults ❌ |
| RPCA λ=0.995 | controlled | moderate | Acceptable ⚠️ |
| RPCA λ=0.999 | controlled (~12%) | **81%+ sustained** | Best trade-off ✅ |

**Core insight:** λ=0.999 (W_eff=1000) prevents the model from absorbing the fault as new normal — sustained, reliable detection at ~2–4× nominal false alarm rate during fast drift phases. For safety-critical applications: reliability > raw sensitivity.

---

## Reference

W. Li, H.H. Yue, S. Valle-Cervantes, S.J. Qin  
*Recursive PCA for Adaptive Process Monitoring*  
Journal of Process Control, Vol. 10, pp. 471–486, 2000  
https://doi.org/10.1016/S0959-1524(00)00022-6
