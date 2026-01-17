# ufekp – A Unified Framework for Efficient Kernel and Polynomial Interpolation

This repository contains the MATLAB implementation that accompanies the paper

> **M. Belianovich, G. E. Fasshauer, A. Narayan, V. Shankar**  
> *A Unified Framework for Efficient Kernel and Polynomial Interpolation* (2025)  
> arXiv:2507.12629

The code implements the unified interpolation scheme that combines compactly-supported positive-definite kernels with multivariate polynomials. It reproduces the 1D/2D/3D Euclidean experiments and the manifold experiments (sphere, torus, and hemisphere), including the ablation and PHS+poly comparisons shown in the paper.

---

## 1. Repository layout

The current layout is:

```text
ufekp/
  code/
    core/            % core numerical linear algebra + interpolation routines
    utils/           % node generators, target functions, helper utilities
  data/
    euclidean/       % data used to produce results on the euclidean domain
    manifolds/       % data used to produce results on the manifolds
  experiments/
    1D/              % drivers for 1D tests on [-1,1]
    2D/              % drivers for 2D disk tests
    2D ablation/     % 2D edge vs interior ablation
    3D/              % drivers for 3D ball tests
    All Dimensions/  % drivers for all dimensions (1D, 2D, 3D) tests
    manifolds/       % drivers for sphere/torus/boundary manifold tests
    plotting/        % plotting scripts that generate the paper-style figures
  paper/             % includes the current pdf of the paper
  results/           % contains all results included and not included in the paper
  README.md
```

---

## 2. Getting Started

### Dependencies
- MATLAB R2020b+  
- `export_fig` installed on path.
- make sure that you are in the repository while using MATLAB, and are not out of the folder, especially to reproduce plots.

### Setup
```matlab
addpath(genpath('code'));
addpath(genpath('experiments'));
addpath(genpath('data'));
```

### Quick Test
#### 1D convergence test
```matlab
cd experiments/1D
Test1D
TestwPHS1D
```

#### 2D disk test 
```matlab
cd experiments/2D
Test2D
TestwPHS2D
```

#### 3D ball test 
```matlab
cd experiments/3D
Test3D
TestwPHS3D
```

#### manifolds
```matlab
cd experiments/manifolds
ManifoldWithBoundaryTest
ManifoldWithBoundaryTestwPHS
SphereTest
SphereTestPHS
TorusTest
TorusTestPHS
```

Plotting can be reproduced using plotting codes in the `experiments/plotting` folder. All the collected results are saved under appropriate function names in folder `results/`.

---

## 3. Core numerical routines

The unified interpolation framework is implemented in the following core files (currently in `Approximation/CleanLagrangeApprox/` and `Approximation/`).

### 3.1 Unified solvers (CSRBF + polynomial)

- **CSRBFGen.m**  
  Main unified solver on Euclidean domains. Given nodes, function values, a CSRBF choice, polynomial degree, and a target condition number \(K_t\), it:
  1. Builds the (sparse) CSRBF Gram matrix using `DistanceMatrixCSRBFwt.m` and `Wendland.m`.
  2. Computes a Cholesky factorization \(A = LL^\top\).
  3. Builds the polynomial matrix \(P\) (Legendre/Jacobi-type basis).
  4. Forms \(B = L^{-1}P\) and applies a rank-revealing QR factorization (via `rand_qr.m`).
  5. Solves the unified block system for kernel and polynomial coefficients.
  6. Returns errors, condition numbers, and timing data.

- **CSRBFGenManifold.m**  
  Unified solver on embedded manifolds (sphere, torus, hemisphere). Uses the same algebra as `CSRBFGen.m`, but nodes lie on a manifold and distances may be warped or intrinsic.

- **CSRBFDiag.m**  
  Implements the “Diag” limit of the unified framework: a polynomial least-squares method obtained when the kernel block is negligible compared to the polynomial block. Used as the **Diag** baseline.

### 3.2 Polynomial baselines and PHS baselines

- **PLS.m**, **PLSManifold.m**  
  Pure polynomial least-squares in Euclidean domains and on manifolds, respectively. Used as the **PLS** baselines.

- **PHS_poly.m**, **PHS_poly2.m**, **PHS_poly_man.m**  
  Polyharmonic spline + polynomial (**PHS+poly** baseline) implementations used in PHS comparison experiments, in both Euclidean and manifold settings.

### 3.3 Kernels, distance matrices, and conditioning tools

- **Wendland.m**, **WendlandShifter.m**, **WendlandShiftedDers.m**  
  Evaluate Wendland CSRBFs and their shifted/scaled variants and derivatives.

- **DistanceMatrixCSRBF.m**, **DistanceMatrixCSRBFwt.m**, **DistanceMatrixCSRBFwtwarped.m**  
  Construct CSRBF distance/Gram matrices (plain, weighted, or warped). These are the only routines that directly form the kernel matrices.

- **RunCond.m**, **rbf_effective_cond.m**, **eval_condition_numbers.m**, **ratio_eval.m**  
  Utilities for:
  - searching for shape/support parameters with \(\mathrm{cond}(A) \approx K_t\),
  - computing effective condition numbers,
  - comparing error vs cost across methods.

### 3.4 Polynomial / Jacobi machinery

The following files provide the building blocks for polynomial spaces and their derivatives:

- **chebyshev_eval.m**, **chebyshev_recurrence.m**  
- **jacobi_recurrence.m**  
- **jacobi_connection_coeffs.m**, **jacobi_iterated_connection_coeffs.m**  
- **jacobi_derivative_coefficients.m**, **jacobi_sparse_diffmat_inv.m**  
- **mpoly_eval.m**, **poly_eval.m**  
- **total_degree_indices.m**, **hyperbolic_cross_indices.m**  
- **mjacobi_symm_faster_differentiation.m**, **mindex_derivative_analysis.m**  
- **rand_qr.m**  

They are used to:
- construct multi-index sets for polynomial spaces,
- evaluate multivariate Jacobi/Legendre-type polynomial bases,
- differentiate polynomials,
- perform stable rank-revealing factorizations.

In the **target layout**, all of these routines reside under:

```text
code/core/
```

---

## 4. Utilities and geometry

Utility functions for node generation, target functions, and diagnostics are currently located in `code/utils/`.

### 4.1 Geometry and node generation

- **spiral_points.m**  
  Generates quasi-uniform nodes on the sphere using a spiral-type mapping.

- **computeHexagonTorusNodes.m**  
  Generates hexagonally-structured nodes on the torus.

- **hemispherePts.m**  
  Generates node sets on a hemispherical patch.

- **half_torus_points.m**  
  Generates node sets on a half-torus.

- **cylinder_points.m**  
  Generates node sets on a cylindrical surface.

These functions are used by the manifold drivers (SphereTest.m, TorusTest.m, ManifoldWithBoundaryTest.m) to construct the node sets used in the manifold experiments.

### 4.2 Target functions and smooth bumps

- **bump_cinf.m**  
  Defines a compactly supported \(C^{\inf} \) bump function used in the manifold-with-boundary experiments and some Euclidean tests.

### 4.3 Miscellaneous utilities

- **compareSparsity.m**  
  Compares sparsity patterns and fill-in of the different systems (PLS, Diag, unified FS) for given node sets and polynomial degrees.

---

## 5. Data: nodes, results, and where they are used

### 5.1 Euclidean data

All Euclidean node sets are stored under `data/euclidean/` folder.

### 5.2 Manifold data

Manifold node sets are likewise currently in `data/manifolds/` folder.

### 5.3 Result files

All drivers write their outputs into `/results/<function_name>/<regime>/results_<function_name>.mat`, where

- function_name is a label like abs_1d, exp_1d, exp_p_2d, exp_p_3d, boundary_test, etc., and

- regime is one of high, low, fixed, mixed (corresponding to the different point/degree regimes in the paper).

Each results_*.mat file typically contains:

- **dim** – spatial dimension (1, 2, or 3)

- **sN** – \( N^{\frac{1}{d}} \) used on the x-axis,

- **el2_poly** – relative \(\ell^2\)error for PLS,

- **el2_diag** – relative \(\ell^2\)error for Diag,

- **el2_fs1**, **el2_fs2**, **el2_fs3** – FS errors for target condition numbers \(K_t = 10^{12}, 10^{8}, 10^{4}\),

- **el2_fs_eff**, **el2_fs_eval** – FS variants tuned by effective condition and evaluation cost,

- (1D only) **el2_fc1**, **el2_fc2**, **el2_fc3** – fixed condition (FC) variants,

- (PHS comparisons) **el2_phs** and timing fields a_time_*, e_time_*.

---

## 6. Experiments and drivers

The drivers that reproduce the paper’s experiments are located in `experiments/`.

### 5.1 1D experiments

- **Test1D.m**

- - Runs 1D tests on \([-1, 1]\) for a chosen target function (abs_1d, rk_1d).

- - Compares PLS, Diag, FS, and FC variants.

- - Writes `results_<name>.mat` into the appropriate `results/<name>/...` folder.

- **TestwPHS1D.m**

- - Runs 1D tests on \([-1, 1]\) for a chosen target function (abs_1d, rk_1d).

- - Compares PLS, Diag, FS, FC, and PHS+poly variants.

- - Writes `results_<name>.mat` into the appropriate `results/<name>/...` folder.

### 5.2 2D experiments

- **Test2D.m**

- - Runs 2D disk experiments for a chosen target function (xy_p_2d, exp_p_2d).

- - Compares PLS, Diag, and FS variants.

- - Writes `results_<name>.mat` into the appropriate `results/<name>/...` folder.

- **TestwPHS2D.m**

- - Runs 2D disk experiments for a chosen target function (xy_p_2d, exp_p_2d).

- - Compares PLS, Diag, FS, and PHS+poly variants.

- - Writes `results_<name>.mat` into the appropriate `results/<name>/...` folder.

#### 2D ablation

- **DriverRun2DAbl.m**, **Ablcond2D.m**

- - Perform 2D edge vs interior ablation experiments.

- - Vary how nodes are concentrated near the boundary and record changes in error and conditioning.

- - Store output in dedicated ablation result folders.

These ablation results are later visualized via `plot_edgeabl.m` and `plot_edgeabl_all.m`.

### 5.3 3D experiments

- **Test3D.m**

- - Runs 3D ball experiments for a chosen target function (xy_p_3d, exp_p_3d).

- - Compares PLS, Diag, and FS variants.

- - Writes `results_<name>.mat` into the appropriate `results/<name>/...` folder.

- **TestwPHS3D.m**

- - Runs 3D ball experiments for a chosen target function (xy_p_3d, exp_p_3d).

- - Compares PLS, Diag, FS, and PHS+poly variants.

- - Writes `results_<name>.mat` into the appropriate `results/<name>/...` folder.

### 5.4 All dimensions experiments

- **TestAll.m**

- - Runs experiments for a chosen dimension (1, 2, 3), and a chosen target function for that dimension.

- - Compares PLS, Diag, FS, and FC variants.

- - Writes `results_<name>.mat` into the appropriate `results/<name>/...` folder.

- **TestwPHS.m**

- - Runs 3D ball experiments for a chosen target function (xy_p_3d, exp_p_3d).

- - Compares PLS, Diag, FS, FC, and PHS+poly variants.

- - Writes `results_<name>.mat` into the appropriate `results/<name>/...` folder.

### 5.5 Manifold experiments

- **SphereTest.m**, **SphereTestPHS.m**

- - Unified FS and PHS+poly comparisons on the sphere.

- - Writes `results_<name>.mat` into the appropriate `results/<name>/...` folder.

- **TorusTest.m**, **TorusTestPHS.m**

- - Unified FS and PHS+poly comparisons on the torus.

- - Writes `results_<name>.mat` into the appropriate `results/<name>/...` folder.

- **ManifoldWithBoundaryTest.m**, **ManifoldWithBoundaryTestPHS.m**

- - Unified FS and PHS+poly comparisons on the hemisphere.

- - Writes `results_<name>.mat` into the appropriate `results/<name>/...` folder.

--- 

## 7. Plotting and figure generation

All plotting scripts are located in the `experiments/plotting/` directory. These scripts read precomputed results stored under the `results/` (or `results_ablation/`) directory and reproduce the figures shown in the paper.

### 7.1 General conventions

Unless stated otherwise, all plotting routines expect results stored in the canonical format

results/<function_name>/<regime>/results_<function_name>.mat

where

- function_name is a string identifier such as `abs_1d`, `exp_p_2d`, `xy_p_3d`, `boundary_test`, `sphere`, or `torus`;
- regime is one of `high`, `low`, `fixed`, or `mixed`, corresponding to the point/degree regimes used in the paper.

Each results_*.mat file typically contains:
- dim – spatial dimension,
- sN – \(N^{\frac{1}{d}}\), used as the x-axis in Euclidean plots,
- error arrays such as el2_poly, el2_diag, el2_fs1, el2_fs2, el2_fs3,
- optional variants such as el2_fc*, el2_phs,
- timing data (a_time_*, e_time_*) for efficiency plots.

All plotting scripts assume that the repository root is on the MATLAB path and that MATLAB is launched from within the repository.

---

### 7.2 Euclidean error plots

plot_results_l2.m

Usage:
plot_results_l2(function_name, regime, base_results_dir, smoothness_idx, variety)

- Reproduces relative \(\ell^2\) error vs. \(N^{\frac{1}{d}}\) plots for Euclidean domains.
- Compares polynomial least squares (PLS), diagonal limit (Diag), and unified (FS) methods.
- For 1D (or when variety=true), multiple FS curves corresponding to different target condition numbers
  \(K_t = 10^12, 10^8, 10^4\) are shown.

Arguments:
- function_name      : string or cell array of function identifiers
- regime             : typically 'high'
- base_results_dir   : optional, default 'results'
- smoothness_idx     : optional kernel smoothness index (clamped to {1,2,3})
- variety            : logical flag controlling multi-panel vs. combined plots

Expected input:
results/<function_name>/<regime>/results_<function_name>.mat

---

plot_results_l2_wphs.m

Usage:
plot_results_l2_wphs(function_name, regime, base_results_dir, smoothness_idx, variety)

- Same as plot_results_l2, but additionally overlays PHS+poly baselines if present.
- The PHS curve is read from the same results file if the field el2_phs exists.

---

### 7.3 Efficiency (error vs. cost) plots

plot_efficiency.m

Usage:
plot_efficiency(function_name, regime, base_results_dir, smoothness_idx, variety)

- Produces error-vs-cost plots comparing PLS, Diag, FS (and FC in 1D).
- Uses timing data stored in the results files.
- With variety=true, produces separate panels for different target condition numbers.

Expected input:
results/<function_name>/<regime>/results_<function_name>.mat

---

plot_efficiency_wphs.m

Usage:
plot_efficiency_wphs(function_name, regime, base_results_dir, smoothness_idx, variety)

- Same as plot_efficiency, but includes PHS+poly curves if available.

---

### 7.4 Manifold error plots

plot_l2_manifold.m

Usage:
plot_l2_manifold(regime, base_results_dir)

- Reproduces relative \(\ell^2\) error plots for manifold experiments:
  boundary_test (hemisphere with boundary),
  sphere,
  torus.

Expected input:
results/<manifold_name>/<regime>/results_<manifold_name>*.mat

The base results file must contain:
- either sN or Ns (node count),
- el2 (FS error),
- el2_poly (PLS baseline).

If a file containing el2_phs is present in the same folder, the PHS+poly baseline is added automatically.

---

plot_sphere_with_phs_sweep.m

Usage:
plot_sphere_with_phs_sweep(regime, base_results_dir)

- Specialized plot for the sphere showing a sweep over PHS degrees.
- Reads FS and PLS data from standard sphere results files and overlays all files matching

results/sphere/<regime>/results_sphere_phsp_*.mat

Each sweep file must contain el2_phs.

---

### 7.5 2D edge vs. interior ablation plots

plot_edgeabl.m
plot_edgeabl_all.m

These scripts visualize the 2D edge-vs-interior ablation experiments.

Expected directory structure:
results_ablation/<function_name>/edge_vs_interior/abl2d_<function_name>_*.mat

Important note:
These plotting scripts load the node file `DiskPoissonNodesClustered.mat` by filename only.
This file must therefore be either:
- in the current working directory, or
- available on the MATLAB path.

Failure to do so will result in a load error.

--- 