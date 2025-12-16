## Pseudocode: Cross-Modal Trajectory Alignment via Orthogonal Procrustes

### Inputs
- Data matrix \( D \in \mathbb{R}^{n \times p} \)
- Cell labels \( \ell \in \{1,\dots,K\}^n \)
- Number of PCs \( d \)
- Target number of trajectory points \( T \)

---

### Step 1: Construct Two Related Datasets
**Goal:** Create two datasets with controlled overlap.

1. Partition samples:
   \[
   \mathcal{I}_A \cup \mathcal{I}_B = \{1,\dots,n\}, \quad
   \mathcal{I}_A \cap \mathcal{I}_B = \varnothing
   \]

2. Partition features:
   \[
   \mathcal{J}_A \cup \mathcal{J}_B = \{1,\dots,p\}, \quad
   \mathcal{J}_A \cap \mathcal{J}_B = \varnothing
   \]

3. Define datasets:
   \[
   A = D[\mathcal{I}_A, \mathcal{J}_A], \quad
   B = D[\mathcal{I}_B, \mathcal{J}_B]
   \]

---

### Step 2: Standardize Features
For each dataset \( X \in \{A,B\} \) and each feature \( j \):
\[
X_{ij} \leftarrow \frac{X_{ij} - \mu_j}{\sigma_j}
\]

---

### Step 3: Dimensionality Reduction (PCA)
For each standardized dataset \( X \):

1. Compute PCA:
   \[
   Z = X W_d \in \mathbb{R}^{n_X \times d}
   \]

2. Store low-dimensional embedding \( Z \)

---

### Step 4: Trajectory Inference
**Goal:** Estimate a smooth principal curve in latent space.

For each embedding \( Z \):

1. Fit lineage-constrained trajectory model:
   \[
   \gamma(t) \subset \mathbb{R}^d
   \]

2. Extract ordered curve points:
   \[
   C = \{\gamma(t_1), \dots, \gamma(t_m)\}
   \]

---

### Step 5: Arc-Length Resampling
**Goal:** Ensure trajectories have equal length.

Given curve \( C = \{c_1,\dots,c_m\} \):

1. Compute cumulative arc length:
   \[
   s_1 = 0, \quad
   s_k = \sum_{i=2}^k \|c_i - c_{i-1}\|
   \]

2. Normalize:
   \[
   \tilde{s}_k = s_k / s_m
   \]

3. Interpolate curve at:
   \[
   \tilde{s}^\star = \left\{ \frac{0}{T-1}, \dots, \frac{T-1}{T-1} \right\}
   \]

4. Obtain resampled curves:
   \[
   X, Y \in \mathbb{R}^{T \times d}
   \]

---

### Step 6: Centering and Normalization
1. Compute means:
   \[
   \bar{X}, \bar{Y}
   \]

2. Center:
   \[
   X \leftarrow X - \bar{X}, \quad
   Y \leftarrow Y - \bar{Y}
   \]

3. (Optional) Normalize Frobenius norm:
   \[
   X \leftarrow \frac{X}{\|X\|_F}, \quad
   Y \leftarrow \frac{Y}{\|Y\|_F}
   \]
Where 
   \[
    \|A\|_F 
    \;=\; 
    \sqrt {tr(A\top A)}
   \]
or equivalently
    \[
    \|A\|_F
    \;=\;
    \sum_{i=1}^{m}\sum_{j=1}^{n} A_{ij}^2.
    \]
This essentially unit normalizes the cummulative variance as opposed to normalizing each column independently. 

---

### Step 7: Orthogonal Procrustes Alignment
**Goal:** Align \( Y \) to \( X \).

1. Cross-covariance:
   \[
   M = Y^\top X
   \]

2. Singular value decomposition:
   \[
   M = U \Sigma V^\top
   \]

3. Optimal rotation:
   \[
   R^\star = U V^\top
   \]

4. (Optional) Scaling:
   \[
   c = \frac{\operatorname{tr}(\Sigma)}{\|X\|_F^2}
   \]

5. Apply transformation:
   \[
   Y^\star = c Y R^\star
   \]

6. Translation:
   \[
   b = \bar{X} - c \bar{Y} R^\star
   \]

---

### Step 8: Alignment Error
Residual sum of squares:
\[
\mathrm{RSS}
= \|X\|_F^2 + c^2 \|Y\|_F^2 - 2c \operatorname{tr}(\Sigma)
\]

---

### Results
- Aligned trajectory \( Y^\star + b \)
- Rotation \( R^\star \)
- Scale \( c \)
- Translation \( b \)
- Alignment error (RSS)

