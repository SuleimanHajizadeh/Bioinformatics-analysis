---
type: concept
title: AlphaFold Deep Learning Structural Biology & Model Validation
tags: [structural-biology, alphafold2, plddt, pae, protein-structure, biophysics]
created: 2026-08-17
updated: 2026-08-17
sources:
  - "[[tnbc-akt1-alphafold-docking-manuscript]]"
aliases: [AlphaFold Modeling, pLDDT Confidence, PAE Matrix, Structural Prediction]
---

# AlphaFold Deep Learning Structural Biology & Model Validation

## Overview
**AlphaFold** (developed by Google DeepMind) predicts 3D atomic coordinates directly from primary amino acid sequences using an evoformer neural architecture and structural module with invariant point attention (IPA). Accurate utilization in drug discovery requires quantitative model validation using intrinsic confidence metrics.

## Key Confidence & Geometry Metrics

### 1. Predicted Local Distance Difference Test (pLDDT)
Measures per-residue local confidence on a scale of $0\text{--}100$:
- $\text{pLDDT} > 90$: High accuracy; side-chain rotamers reliable for drug docking and pocket analysis.
- $70 \le \text{pLDDT} \le 90$: Good backbone confidence.
- $50 \le \text{pLDDT} < 70$: Low confidence; indicates intrinsic disorder or flexible loop regions.
- $\text{pLDDT} < 50$: Intrinsically Disordered Region (IDR).

### 2. Predicted Aligned Error (PAE)
- $N \times N$ matrix indicating expected positional error (in Ångströms) of residue $y$ when aligned on residue $x$.
- Low cross-domain error indicates rigid inter-domain packing; high error indicates flexible domain-domain movement.

### 3. Stereochemical Validation (Ramachandran Plot)
- Validates backbone dihedral angles ($\phi, \psi$).
- In high-quality models (such as modeled [[akt1-kinase]]), $> 90\%$ of non-glycine residues reside within energetically favored regions.

## Related Concepts & Entities
- [[molecular-docking-virtual-screening]]
- [[akt1-kinase]]
- [[alphafold2]]

## Citations & Sources
- [[tnbc-akt1-alphafold-docking-manuscript]]
