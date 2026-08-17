---
type: concept
title: Molecular Docking, Binding Affinity & In Silico Drug Discovery
tags: [drug-discovery, molecular-docking, autodock-vina, biophysics, binding-affinity]
created: 2026-08-17
updated: 2026-08-17
sources:
  - "[[tnbc-akt1-alphafold-docking-manuscript]]"
aliases: [Molecular Docking, Virtual Screening, AutoDock Vina, Free Energy of Binding]
---

# Molecular Docking, Binding Affinity & In Silico Drug Discovery

## Overview
**Molecular Docking** computationally models the quaternary interaction of a small molecule (ligand) within the binding pocket of a target macromolecule (receptor protein). Scoring functions estimate the Gibbs free energy of binding ($\Delta G_{\text{bind}}$) in $\text{kcal/mol}$.

## Physics & Mechanics in AutoDock Vina

### 1. Scoring Function & Free Energy
$$\Delta G_{\text{bind}} = \Delta G_{\text{vdW}} + \Delta G_{\text{hbond}} + \Delta G_{\text{elec}} + \Delta G_{\text{desolv}} + \Delta G_{\text{tors}}$$
- **Negative $\Delta G$**: Exergonic, spontaneous binding. Values more negative than $-8.0\text{ kcal/mol}$ typically indicate nanomolar/sub-micromolar affinity candidates.

### 2. Allosteric vs. Orthosteric (ATP-Competitive) Inhibition
- **ATP-Competitive (Orthosteric)**: Binds the catalytic ATP pocket (e.g., Ipatasertib). Often suffers from off-target toxicity due to high sequence conservation across the AGC kinase family.
- **Allosteric Inhibition**: Binds a distant regulatory cavity (e.g., the cavity between the PH domain and catalytic kinase domain in [[akt1-kinase]] targeted by [[mk-2206]]). Locks the enzyme in an inactive, "PH-in" closed conformation.

## Analytical Application
- In the TNBC manuscript, [[mk-2206]] achieved $\Delta G = -9.4 \text{ kcal/mol}$ with key stabilizing hydrogen bonds to Trp80 and Asp292.

## Related Concepts & Entities
- [[akt1-kinase]]
- [[mk-2206]]
- [[autodock-vina]]
- [[alphafold-structural-modeling]]

## Citations & Sources
- [[tnbc-akt1-alphafold-docking-manuscript]]
