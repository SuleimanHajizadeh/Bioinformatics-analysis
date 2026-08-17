---
type: entity
title: AlphaFold2 Neural Protein Structure Prediction System
tags: [tools, ai, deep-learning, structural-biology, deepmind]
created: 2026-08-17
updated: 2026-08-17
sources:
  - "[[tnbc-akt1-alphafold-docking-manuscript]]"
aliases: [AlphaFold2, AF2, Jumper et al. 2021]
---

# AlphaFold2 Neural Protein Structure Prediction System

## Overview
**AlphaFold2** is an artificial intelligence system developed by Google DeepMind that predicts 3D macromolecular structures from primary amino acid sequences with atomic-level accuracy.

## Architecture
- **Multiple Sequence Alignment (MSA) & Pair Representation**: Jointly processed through 48 Evoformer neural blocks.
- **Structure Module**: Invariant Point Attention (IPA) directly predicts 3D rotation and translation coordinates of protein backbones without template dependency.

## Role in Repository
- Generated high-accuracy structural coordinates for [[akt1-kinase]] (mean pLDDT 92.4) enabling virtual screening and molecular docking.

## Related Concepts & Entities
- [[alphafold-structural-modeling]]
- [[akt1-kinase]]
- [[autodock-vina]]

## Citations & Sources
- [[tnbc-akt1-alphafold-docking-manuscript]]
