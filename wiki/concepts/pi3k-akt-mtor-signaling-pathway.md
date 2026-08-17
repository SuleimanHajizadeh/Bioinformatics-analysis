---
type: concept
title: PI3K-Akt-mTOR Signaling Axis in Oncology & Drug Resistance
tags: [pathway, signaling, pi3k, akt, mtor, oncology, tnbc, apoptosis]
created: 2026-08-17
updated: 2026-08-17
sources:
  - "[[tnbc-akt1-alphafold-docking-manuscript]]"
  - "[[statistical-methods-primer]]"
aliases: [PI3K-Akt Pathway, Akt Signaling, mTORC1 Cascade]
---

# PI3K-Akt-mTOR Signaling Axis in Oncology & Drug Resistance

## Overview
The **PI3K-Akt-mTOR pathway** is an intracellular signaling cascade critical for cell survival, nutrient sensing, proliferation, protein synthesis, and metabolic plasticity. In [[triple-negative-breast-cancer]], genetic lesions (*PIK3CA* activating mutations, *PTEN* loss-of-function, *INPP4B* downregulation) cause constitutive pathway hyperactivation, driving therapeutic chemoresistance.

## Molecular Cascade & Enzymatic Mechanism

```mermaid
graph TD
    RTK[Receptor Tyrosine Kinases / EGFR] -->|Activation| PI3K[Class IA PI3K p110/p85]
    PI3K -->|Phosphorylates PIP2 to| PIP3[PIP3 Second Messenger]
    PTEN[PTEN Phosphatase] -.->|Dephosphorylates PIP3 to PIP2| PIP3
    PIP3 -->|Recruits via PH Domains| PDK1[PDK1 Kinase] & AKT[[akt1-kinase]]
    PDK1 -->|Phosphorylates Thr308| AKT
    mTORC2[mTORC2 Complex] -->|Phosphorylates Ser473| AKT
    AKT -->|Inhibits via Phosphorylation| TSC[TSC1/TSC2 Complex]
    TSC -.->|Inactivates Rheb| RHEB[Rheb-GTP]
    RHEB -->|Activates| mTORC1[mTORC1 Complex]
    mTORC1 --> S6K1[p70S6K: Translation] & EBP1[4E-BP1: Cap-dependent Translation]
    AKT -->|Phosphorylates & Inactivates| BAD[BAD: Anti-Apoptosis] & FOXO[FOXO3a: Cell Cycle Arrest] & GSK3[GSK-3beta: Metabolism]
```

## Therapeutic Targets & Inhibitor Classes
1. **Pan-AKT / Isoform-Selective Allosteric Inhibitors**: [[mk-2206]] binds the interface between the PH and kinase domains, preventing membrane recruitment.
2. **ATP-Competitive Orthosteric Inhibitors**: [[ipatasertib]] (GDC-0068) and Capivasertib (AZD5363) bind the catalytic cleft.
3. **Dual PI3K/mTOR Inhibitors**: Dactolisib (BEZ235), Gedatolisib.

## Related Concepts & Entities
- [[akt1-kinase]]
- [[mk-2206]]
- [[ipatasertib]]
- [[triple-negative-breast-cancer]]
- [[transcriptome-to-structure-drug-discovery]]

## Citations & Sources
- [[tnbc-akt1-alphafold-docking-manuscript]]
