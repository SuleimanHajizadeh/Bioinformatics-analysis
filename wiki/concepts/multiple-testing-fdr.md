---
type: concept
title: Multiple Testing Correction and False Discovery Rate (FDR)
tags: [statistics, fdr, benjamini-hochberg, hypothesis-testing]
created: 2026-08-17
updated: 2026-08-17
sources:
  - "[[statistical-methods-primer]]"
aliases: [FDR, Benjamini-Hochberg, BH Procedure, Multiple Hypothesis Testing]
---

# Multiple Testing Correction and False Discovery Rate (FDR)

## Overview
In genomic experiments assessing tens of thousands of genes simultaneously, testing each hypothesis at a naive threshold ($\alpha = 0.05$) yields hundreds to thousands of false positives. **False Discovery Rate (FDR)** control methods adjust statistical significance to constrain the expected proportion of false positives among rejected hypotheses.

## Key Mechanisms & Mathematical Formulation

### The Multiple Testing Dilemma
For $m = 20,000$ independent tests at nominal $\alpha = 0.05$:
$$\text{Expected False Positives} = m \cdot \alpha = 20,000 \times 0.05 = 1,000$$

### FWER vs. FDR
- **Family-Wise Error Rate (FWER - Bonferroni)**: Controls $P(\text{False Positives} \ge 1) \le \alpha$. Significance threshold is $\alpha_{\text{Bonf}} = \alpha / m$. Highly conservative; introduces high false negative rates in exploratory omics.
- **False Discovery Rate (FDR - Benjamini-Hochberg)**: Controls $\text{E}[V/R]$ where $V$ is false discoveries and $R$ is total rejected null hypotheses.

### Benjamini-Hochberg (BH) Step-Up Algorithm
1. Sort raw p-values: $P_{(1)} \le P_{(2)} \le \dots \le P_{(m)}$.
2. Find largest index $k$ such that:
   $$P_{(k)} \le \frac{k}{m} q$$
3. Reject $H_0$ for all $i \le k$.
4. Adjusted p-value (q-value):
   $$q_i = \min \left( \min_{j \ge i} \left\{ \frac{m \cdot P_{(j)}}{j} \right\}, 1 \right)$$

## Biological & Analytical Relevance
- Standard requirement for differential expression calls in [[deseq2]].
- Threshold typically set to $q < 0.05$ or $q < 0.01$.

## Related Concepts & Entities
- [[negative-binomial-model]]
- [[deseq2]]
- [[transcriptomic-analysis-framework]]

## Citations & Sources
- [[statistical-methods-primer]]
