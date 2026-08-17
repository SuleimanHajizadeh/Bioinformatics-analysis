---
type: concept
title: Negative Binomial Distribution in RNA-Seq
tags: [statistics, transcriptomics, count-data, overdispersion]
created: 2026-08-17
updated: 2026-08-17
sources:
  - "[[statistical-methods-primer]]"
aliases: [NB Model, Gamma-Poisson, Overdispersion Modeling]
---

# Negative Binomial Distribution in RNA-Seq

## Overview
The **Negative Binomial (NB)** distribution (or Gamma-Poisson mixture) is the standard statistical distribution used to model raw sequencing read counts per gene in high-throughput RNA-Seq experiments.

## Key Mechanisms & Mathematical Formulation

### Poisson Limitation & Overdispersion
Raw sequencing counts $k \in \{0, 1, 2, \dots\}$ were initially modeled with a Poisson distribution:
$$P(X = k) = \frac{\lambda^k e^{-\lambda}}{k!}, \quad \text{E}[X] = \text{Var}(X) = \lambda$$

However, biological replicates exhibit non-identical expression rates, violating the equidispersion property ($\text{Var}(X) > \text{E}[X]$).

### Compound Gamma-Poisson Formulation
By allowing the Poisson mean $\lambda$ to follow a Gamma distribution $\text{Gamma}(\alpha, \beta)$, the marginal distribution of counts becomes Negative Binomial:
$$X \sim \text{NB}(\mu, \alpha)$$
$$\text{Var}(X) = \mu + \alpha \mu^2$$

where:
- $\mu$: Mean expression level across replicates.
- $\alpha$: Biological dispersion parameter ($BCV = \sqrt{\alpha}$).

As $\alpha \to 0$, the NB distribution collapses back to the Poisson distribution.

## Biological / Analytical Relevance
- Accurately captures technical variability (Poisson noise) and biological heterogeneity (Gamma spread).
- Implemented in major differential expression packages including [[deseq2]], `edgeR`, and `limma-voom`.

## Related Concepts & Entities
- [[multiple-testing-fdr]]
- [[variance-stabilizing-transformation]]
- [[deseq2]]

## Citations & Sources
- [[statistical-methods-primer]]
