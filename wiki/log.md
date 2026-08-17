# Wiki Activity Log

This is an append-only chronological log of all ingestion, query synthesis, and maintenance operations.

Each entry follows the parseable format:
`## [YYYY-MM-DD] <action> | <Title>`

---

## [2026-08-17] init | LLM Wiki Architecture Setup
- Initialized 3-layer architecture (`AGENTS.md`, `raw/`, `wiki/`).
- Configured folder structure: `wiki/concepts/`, `wiki/entities/`, `wiki/summaries/`, `wiki/synthesis/`, `raw/assets/`.
- Configured Obsidian-compatible wikilinks and YAML frontmatter specifications in [AGENTS.md](file:///Users/macbookairm2/Documents/GitHub/Bioinformatics-analysis/AGENTS.md).

## [2026-08-17] ingest | Statistical Methodology Primer
- **Source**: `docs/STATISTICAL_METHODS.md`
- **Summary**: [[statistical-methods-primer]]
- **Entities/Concepts Touched/Created**:
  - Concept: [[negative-binomial-model]]
  - Concept: [[multiple-testing-fdr]]
  - Concept: [[variance-stabilizing-transformation]]
  - Concept: [[wgcna-coexpression-networks]]
  - Entity: [[deseq2]]
  - Entity: [[triple-negative-breast-cancer]]
  - Synthesis: [[transcriptomic-analysis-framework]]
- **Updated Catalog**: [[index|wiki/index.md]]
