# LLM Wiki Schema & Operational Guidelines

This document defines the operational rules, directory layout, frontmatter formats, and standard workflows for maintaining this persistent, LLM-driven knowledge base.

---

## 1. Core Principles

1. **Persistent Knowledge Base**: Knowledge is compiled incrementally into interlinked markdown files in `wiki/`, compounding over time rather than being re-derived from scratch on every query.
2. **Three Layers**:
   - `raw/`: Immutable source materials (papers, notes, transcripts, datasets). The LLM reads but never alters these files.
   - `wiki/`: LLM-maintained markdown knowledge base (summaries, concepts, entities, syntheses).
   - `AGENTS.md`: Schema, conventions, and operational workflows.
3. **Obsidian Compatibility**: Use standard Obsidian-style internal links: `[[page-name]]` or `[[page-name|Display Text]]`.
4. **Active Bookkeeping**: Every ingestion updates relevant entity/concept pages, `wiki/index.md`, and appends an entry to `wiki/log.md`.

---

## 2. Directory Structure

```text
.
├── AGENTS.md                  # Wiki Schema & Operational Directives
├── raw/                       # Immutable Raw Sources (PDFs, markdown, text, data)
│   └── assets/                # Downloaded figures, diagrams, and images
└── wiki/                      # Persistent Knowledge Base
    ├── index.md               # Dynamic catalog & overview of all wiki topics/entities
    ├── log.md                 # Append-only chronological log of operations
    ├── concepts/              # Core conceptual synthesis & theoretical foundations
    ├── entities/              # Specific tools, algorithms, databases, genes, diseases
    ├── summaries/             # Detailed per-source summaries with back-references
    └── synthesis/             # Cross-source syntheses, comparative analyses, theses
```

---

## 3. Page Templates & Frontmatter Conventions

All wiki markdown files must start with structured YAML frontmatter.

### Concept / Entity / Synthesis Template
```markdown
---
type: concept | entity | synthesis
title: Page Title
tags: [tag1, tag2]
created: YYYY-MM-DD
updated: YYYY-MM-DD
sources:
  - "[[source-summary-name]]"
aliases: [Alternative Name, Acronym]
---

# Title

## Overview
Brief 2-3 sentence high-level definition and context.

## Key Mechanisms / Properties
Core details, mathematical formulations, or characteristics.

## Biological / Analytical Relevance
Context in bioinformatics, clinical relevance, or computational pipelines.

## Related Concepts & Entities
- [[Related-Concept-1]]
- [[Related-Entity-1]]

## Citations & Sources
- [[source-summary-name]]
```

### Source Summary Template (`wiki/summaries/`)
```markdown
---
type: summary
title: "Title of Paper/Document"
authors: "Author Names"
year: YYYY
raw_source: "raw/filename.ext"
ingested: YYYY-MM-DD
tags: [tag1, tag2]
---

# Title of Document

## Core Thesis / Executive Summary
High-level summary of findings and conclusions.

## Key Takeaways
- Point 1
- Point 2

## Methodological Details
Techniques, algorithms, cohorts, parameters.

## Touched Wiki Pages
- Concept: [[Concept-Name]]
- Entity: [[Entity-Name]]
```

---

## 4. Standard Workflows

### 4.1 Ingestion (`/ingest <source>` or on new raw document)
1. **Read & Extract**: Thoroughly read the raw document from `raw/` or provided input.
2. **Create Summary**: Generate `wiki/summaries/<source-slug>.md` using the summary template.
3. **Cross-Link & Update Entities/Concepts**:
   - Inspect existing pages in `wiki/concepts/` and `wiki/entities/`.
   - Update existing pages with new findings, noting any contradictions or consensus.
   - Create new concept/entity pages if an important subject lacks one.
4. **Update Catalog**: Add new/updated pages with 1-line descriptions to `wiki/index.md`.
5. **Log Operation**: Append entry to `wiki/log.md` with format:
   `## [YYYY-MM-DD] ingest | <Source Title>`

### 4.2 Query & Synthesis
1. **Search & Read**: Consult `wiki/index.md`, read relevant concept/entity/summary pages.
2. **Synthesize Response**: Provide a comprehensive answer with `[[wikilinks]]` and raw citations.
3. **File Back**: If the query produces a substantial novel synthesis, comparison table, or thesis, file it directly as a new page in `wiki/synthesis/<topic-slug>.md` and log it in `wiki/log.md`.

### 4.3 Wiki Lint & Maintenance (`/lint`)
1. **Check Broken Links**: Identify `[[wikilinks]]` pointing to non-existent files.
2. **Identify Orphans**: Locate pages with zero inbound links.
3. **Resolve Conflicts**: Highlight contradicting claims across summaries.
4. **Suggest Extensions**: Propose missing entity/concept pages or targeted web searches for knowledge gaps.
