# Raw Sources Directory

This directory stores immutable raw source materials (papers, markdown notes, datasets, transcripts, PDFs, screenshots).

## Rules
- **Immutable**: Raw sources are read by the LLM agent but are never modified.
- **Organization**:
  - Store text, markdown clippings, and PDFs in `raw/`
  - Store images, figures, and charts in `raw/assets/`
- **Clipping**: When using Obsidian Web Clipper or downloading web articles, place markdown and assets here, then trigger `/ingest`.
