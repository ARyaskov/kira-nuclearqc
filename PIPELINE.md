# Pipeline Overview

This document summarizes the deterministic pipeline stages implemented in `kira-nuclearqc`.

## Stage 1: Input Discovery and Gene Index
- Discover `matrix.mtx(.gz)`, `features.tsv(.gz)` or `genes.tsv`, `barcodes.tsv(.gz)`
- In `--run-mode pipeline`, deterministically resolve shared cache name and prefer reading `kira-organelle.bin` / `<PREFIX>.kira-organelle.bin`
- Parse features and barcodes
- Detect species (Human / Mouse / Unknown)
- Build deterministic gene index
- Optional metadata join

## Stage 2: Expression Access and Normalization
- Read MTX numeric values in CSC order
- Or read CSC from shared `kira-organelle.bin` backend
- Validate dimensions vs features/barcodes
- Compute per-cell `libsize` and `nnz`
- Optional log-normalization
- Deterministic cache for normalized values

## Stage 3: Panels
- Built-in panel definitions
- Species-aware mapping of panel genes
- Per-cell panel sums, detected counts, and coverage

## Stage 4: Axes
- TBI, RCI, PDS, TRS, NSAI
- Deterministic entropy and ratio calculations
- Nuclear Genome Stability per-cell metrics:
  - core panel trimmed means: `replication_core`, `ddr_core`, `hr_core`, `nhej_core`, `sphase_core`, `senescence_core`
  - robust within-sample z-scores (median/MAD)
  - derived scores: `RSS`, `DDR`, `RB`, `CDS`, `SAS`
  - flags: `replication_stress_high`, `checkpoint_addicted`, `senescent_like`, `genomic_instability_risk`

## Stage 5: Composite Scores
- NPS, CI, RLS
- Confidence model
- Driver attribution per score

## Stage 6: Regimes and Flags
- Deterministic regime classification
- Quality/confounder flags

## Stage 7: Reporting
- Standalone mode outputs:
  - `nuclearqc.tsv` (cell or sample mode)
  - `summary.json`
    - additive `genome_stability` block with panel version, thresholds, global and cluster summaries
  - `report.txt`
  - `panels_report.tsv`
- Pipeline mode outputs in `<out>/kira-nuclearqc/`:
  - `nuclearqc.tsv`
  - `summary.json`
  - `report.txt`
  - `panels_report.tsv`
  - `pipeline_step.json`

## CLI Run Mode
- `--run-mode standalone|pipeline` (default `standalone`)
- `pipeline` mode uses shared cache when present, otherwise logs WARN and falls back to 10x inputs

## Determinism Guarantees
- Stable ordering and formatting
- No parallelism
- Compile-time SIMD selection only
