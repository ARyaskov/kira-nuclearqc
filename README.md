# kira-nuclearqc

Deterministic CLI for nuclear state and transcriptional plasticity analysis from 10x scRNA-seq MTX inputs.

## Status
- Stages 1–7 implemented with deterministic outputs
- SIMD backend selection at compile time (AVX2 / NEON / scalar)
- Includes Nuclear Genome Stability metrics (replication stress, DDR, repair balance, checkpoint dependency, senescence proxies)

## Build
```bash
cargo build --release
```

## Test
```bash
cargo test -q
```

## Usage
```bash
kira-nuclearqc run --input <dir> --out <outdir> [--mode cell|sample] [--run-mode standalone|pipeline]
```

## Outputs
- `nuclearqc.tsv`
- `summary.json`
- `report.txt`
- `panels_report.tsv`

`nuclearqc.tsv` now includes additive per-cell genome-stability columns:
- cores: `replication_core`, `ddr_core`, `hr_core`, `nhej_core`, `sphase_core`, `senescence_core`
- derived: `RSS`, `DDR`, `RB`, `CDS`, `SAS`
- flags: `replication_stress_high`, `checkpoint_addicted`, `senescent_like`, `genomic_instability_risk`

`summary.json` includes additive `genome_stability` global/cluster summaries with panel coverage audits and deterministic thresholds.

### Run Modes
- `standalone` (default): reads standard 10x inputs and writes outputs directly into `--out`.
- `pipeline`: prefers shared cache `<PREFIX>.kira-organelle.bin` (or `kira-organelle.bin`) from `--input`; writes outputs into `--out/kira-nuclearqc/` and emits `pipeline_step.json`.

If `--run-mode pipeline` is used and shared cache is not found, the tool logs a warning and falls back to 10x MTX reading.

## Shared Cache
- Cache format specification: [kira-shared-sc-cache/CACHE_FILE.md](https://github.com/ARyaskov/kira-shared-sc-cache/blob/main/CACHE_FILE.md)

## Determinism
- Stable ordering for panels, regimes, and outputs
- Fixed numeric formatting
- No runtime CPU feature detection; SIMD is selected at compile time

## License

MIT

See `LICENSE`.
