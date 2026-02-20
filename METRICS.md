# kira-nuclearqc Metrics Specification

This document defines canonical metrics produced by `kira-nuclearqc`, including formulas, thresholds, and constants from runtime code.

Scope:
- stage-2 normalization inputs to scoring
- stage-4 axis metrics (core + DDR)
- stage-5 composite scores and confidence
- stage-6 regime classification and flag thresholds
- output metric keys in `nuclearqc.tsv` / `summary.json`

## Canonical Conventions

1. Determinism
- No stochastic operations in metric computation.
- Stable formulas and fixed constants in code.

2. Sample axis semantics
- Metrics are computed per cell (barcode) in `--mode cell`.
- In `--mode sample`, per-cell metrics are aggregated into sample-level quantiles/fractions.

3. Value domains
- Most axes/composites are clamped to `[0, 1]` via `clip01(x)`.
- `drbi` is rescaled from `[-1, 1]` to `[0, 1]` with `(x + 1)/2`.

4. Directionality
- Higher is worse (more rigid/stressed): `tbi`, `rci`, `pds`, `trs`, `nsai`, `rss`, `cci`, `trci`, `ci`.
- Higher is more adaptive/plastic: `nps`, `rls`.
- `drbi` interpretation:
  - higher -> HR-dominant repair bias
  - lower -> NHEJ-dominant repair bias

## Notation

Let:
- `clip01(x) = min(max(x, 0), 1)`
- `rescale01(x, min, max) = clip01((x - min)/(max - min))` for `max > min`, else `0`
- `H(values)` = Shannon entropy
- `H_norm(values)` = normalized entropy by `ln(k)` where `k` is number of non-zero elements

Panel primitives per cell:
- `panel_sum[p]` = sum of expression values for genes mapped to panel `p`
- `panel_coverage[p]` = detected_genes_in_panel / panel_size
- `program_sum` = sum of all `Program` panel sums
- `stress_sum` = sum of all `Stress` panel sums
- `dev_sum` = sum of all `Developmental` panel sums

Expression preprocessing:
- if `--normalize`: `x_norm = ln(1 + count/libsize * 10000)`
- else: raw counts

## Stage-4 Axis Metrics

### Relative activation signal

For raw vector `v` across cells:
- `p70 = quantile(v, rel_p70)` with index `ceil((n-1)*rel_p70)`
- `p85 = quantile(v, rel_p85)` with index `ceil((n-1)*rel_p85)`
- `v_rel = clip01((v - p70)/(p85 - p70))` if `p85 > p70`, else `0`

Defaults:
- `rel_p70 = 0.70`
- `rel_p85 = 0.85`

### `a1_tbi` (Transcriptional Balance Index)

- `frac = expressed_genes / n_genes_mappable`
- `frac_norm = rescale01(frac, frac_rescale_min, frac_rescale_max)`
- `gene_entropy_norm = H_norm(non-zero gene values in cell)`
- `panel_entropy_norm = H_norm(program panel sums)`
- `tbi = clip01(tbi_w1*frac_norm + tbi_w2*gene_entropy_norm + tbi_w3*panel_entropy_norm)`

Defaults:
- `frac_rescale_min = 0.05`
- `frac_rescale_max = 0.60`
- `tbi_w1 = 0.4`
- `tbi_w2 = 0.4`
- `tbi_w3 = 0.2`

### `a2_rci` (Regulatory Coordination Index)

Input set: TF + chromatin panel sums.

- if `sum_tf < tf_min_sum`: `rci = 0` and `low_tf_signal = true`
- `entropy_norm = H_norm(tf+chromatin panel sums)`
- `anti_dom = 1 - max_tf/sum_tf`
- `rci = 0.5*entropy_norm + 0.5*anti_dom`
- final: `clip01(rci)`

Default:
- `tf_min_sum = 1.0` (immune profile: `0.5`)

### `a3_pds` (Program Dominance Score)

Input set: program panel sums.

- if `program_sum < program_min_sum` or `program_sum == 0`: `pds = 0`
- `max_share = top1/program_sum`
- `top3_share = (top1+top2+top3)/program_sum`
- `pds = clip01(0.7*max_share + 0.3*top3_share)`

Default:
- `program_min_sum = 1.0` (immune profile: `0.5`)

### `a4_trs` (Terminal Rigidity Score)

- `trs = clip01(trs_a*(1 - tbi) + trs_b*(1 - rci) + trs_c*pds)`

Defaults:
- `trs_a = 0.4`
- `trs_b = 0.3`
- `trs_c = 0.3`

### `a5_nsai` (Nuclear Stress Adaptation Index)

- if `program_sum < program_min_sum` or `program_sum == 0`: `nsai = 0`
- `stress_ratio = stress_sum/program_sum`
- `dev_ratio = dev_sum/program_sum`
- `nsai = clip01(stress_ratio - dev_ratio + stress_boost)`

Default:
- `stress_boost = 0.0`

### `a6_iaa`, `a7_dfa`, `a8_cea` (immune-aware program axes)

Raw axis values are panel sums from:
- `immune_activation`
- `differentiation_flux`
- `clonal_engagement`

Activation modes:
- `Absolute`: `clip01(raw)`
- `Relative`: `relative_score(raw)` from p70/p85 transform above
- `Hybrid`: `clip01(0.5*clip01(raw) + 0.5*relative_score(raw))`

Profiles:
- default strict profile: `Absolute`
- immune-aware profile: `Hybrid`

## DDR Metrics

DDR uses normalized relative inputs from:
- `replication_stress_genes`
- `checkpoint_activation`
- `replication_fork_stability`
- `dna_repair_hr`
- `dna_repair_nhej`
- `chromatin_compaction`
- `chromatin_open_state`
- `transcriptional_activity_proxy = tbi`

### `rss`
- `fork_instability = clip01(1 - fork_stability)`
- `rss_raw = 0.40*replication_stress + 0.30*checkpoint + 0.20*fork_instability - 0.20*fork_stability`
- `rss = clip01(rss_raw)`

### `drbi`
- `drbi_raw = hr - nhej`
- `drbi = clip01((drbi_raw + 1)/2)`

### `cci`
- `cci_raw = 0.50*compaction - 0.40*open`
- `cci = clip01(cci_raw)`

### `trci`
- `trci_raw = 0.35*replication_stress + 0.35*transcription - 0.25*fork_stability`
- `trci = clip01(trci_raw)`

Additional derived diagnostic:
- `axis_variance` = population variance of 12 axes (`tbi,rci,pds,trs,nsai,iaa,dfa,cea,rss,drbi,cci,trci`)

## Stage-5 Composite Scores

### `c1_nps`
- `nps = clip01(0.45*tbi + 0.35*rci - 0.20*pds - 0.20*trs)`

### `c2_ci`
- base: `ci = clip01(0.55*trs + 0.45*pds - 0.15*tbi)`
- with DDR enabled (default in current pipeline): `ci = clip01(ci + 0.15*cci)`

### `c3_rls`

Immune-aware mode:
- `axis_var_norm = clip01(axis_variance / 0.05)`
- `rigid_commit = max(trs, pds)`
- `rls = clip01(0.35*tbi + 0.20*dfa + 0.20*iaa + 0.15*nsai + 0.10*axis_var_norm - 0.30*rigid_commit)`
- floor rule: if not `allow_zero` and any of `p90(iaa), p90(dfa), p90(nsai) >= 0.8`, then `rls = max(rls, 0.1)`
- `allow_zero` when `tbi<0.2 && dfa<0.2 && iaa<0.2 && nsai<0.2 && axis_var_norm<0.05 && confidence>=0.6`
- with DDR enabled: `rls = clip01(rls - 0.25*rss - 0.20*trci)`

Strict mode (`--strict-nuclear`) legacy:
- `rls_base = clip01(0.45*tbi + 0.35*rci - 0.25*pds - 0.15*nsai)`
- `rls = rls_base * confidence`
- with DDR enabled: same subtraction `-0.25*rss -0.20*trci` and clamp

## Confidence Score

### Immune-aware confidence (default)

Components:
- `panel_coverage_score = 0` if key panels missing, else `clip01(key_panel_coverage_median / 0.6)`
- `expression_support_score = clip01(sqrt(panel_nonzero_fraction))`
- `axis_structure_score = clip01(axis_variance / 0.05)`
- `consistency_score = clip01(1 - penalty)`, where:
  - `penalty = max(0, trs+tbi-1.2) + max(0, pds+tbi-1.2) + max(0, trs+rci-1.3)`

Final:
- `confidence = clip01(0.30*panel_coverage + 0.25*expr_support + 0.25*axis_structure + 0.20*consistency)`
- fallback to `0` when coverage/nonzero unavailable and `axis_structure_score == 0`
- minimum floor: if key panels are present and `axis_structure_score >= 0.2`, then `confidence = max(confidence, 0.2)`

### Strict legacy confidence

- `a = 0.5*key_panel_coverage_median`
- `b = 0.3*expr_frac`
- `c = 0.2*(1 - ambient_rna_risk)`
- `confidence = clip01(a*b*c)`

## Regime Classification

Regime order is strict and first-match wins:

1. `TranscriptionallyCollapsed`
- `expressed_genes < min_expr_genes`
- OR `(tbi < 0.15 && gene_entropy < 0.10 && program_sum < program_min_sum)`

2. `RigidDegenerative`
- `trs >= 0.75 && nsai >= 0.55 && rci <= 0.35`

3. `CommittedState`
- `trs >= 0.70 && pds >= 0.60 && tbi <= 0.45 && nsai < 0.55`

4. `StressAdaptive`
- `nsai >= 0.65 && rci >= 0.35 && (tbi >= 0.35 || pds <= 0.60)`

5. `PlasticAdaptive`
- `nps >= 0.60 && trs <= 0.45 && pds <= 0.50`

6. `TransientAdaptive` (immune-aware mode only)
- `(nps >= 0.45 || iaa >= 0.35 || dfa >= 0.35) && trs <= 0.55 && pds <= 0.65`

7. `Unclassified`

## Flag Thresholds

- `LowExprGenes`: `expressed_genes < min_expr_genes`
- `LowPanelCoverage`: `key_panel_coverage_median < 0.4`
- `MissingKeyPanels`: any panel has `panel_size_mappable == 0`
- `HighProgramDominance`: `pds > 0.75`
- `HighStressBias`: `nsai > 0.75`
- `LowTfSignal`: `sum_tf_panels < tf_min_sum`
- `AmbientRnaRisk`: input ambient flag true
- `CellCycleConfounder`: `proliferation_program_share > 0.5`
- `LowConfidence`: `confidence < confidence_low` and (`strict mode` OR `axis_variance < 0.01`)
- `HighReplicationStress`: `rss > 0.70`
- `HrDominantRepair`: `drbi > 0.75`
- `NhejDominantRepair`: `drbi < 0.25`
- `ChromatinHypercompact`: `cci > 0.75`
- `HighTrConflict`: `trci > 0.70`
- `ModelLimitation`: activation mode not absolute OR `iaa>0` OR `dfa>0` OR `cea>0`
- `BiologicalSilence`: set only when no `ModelLimitation` and `confidence >= confidence_low`

## Default Constant Profiles

`default_v1`:
- `expr_min=0.0`
- `min_expr_genes=10`
- `frac_rescale_min=0.05`
- `frac_rescale_max=0.60`
- `tf_min_sum=1.0`
- `program_min_sum=1.0`
- `tbi_w1=0.4, tbi_w2=0.4, tbi_w3=0.2`
- `trs_a=0.4, trs_b=0.3, trs_c=0.3`
- `stress_boost=0.0`
- `activation_mode=Absolute`
- `rel_p70=0.70, rel_p85=0.85`
- `confidence_low=0.4`
- `scoring_mode=StrictBulk`

`immune_v1` overrides:
- `activation_mode=Hybrid`
- `min_expr_genes=5`
- `tf_min_sum=0.5`
- `program_min_sum=0.5`
- `scoring_mode=ImmuneAware`

Current CLI default path uses `immune_v1` unless `--strict-nuclear` is passed.

## Output Metric Keys

Per-cell TSV (`nuclearqc.tsv`) metric columns:
- Axes: `a1_tbi`, `a2_rci`, `a3_pds`, `a4_trs`, `a5_nsai`, `a6_iaa`, `a7_dfa`, `a8_cea`
- DDR: `rss`, `drbi`, `cci`, `trci`
- Composites: `c1_nps`, `c2_ci`, `c3_rls`
- Confidence: `confidence`

Summary JSON (`summary.json`) key aggregates:
- composites medians: `nps_median`, `ci_median`, `rls_median`
- tails: `trs_ge_0_75`, `nps_ge_0_60`, `rls_le_0_35`
- DDR distributions: `rss`, `drbi`, `cci`, `trci` (`median`, `p90`, `p99`)
- confidence QC: `low_confidence_fraction`, `confidence_median`, `confidence_p10`
