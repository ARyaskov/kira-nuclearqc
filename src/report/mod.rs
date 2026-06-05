use crate::metrics::genome_stability::aggregate::GenomeStabilitySummary;

pub mod json;
pub mod text;

#[derive(Debug, Clone)]
pub struct NamedStats {
    pub name: &'static str,
    pub median: f32,
    pub p90: f32,
    pub p99: f32,
}

#[derive(Debug, Clone)]
pub struct RegimeStat {
    pub name: &'static str,
    pub count: usize,
    pub fraction: f32,
}

#[derive(Debug, Clone)]
pub struct SummaryData {
    pub tool_name: String,
    pub tool_version: String,
    pub git_hash: Option<String>,
    pub simd_backend: String,
    pub run_mode: String,
    pub resolution: String,

    pub n_cells: usize,
    pub n_genes_raw: usize,
    pub n_genes_mappable: usize,
    pub species: String,

    pub normalize: bool,
    pub scale: f32,
    pub log1p: bool,
    pub axis_activation_mode: String,
    pub confidence_breakdown: Option<[f32; 4]>,
    pub scoring_mode: String,

    pub confidence_median: f32,
    pub confidence_p10: f32,
    pub low_confidence_fraction: f32,
    pub low_expr_fraction: f32,

    pub axes: Vec<NamedStats>,
    pub ddr_metrics: Vec<NamedStats>,
    pub composites: Vec<NamedStats>,

    pub regimes: Vec<RegimeStat>,

    pub trs_ge_0_75: f32,
    pub nps_ge_0_60: f32,
    pub rls_le_0_35: f32,

    pub missing_genes_by_panel: Vec<(String, Vec<String>)>,
    pub rls_contributors_top: Vec<String>,
    pub genome_stability: GenomeStabilitySummary,
}

#[derive(Debug, Clone)]
pub struct ReportContext {
    pub n_cells: usize,
    pub regimes: Vec<RegimeStat>,
    pub nps_median: f32,
    pub ci_median: f32,
    pub nsai_median: f32,
    pub rls_median: f32,
    pub low_confidence_fraction: f32,
    pub low_expr_fraction: f32,
    pub ambient_rna_fraction: f32,
    pub cell_cycle_fraction: f32,
    pub immune_note: bool,
    pub confidence_breakdown: Option<[f32; 4]>,
    pub rls_contributors_top: Vec<String>,
    pub rls_tail_fraction: f32,
    pub immune_tail_note: bool,
    pub scoring_mode: String,
    pub axis_activation_mode: String,
    pub confidence_model: String,
}

pub fn format_f32_6(v: f32) -> String {
    format!("{:.6}", v)
}

/// Canonical quantile rule: input sorted ascending, index `ceil((n-1) * p)`,
/// empty input returns 0.0. All median/percentile helpers route through this.
/// Callers with possible NaN must filter first (sorting non-finite is UB).
pub fn quantile_sorted(sorted: &[f32], p: f32) -> f32 {
    if sorted.is_empty() {
        return 0.0;
    }
    let n = sorted.len();
    let idx = (((n - 1) as f32) * p).ceil() as usize;
    let idx = idx.min(n - 1);
    sorted[idx]
}

pub fn quantile_indexed(values: &[f32], p: f32) -> f32 {
    if values.is_empty() {
        return 0.0;
    }
    let mut sorted = values.to_vec();
    sort_floats(&mut sorted);
    quantile_sorted(&sorted, p)
}

#[inline]
pub fn sort_floats(values: &mut [f32]) {
    values.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
}

/// Sort once, then read out (median, p90, p99).
pub fn triple_quantiles(values: &[f32]) -> (f32, f32, f32) {
    if values.is_empty() {
        return (0.0, 0.0, 0.0);
    }
    let mut sorted = values.to_vec();
    sort_floats(&mut sorted);
    (
        quantile_sorted(&sorted, 0.5),
        quantile_sorted(&sorted, 0.90),
        quantile_sorted(&sorted, 0.99),
    )
}

/// Sort once, read median and p10 (used by panels_report).
pub fn median_p10(values: &[f32]) -> (f32, f32) {
    if values.is_empty() {
        return (0.0, 0.0);
    }
    let mut sorted = values.to_vec();
    sort_floats(&mut sorted);
    (
        quantile_sorted(&sorted, 0.5),
        quantile_sorted(&sorted, 0.10),
    )
}

pub fn median(values: &[f32]) -> f32 {
    quantile_indexed(values, 0.5)
}

pub fn p10(values: &[f32]) -> f32 {
    quantile_indexed(values, 0.10)
}

pub fn p90(values: &[f32]) -> f32 {
    quantile_indexed(values, 0.90)
}

pub fn p99(values: &[f32]) -> f32 {
    quantile_indexed(values, 0.99)
}

pub fn bool_fraction(values: &[bool]) -> f32 {
    if values.is_empty() {
        return 0.0;
    }
    let mut count = 0usize;
    for &v in values {
        if v {
            count += 1;
        }
    }
    count as f32 / values.len() as f32
}

#[cfg(test)]
#[path = "../../tests/src_inline/report/mod.rs"]
mod tests;
