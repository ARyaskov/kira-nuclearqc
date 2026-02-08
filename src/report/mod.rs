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
    pub composites: Vec<NamedStats>,

    pub regimes: Vec<RegimeStat>,

    pub trs_ge_0_75: f32,
    pub nps_ge_0_60: f32,
    pub rls_le_0_35: f32,

    pub missing_genes_by_panel: Vec<(String, Vec<String>)>,
    pub rls_contributors_top: Vec<String>,
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

pub fn quantile_indexed(values: &[f32], p: f32) -> f32 {
    if values.is_empty() {
        return 0.0;
    }
    let mut sorted = values.to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let n = sorted.len();
    let idx = ((n - 1) as f32 * p).ceil() as usize;
    sorted[idx]
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
mod tests {
    use super::*;

    #[test]
    fn test_quantiles() {
        let v = vec![1.0f32, 2.0, 3.0, 4.0, 5.0];
        assert_eq!(median(&v), 3.0);
        assert_eq!(p90(&v), 5.0);
        assert_eq!(p99(&v), 5.0);
    }
}
