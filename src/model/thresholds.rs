#[derive(Debug, Clone)]
pub struct ThresholdProfile {
    pub expr_min: f32,
    pub min_expr_genes: u32,
    pub frac_rescale_min: f32,
    pub frac_rescale_max: f32,
    pub tf_min_sum: f32,
    pub program_min_sum: f32,
    pub tbi_w1: f32,
    pub tbi_w2: f32,
    pub tbi_w3: f32,
    pub trs_a: f32,
    pub trs_b: f32,
    pub trs_c: f32,
    pub stress_boost: f32,
    pub activation_mode: AxisActivationMode,
    pub rel_p70: f32,
    pub rel_p85: f32,
    pub confidence_low: f32,
    pub scoring_mode: NuclearScoringMode,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AxisActivationMode {
    Absolute,
    Relative,
    Hybrid,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NuclearScoringMode {
    ImmuneAware,
    StrictBulk,
}

impl ThresholdProfile {
    pub fn default_v1() -> Self {
        Self {
            expr_min: 0.0,
            min_expr_genes: 10,
            frac_rescale_min: 0.05,
            frac_rescale_max: 0.60,
            tf_min_sum: 1.0,
            program_min_sum: 1.0,
            tbi_w1: 0.4,
            tbi_w2: 0.4,
            tbi_w3: 0.2,
            trs_a: 0.4,
            trs_b: 0.3,
            trs_c: 0.3,
            stress_boost: 0.0,
            activation_mode: AxisActivationMode::Absolute,
            rel_p70: 0.70,
            rel_p85: 0.85,
            confidence_low: 0.4,
            scoring_mode: NuclearScoringMode::StrictBulk,
        }
    }

    pub fn immune_v1() -> Self {
        let mut base = Self::default_v1();
        base.activation_mode = AxisActivationMode::Hybrid;
        base.min_expr_genes = 5;
        base.tf_min_sum = 0.5;
        base.program_min_sum = 0.5;
        base.scoring_mode = NuclearScoringMode::ImmuneAware;
        base
    }
}
