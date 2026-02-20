#[derive(Debug, Clone)]
pub struct Axes {
    pub tbi: Vec<f32>,
    pub rci: Vec<f32>,
    pub pds: Vec<f32>,
    pub trs: Vec<f32>,
    pub nsai: Vec<f32>,
    pub iaa: Vec<f32>,
    pub dfa: Vec<f32>,
    pub cea: Vec<f32>,
    pub rss: Vec<f32>,
    pub drbi: Vec<f32>,
    pub cci: Vec<f32>,
    pub trci: Vec<f32>,
}

#[derive(Debug, Clone, Default)]
pub struct AxisDrivers {
    pub expressed_genes: u32,
    pub gene_entropy: f32,
    pub panel_entropy: f32,
    pub max_program_share: f32,
    pub tf_entropy: f32,
    pub stress_ratio: f32,
    pub dev_ratio: f32,
    pub iaa_raw: f32,
    pub dfa_raw: f32,
    pub cea_raw: f32,
    pub axis_variance: f32,
}

#[derive(Debug, Clone, Default)]
pub struct AxisFlags {
    pub low_tf_signal: bool,
}

pub fn clip01(x: f32) -> f32 {
    if x < 0.0 {
        0.0
    } else if x > 1.0 {
        1.0
    } else {
        x
    }
}
