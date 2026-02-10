pub mod defs;
pub mod loader;
pub mod mapping;

pub use defs::PanelGroup;

#[derive(Debug, Clone)]
pub struct Panel {
    pub id: &'static str,
    pub name: &'static str,
    pub group: PanelGroup,
    pub genes: Vec<u32>,
    pub missing: Vec<String>,
}

#[derive(Debug, Clone)]
pub struct PanelSet {
    pub panels: Vec<Panel>,
}

#[derive(Debug, Clone)]
pub struct PanelScores {
    pub panel_sum: Vec<Vec<f32>>,
    pub panel_detected: Vec<Vec<u32>>,
    pub panel_coverage: Vec<Vec<f32>>,
}

#[derive(Debug, Clone)]
pub struct PanelAudit {
    pub panel_id: String,
    pub panel_size_defined: usize,
    pub panel_size_mappable: usize,
    pub missing_genes: Vec<String>,
}

#[cfg(test)]
#[path = "../../tests/src_inline/panels/tests.rs"]
mod tests;
