use crate::input::{InputBundle, InputError};
use crate::panels::loader::load_panels;
use crate::panels::{PanelAudit, PanelScores, PanelSet};
use crate::pipeline::stage2_normalize::ExprAccessor;

#[derive(Debug)]
pub struct Stage3Output {
    pub panels: PanelSet,
    pub scores: PanelScores,
    pub audits: Vec<PanelAudit>,
}

pub fn run_stage3(
    bundle: &InputBundle,
    accessor: &dyn ExprAccessor,
) -> Result<Stage3Output, InputError> {
    let (panel_set, audits) = load_panels(bundle.species, &bundle.gene_index);
    let scores = score_panels(accessor, &panel_set);
    Ok(Stage3Output {
        panels: panel_set,
        scores,
        audits,
    })
}

pub fn score_panels(accessor: &dyn ExprAccessor, panel_set: &PanelSet) -> PanelScores {
    let n_cells = accessor.n_cells();
    let n_panels = panel_set.panels.len();

    let mut gene_to_panels: Vec<Vec<usize>> = vec![Vec::new(); accessor.n_genes()];
    for (panel_idx, panel) in panel_set.panels.iter().enumerate() {
        for &gene_id in &panel.genes {
            let idx = gene_id as usize;
            if idx < gene_to_panels.len() {
                gene_to_panels[idx].push(panel_idx);
            }
        }
    }

    let panel_sizes: Vec<usize> = panel_set.panels.iter().map(|p| p.genes.len()).collect();

    let mut panel_sum = Vec::with_capacity(n_cells);
    let mut panel_detected = Vec::with_capacity(n_cells);
    let mut panel_coverage = Vec::with_capacity(n_cells);

    for cell in 0..n_cells {
        let mut sums = vec![0f64; n_panels];
        let mut detected = vec![0u32; n_panels];

        accessor.for_cell(cell, &mut |gene_id, value| {
            if value == 0.0 {
                return;
            }
            let panels = &gene_to_panels[gene_id as usize];
            if panels.is_empty() {
                return;
            }
            for &p in panels {
                sums[p] += value as f64;
                if value > 0.0 {
                    detected[p] += 1;
                }
            }
        });

        let mut sums_f32 = Vec::with_capacity(n_panels);
        let mut coverage = Vec::with_capacity(n_panels);
        for p in 0..n_panels {
            sums_f32.push(sums[p] as f32);
            let size = panel_sizes[p];
            if size == 0 {
                coverage.push(0.0);
            } else {
                coverage.push(detected[p] as f32 / size as f32);
            }
        }

        panel_sum.push(sums_f32);
        panel_detected.push(detected);
        panel_coverage.push(coverage);
    }

    PanelScores {
        panel_sum,
        panel_detected,
        panel_coverage,
    }
}

#[cfg(test)]
#[path = "../../tests/src_inline/pipeline/stage3_panels.rs"]
mod tests;
