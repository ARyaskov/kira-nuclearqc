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
mod tests {
    use super::*;
    use std::fs::{self, File};
    use std::io::{BufWriter, Write};
    use std::path::{Path, PathBuf};
    use std::sync::atomic::{AtomicUsize, Ordering};

    use crate::input::load_input;
    use crate::pipeline::stage2_normalize::{Stage2Params, build_expr_accessor};

    static DIR_COUNTER: AtomicUsize = AtomicUsize::new(0);

    fn make_temp_dir() -> PathBuf {
        let mut dir = std::env::temp_dir();
        let id = DIR_COUNTER.fetch_add(1, Ordering::SeqCst);
        dir.push(format!(
            "kira_nuclearqc_stage3_{}_{}",
            std::process::id(),
            id
        ));
        fs::create_dir_all(&dir).unwrap();
        dir
    }

    fn write_file(path: &Path, contents: &str) {
        let mut f = BufWriter::new(File::create(path).unwrap());
        f.write_all(contents.as_bytes()).unwrap();
    }

    fn write_mtx(path: &Path, rows: usize, cols: usize, entries: &[(usize, usize, i64)]) {
        let mut out = String::new();
        out.push_str("%%MatrixMarket matrix coordinate integer general\n");
        out.push_str("% generated\n");
        out.push_str(&format!("{} {} {}\n", rows, cols, entries.len()));
        for (r, c, v) in entries {
            out.push_str(&format!("{} {} {}\n", r, c, v));
        }
        write_file(path, &out);
    }

    fn setup_bundle(
        dir: &Path,
        rows: usize,
        cols: usize,
        entries: &[(usize, usize, i64)],
    ) -> InputBundle {
        let mtx_path = dir.join("matrix.mtx");
        let features_path = dir.join("features.tsv");
        let barcodes_path = dir.join("barcodes.tsv");

        let mut feats = String::new();
        let symbols = ["ACTB", "GAPDH", "SOX2", "FOS", "MKI67"];
        for i in 0..rows {
            feats.push_str(&format!("G{}\t{}\tGene Expression\n", i + 1, symbols[i]));
        }
        write_file(&features_path, &feats);

        let mut bcs = String::new();
        for i in 0..cols {
            bcs.push_str(&format!("CELL-{}\n", i + 1));
        }
        write_file(&barcodes_path, &bcs);

        write_mtx(&mtx_path, rows, cols, entries);

        load_input(dir, None).unwrap()
    }

    #[test]
    fn test_panel_scoring_simple() {
        let dir = make_temp_dir();
        let bundle = setup_bundle(&dir, 5, 2, &[(1, 1, 2), (2, 1, 1), (3, 2, 3), (4, 2, 4)]);

        let accessor = build_expr_accessor(
            &bundle,
            &Stage2Params {
                normalize: false,
                cache_normalized: false,
                cache_path: None,
            },
        )
        .unwrap();

        let output = run_stage3(&bundle, accessor.as_ref()).unwrap();
        let panels = &output.panels.panels;

        let hk_idx = panels
            .iter()
            .position(|p| p.id == "housekeeping_core")
            .unwrap();
        let tf_idx = panels.iter().position(|p| p.id == "tf_basic").unwrap();
        let stress_idx = panels
            .iter()
            .position(|p| p.id == "stress_response")
            .unwrap();

        assert_eq!(output.scores.panel_sum[0][hk_idx], 3.0);
        assert_eq!(output.scores.panel_sum[1][tf_idx], 3.0);
        assert_eq!(output.scores.panel_sum[1][stress_idx], 4.0);
    }

    #[test]
    fn test_determinism() {
        let dir = make_temp_dir();
        let bundle = setup_bundle(&dir, 5, 2, &[(1, 1, 2), (2, 1, 1), (3, 2, 3), (4, 2, 4)]);

        let accessor = build_expr_accessor(
            &bundle,
            &Stage2Params {
                normalize: false,
                cache_normalized: false,
                cache_path: None,
            },
        )
        .unwrap();

        let a = run_stage3(&bundle, accessor.as_ref()).unwrap();
        let b = run_stage3(&bundle, accessor.as_ref()).unwrap();

        assert_eq!(a.scores.panel_sum, b.scores.panel_sum);
        assert_eq!(a.scores.panel_detected, b.scores.panel_detected);
        assert_eq!(a.scores.panel_coverage, b.scores.panel_coverage);
    }
}
