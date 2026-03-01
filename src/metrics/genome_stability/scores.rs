use std::collections::BTreeMap;

use crate::input::{GeneIndex, Species};
use crate::panels::mapping::{build_symbol_map, map_symbol};
use crate::pipeline::stage2_normalize::ExprAccessor;

use super::panels::{
    GENOME_STABILITY_PANEL_VERSION, GenomePanelId, PANEL_MIN_GENES, PANELS, PanelDef, TRIM_FRACTION,
};

const EPS: f32 = 1e-9;
const CONSISTENT_ZERO_MAD_VALUE: f32 = 0.0;

#[derive(Debug, Clone)]
pub struct GenomePanelAudit {
    pub panel_id: String,
    pub panel_size_defined: usize,
    pub panel_size_mappable: usize,
    pub missing_genes: Vec<String>,
}

#[derive(Debug, Clone)]
pub struct RobustNormStat {
    pub name: &'static str,
    pub median: f32,
    pub mad: f32,
}

#[derive(Debug, Clone)]
pub struct GenomeStabilityCellScores {
    pub replication_core: Vec<f32>,
    pub ddr_core: Vec<f32>,
    pub hr_core: Vec<f32>,
    pub nhej_core: Vec<f32>,
    pub sphase_core: Vec<f32>,
    pub senescence_core: Vec<f32>,

    pub rss: Vec<f32>,
    pub ddr: Vec<f32>,
    pub rb: Vec<f32>,
    pub cds: Vec<f32>,
    pub sas: Vec<f32>,

    pub replication_stress_high: Vec<bool>,
    pub ddr_high: Vec<bool>,
    pub checkpoint_addicted: Vec<bool>,
    pub senescent_like: Vec<bool>,
    pub genomic_instability_risk: Vec<bool>,
}

#[derive(Debug, Clone)]
pub struct GenomeStabilityComputation {
    pub panel_version: &'static str,
    pub panel_audits: Vec<GenomePanelAudit>,
    pub norm_stats: Vec<RobustNormStat>,
    pub cells: GenomeStabilityCellScores,
}

#[derive(Debug, Clone)]
struct ResolvedPanel {
    id: GenomePanelId,
    mapped_gene_ids: Vec<u32>,
    missing_genes: Vec<String>,
    panel_size_defined: usize,
}

pub fn compute_genome_stability(
    accessor: &dyn ExprAccessor,
    gene_index: &GeneIndex,
    species: Species,
) -> GenomeStabilityComputation {
    let resolved = resolve_panels(gene_index, species);
    let n_cells = accessor.n_cells();

    let mut gene_to_mask: BTreeMap<u32, u8> = BTreeMap::new();
    for (panel_idx, panel) in resolved.iter().enumerate() {
        for &gene_id in &panel.mapped_gene_ids {
            let entry = gene_to_mask.entry(gene_id).or_insert(0);
            *entry |= 1u8 << panel_idx;
        }
    }

    let mut replication_core = vec![f32::NAN; n_cells];
    let mut ddr_core = vec![f32::NAN; n_cells];
    let mut hr_core = vec![f32::NAN; n_cells];
    let mut nhej_core = vec![f32::NAN; n_cells];
    let mut sphase_core = vec![f32::NAN; n_cells];
    let mut senescence_core = vec![f32::NAN; n_cells];

    let mut buffers = vec![Vec::<f32>::new(); resolved.len()];
    for (idx, panel) in resolved.iter().enumerate() {
        buffers[idx].reserve(panel.mapped_gene_ids.len());
    }

    for cell in 0..n_cells {
        for values in &mut buffers {
            values.clear();
        }

        accessor.for_cell(cell, &mut |gene_id, value| {
            if let Some(mask) = gene_to_mask.get(&gene_id) {
                for panel_idx in 0..resolved.len() {
                    if (*mask & (1u8 << panel_idx)) != 0 {
                        buffers[panel_idx].push(value);
                    }
                }
            }
        });

        for (panel_idx, panel) in resolved.iter().enumerate() {
            let score = trimmed_mean(&mut buffers[panel_idx], PANEL_MIN_GENES);
            match panel.id {
                GenomePanelId::ReplicationStress => replication_core[cell] = score,
                GenomePanelId::Ddr => ddr_core[cell] = score,
                GenomePanelId::Hr => hr_core[cell] = score,
                GenomePanelId::Nhej => nhej_core[cell] = score,
                GenomePanelId::SPhase => sphase_core[cell] = score,
                GenomePanelId::Senescence => senescence_core[cell] = score,
            }
        }
    }

    let (z_replication_core, replication_stat) =
        robust_zscore(&replication_core, "replication_core");
    let (z_ddr_core, ddr_stat) = robust_zscore(&ddr_core, "ddr_core");
    let (z_hr_core, hr_stat) = robust_zscore(&hr_core, "hr_core");
    let (z_nhej_core, nhej_stat) = robust_zscore(&nhej_core, "nhej_core");
    let (z_sphase_core, sphase_stat) = robust_zscore(&sphase_core, "sphase_core");
    let (z_senescence_core, senescence_stat) = robust_zscore(&senescence_core, "senescence_core");

    let mut rss = vec![f32::NAN; n_cells];
    let mut ddr = vec![f32::NAN; n_cells];
    let mut rb = vec![f32::NAN; n_cells];
    let mut cds = vec![f32::NAN; n_cells];
    let mut sas = vec![f32::NAN; n_cells];

    let mut replication_stress_high = vec![false; n_cells];
    let mut ddr_high = vec![false; n_cells];
    let mut checkpoint_addicted = vec![false; n_cells];
    let mut senescent_like = vec![false; n_cells];
    let mut genomic_instability_risk = vec![false; n_cells];

    for cell in 0..n_cells {
        let zr = z_replication_core[cell];
        let zd = z_ddr_core[cell];
        let zh = z_hr_core[cell];
        let zn = z_nhej_core[cell];
        let zs = z_sphase_core[cell];
        let zsen = z_senescence_core[cell];

        if zr.is_finite() && zs.is_finite() && zd.is_finite() {
            rss[cell] = 0.45 * zr + 0.25 * zs + 0.30 * zd;
        }

        if zd.is_finite() {
            ddr[cell] = zd;
        }

        if zh.is_finite() && zn.is_finite() {
            rb[cell] = zh - zn;
        }

        if zr.is_finite() && zd.is_finite() {
            let mut score = 0.6 * zr + 0.4 * zd;
            if zsen.is_finite() {
                score -= 0.2 * zsen.max(0.0);
            }
            cds[cell] = score;
        }

        if zsen.is_finite() && zd.is_finite() {
            sas[cell] = 0.7 * zsen + 0.3 * zd;
        }

        if rss[cell].is_finite() {
            replication_stress_high[cell] = rss[cell] >= 2.0;
        }
        if ddr[cell].is_finite() {
            ddr_high[cell] = ddr[cell] >= 2.0;
        }
        if cds[cell].is_finite() && sas[cell].is_finite() {
            checkpoint_addicted[cell] = cds[cell] >= 2.0 && sas[cell] < 1.5;
        }
        if sas[cell].is_finite() && zs.is_finite() {
            senescent_like[cell] = sas[cell] >= 2.0 && zs <= -0.5;
        }
        if rss[cell].is_finite() && ddr[cell].is_finite() {
            genomic_instability_risk[cell] = (rss[cell] >= 2.0 && ddr[cell] >= 2.0)
                || (cds[cell].is_finite()
                    && rb[cell].is_finite()
                    && cds[cell] >= 2.0
                    && rb[cell] <= -1.0);
        }
    }

    GenomeStabilityComputation {
        panel_version: GENOME_STABILITY_PANEL_VERSION,
        panel_audits: resolved
            .iter()
            .map(|panel| GenomePanelAudit {
                panel_id: panel.id.as_str().to_string(),
                panel_size_defined: panel.panel_size_defined,
                panel_size_mappable: panel.mapped_gene_ids.len(),
                missing_genes: panel.missing_genes.clone(),
            })
            .collect(),
        norm_stats: vec![
            replication_stat,
            ddr_stat,
            hr_stat,
            nhej_stat,
            sphase_stat,
            senescence_stat,
        ],
        cells: GenomeStabilityCellScores {
            replication_core,
            ddr_core,
            hr_core,
            nhej_core,
            sphase_core,
            senescence_core,
            rss,
            ddr,
            rb,
            cds,
            sas,
            replication_stress_high,
            ddr_high,
            checkpoint_addicted,
            senescent_like,
            genomic_instability_risk,
        },
    }
}

fn resolve_panels(gene_index: &GeneIndex, species: Species) -> Vec<ResolvedPanel> {
    let symbol_map = build_symbol_map(gene_index);
    let mut out = Vec::with_capacity(PANELS.len());

    for panel in PANELS {
        out.push(resolve_panel(panel, species, &symbol_map));
    }

    out
}

fn resolve_panel(
    panel: &PanelDef,
    species: Species,
    symbol_map: &BTreeMap<String, u32>,
) -> ResolvedPanel {
    let mut mapped_gene_ids = Vec::with_capacity(panel.genes.len());
    let mut missing_genes = Vec::new();

    for &gene in panel.genes {
        if let Some(gene_id) = map_symbol(species, gene, symbol_map) {
            mapped_gene_ids.push(gene_id);
        } else {
            missing_genes.push(gene.to_string());
        }
    }

    ResolvedPanel {
        id: panel.id,
        mapped_gene_ids,
        missing_genes,
        panel_size_defined: panel.genes.len(),
    }
}

fn trimmed_mean(values: &mut [f32], min_genes: usize) -> f32 {
    if values.len() < min_genes {
        return f32::NAN;
    }
    values.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let n = values.len();
    let trim = (TRIM_FRACTION * n as f32).floor() as usize;
    let start = trim.min(n);
    let end = n.saturating_sub(trim);
    if start >= end {
        return mean(values);
    }
    mean(&values[start..end])
}

fn robust_zscore(values: &[f32], name: &'static str) -> (Vec<f32>, RobustNormStat) {
    let mut finite = values
        .iter()
        .copied()
        .filter(|v| v.is_finite())
        .collect::<Vec<_>>();
    if finite.is_empty() {
        return (
            vec![f32::NAN; values.len()],
            RobustNormStat {
                name,
                median: 0.0,
                mad: 0.0,
            },
        );
    }

    let median = median_in_place(&mut finite);

    let mut dev = finite
        .iter()
        .map(|v| (v - median).abs())
        .collect::<Vec<_>>();
    let mad = median_in_place(&mut dev);

    let mut out = vec![f32::NAN; values.len()];
    if mad == 0.0 {
        for (idx, value) in values.iter().enumerate() {
            if value.is_finite() {
                out[idx] = CONSISTENT_ZERO_MAD_VALUE;
            }
        }
    } else {
        let denom = 1.4826 * mad + EPS;
        for (idx, value) in values.iter().enumerate() {
            if value.is_finite() {
                out[idx] = (*value - median) / denom;
            }
        }
    }

    (out, RobustNormStat { name, median, mad })
}

fn mean(values: &[f32]) -> f32 {
    if values.is_empty() {
        return f32::NAN;
    }
    let mut sum = 0.0f64;
    for &value in values {
        sum += value as f64;
    }
    (sum / values.len() as f64) as f32
}

fn median_in_place(values: &mut [f32]) -> f32 {
    if values.is_empty() {
        return 0.0;
    }
    values.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let mid = values.len() / 2;
    if values.len() % 2 == 0 {
        (values[mid - 1] + values[mid]) * 0.5
    } else {
        values[mid]
    }
}

#[cfg(test)]
mod tests {
    use super::{median_in_place, trimmed_mean};

    #[test]
    fn trimmed_mean_applies_10_percent_trim() {
        let mut values = vec![0.0, 1.0, 2.0, 3.0, 100.0, 101.0, 102.0, 103.0, 104.0, 105.0];
        let score = trimmed_mean(&mut values, 3);
        assert!((score - 64.5).abs() < 1e-6);
    }

    #[test]
    fn trimmed_mean_respects_min_genes() {
        let mut values = vec![1.0, 2.0];
        let score = trimmed_mean(&mut values, 3);
        assert!(score.is_nan());
    }

    #[test]
    fn median_even_and_odd() {
        let mut odd = vec![3.0, 1.0, 2.0];
        let mut even = vec![4.0, 1.0, 2.0, 3.0];
        assert_eq!(median_in_place(&mut odd), 2.0);
        assert_eq!(median_in_place(&mut even), 2.5);
    }
}
