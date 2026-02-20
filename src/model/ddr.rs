use crate::model::axes::clip01;

#[derive(Debug, Clone, Default)]
pub struct DdrMetrics {
    pub rss: Vec<f32>,
    pub drbi: Vec<f32>,
    pub cci: Vec<f32>,
    pub trci: Vec<f32>,
}

pub fn compute_ddr_metrics(
    replication_stress_norm: &[f32],
    checkpoint_activation_norm: &[f32],
    replication_fork_stability_norm: &[f32],
    hr_norm: &[f32],
    nhej_norm: &[f32],
    chromatin_compaction_norm: &[f32],
    chromatin_open_norm: &[f32],
    transcriptional_activity_proxy: &[f32],
) -> DdrMetrics {
    let n_cells = replication_stress_norm.len();
    let mut rss = vec![0.0; n_cells];
    let mut drbi = vec![0.0; n_cells];
    let mut cci = vec![0.0; n_cells];
    let mut trci = vec![0.0; n_cells];

    for cell in 0..n_cells {
        let replication_stress = replication_stress_norm[cell];
        let checkpoint = checkpoint_activation_norm[cell];
        let fork_stability = replication_fork_stability_norm[cell];
        let fork_instability = clip01(1.0 - fork_stability);
        let hr = hr_norm[cell];
        let nhej = nhej_norm[cell];
        let compaction = chromatin_compaction_norm[cell];
        let open = chromatin_open_norm[cell];
        let transcription = transcriptional_activity_proxy[cell];

        let rss_raw = 0.40 * replication_stress + 0.30 * checkpoint + 0.20 * fork_instability
            - 0.20 * fork_stability;
        rss[cell] = clip01(rss_raw);

        let drbi_raw = hr - nhej;
        drbi[cell] = rescale_to_0_1(drbi_raw);

        let cci_raw = 0.50 * compaction - 0.40 * open;
        cci[cell] = clip01(cci_raw);

        let trci_raw = 0.35 * replication_stress + 0.35 * transcription - 0.25 * fork_stability;
        trci[cell] = clip01(trci_raw);
    }

    DdrMetrics {
        rss,
        drbi,
        cci,
        trci,
    }
}

pub fn rescale_to_0_1(x: f32) -> f32 {
    let v = (x + 1.0) * 0.5;
    clip01(v)
}

#[cfg(test)]
#[path = "../../tests/src_inline/model/ddr.rs"]
mod tests;
