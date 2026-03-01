use std::collections::BTreeMap;

use super::scores::{GenomePanelAudit, GenomeStabilityCellScores, RobustNormStat};

#[derive(Debug, Clone)]
pub struct GenomeStabilityThresholds {
    pub replication_stress_high_rss: f32,
    pub ddr_high: f32,
    pub checkpoint_addicted_cds: f32,
    pub checkpoint_addicted_sas_max: f32,
    pub senescent_like_sas: f32,
    pub senescent_like_sphase_z_max: f32,
    pub gir_rb_nhej_skew_max: f32,
    pub panel_min_genes: usize,
    pub trimmed_mean_fraction: f32,
    pub robust_z_eps: f32,
}

#[derive(Debug, Clone)]
pub struct GenomeMetricDistribution {
    pub name: &'static str,
    pub median: f32,
    pub p10: f32,
    pub p90: f32,
}

#[derive(Debug, Clone)]
pub struct GenomeFlagFraction {
    pub name: &'static str,
    pub fraction: f32,
}

#[derive(Debug, Clone)]
pub struct ClusterGenomeStats {
    pub cluster_id: String,
    pub n_cells: usize,
    pub metrics: Vec<GenomeMetricDistribution>,
    pub flag_fractions: Vec<GenomeFlagFraction>,
}

#[derive(Debug, Clone)]
pub struct TopClusterStat {
    pub cluster_id: String,
    pub value: f32,
}

#[derive(Debug, Clone)]
pub struct GenomeGlobalStats {
    pub robust_norm: Vec<RobustNormStat>,
    pub distributions: Vec<GenomeMetricDistribution>,
    pub flag_fractions: Vec<GenomeFlagFraction>,
    pub missing_fraction: Vec<GenomeFlagFraction>,
    pub top_clusters_by_rss: Vec<TopClusterStat>,
    pub top_clusters_by_checkpoint_addicted: Vec<TopClusterStat>,
}

#[derive(Debug, Clone)]
pub struct GenomeStabilitySummary {
    pub panel_version: &'static str,
    pub thresholds: GenomeStabilityThresholds,
    pub panel_audits: Vec<GenomePanelAudit>,
    pub global_stats: GenomeGlobalStats,
    pub cluster_stats: Vec<ClusterGenomeStats>,
}

pub fn summarize_genome_stability(
    panel_version: &'static str,
    panel_audits: &[GenomePanelAudit],
    cells: &GenomeStabilityCellScores,
    robust_norm: &[RobustNormStat],
    cluster_labels: Option<&[String]>,
) -> GenomeStabilitySummary {
    let cluster_stats = compute_cluster_stats(cells, cluster_labels);

    let top_clusters_by_rss = top_clusters(&cluster_stats, |cluster| {
        metric_value(&cluster.metrics, "rss", |m| m.median)
    });
    let top_clusters_by_checkpoint_addicted = top_clusters(&cluster_stats, |cluster| {
        flag_value(&cluster.flag_fractions, "checkpoint_addicted")
    });

    GenomeStabilitySummary {
        panel_version,
        thresholds: GenomeStabilityThresholds {
            replication_stress_high_rss: 2.0,
            ddr_high: 2.0,
            checkpoint_addicted_cds: 2.0,
            checkpoint_addicted_sas_max: 1.5,
            senescent_like_sas: 2.0,
            senescent_like_sphase_z_max: -0.5,
            gir_rb_nhej_skew_max: -1.0,
            panel_min_genes: 3,
            trimmed_mean_fraction: 0.10,
            robust_z_eps: 1e-9,
        },
        panel_audits: panel_audits.to_vec(),
        global_stats: GenomeGlobalStats {
            robust_norm: robust_norm.to_vec(),
            distributions: vec![
                distribution("rss", &cells.rss),
                distribution("ddr", &cells.ddr),
                distribution("rb", &cells.rb),
                distribution("cds", &cells.cds),
                distribution("sas", &cells.sas),
            ],
            flag_fractions: vec![
                GenomeFlagFraction {
                    name: "replication_stress_high",
                    fraction: bool_fraction(&cells.replication_stress_high),
                },
                GenomeFlagFraction {
                    name: "ddr_high",
                    fraction: bool_fraction(&cells.ddr_high),
                },
                GenomeFlagFraction {
                    name: "checkpoint_addicted",
                    fraction: bool_fraction(&cells.checkpoint_addicted),
                },
                GenomeFlagFraction {
                    name: "senescent_like",
                    fraction: bool_fraction(&cells.senescent_like),
                },
                GenomeFlagFraction {
                    name: "genomic_instability_risk",
                    fraction: bool_fraction(&cells.genomic_instability_risk),
                },
            ],
            missing_fraction: vec![
                missing_fraction("replication_core_missing", &cells.replication_core),
                missing_fraction("ddr_core_missing", &cells.ddr_core),
                missing_fraction("hr_core_missing", &cells.hr_core),
                missing_fraction("nhej_core_missing", &cells.nhej_core),
                missing_fraction("sphase_core_missing", &cells.sphase_core),
                missing_fraction("senescence_core_missing", &cells.senescence_core),
                missing_fraction("rss_missing", &cells.rss),
                missing_fraction("ddr_missing", &cells.ddr),
                missing_fraction("rb_missing", &cells.rb),
                missing_fraction("cds_missing", &cells.cds),
                missing_fraction("sas_missing", &cells.sas),
            ],
            top_clusters_by_rss,
            top_clusters_by_checkpoint_addicted,
        },
        cluster_stats,
    }
}

fn compute_cluster_stats(
    cells: &GenomeStabilityCellScores,
    cluster_labels: Option<&[String]>,
) -> Vec<ClusterGenomeStats> {
    let Some(labels) = cluster_labels else {
        return Vec::new();
    };
    if labels.is_empty() {
        return Vec::new();
    }

    let mut by_cluster: BTreeMap<String, Vec<usize>> = BTreeMap::new();
    for (idx, raw_label) in labels.iter().enumerate() {
        let label = raw_label.trim();
        if label.is_empty() {
            continue;
        }
        by_cluster.entry(label.to_string()).or_default().push(idx);
    }

    let mut out = Vec::with_capacity(by_cluster.len());
    for (cluster_id, idxs) in by_cluster {
        let replication_stress_high = idxs
            .iter()
            .map(|&idx| cells.replication_stress_high[idx])
            .collect::<Vec<_>>();
        let ddr_high = idxs
            .iter()
            .map(|&idx| cells.ddr_high[idx])
            .collect::<Vec<_>>();
        let checkpoint_addicted = idxs
            .iter()
            .map(|&idx| cells.checkpoint_addicted[idx])
            .collect::<Vec<_>>();
        let senescent_like = idxs
            .iter()
            .map(|&idx| cells.senescent_like[idx])
            .collect::<Vec<_>>();
        let genomic_instability_risk = idxs
            .iter()
            .map(|&idx| cells.genomic_instability_risk[idx])
            .collect::<Vec<_>>();

        out.push(ClusterGenomeStats {
            cluster_id,
            n_cells: idxs.len(),
            metrics: vec![
                distribution_idx("rss", &cells.rss, &idxs),
                distribution_idx("ddr", &cells.ddr, &idxs),
                distribution_idx("rb", &cells.rb, &idxs),
                distribution_idx("cds", &cells.cds, &idxs),
                distribution_idx("sas", &cells.sas, &idxs),
            ],
            flag_fractions: vec![
                GenomeFlagFraction {
                    name: "replication_stress_high",
                    fraction: bool_fraction(&replication_stress_high),
                },
                GenomeFlagFraction {
                    name: "ddr_high",
                    fraction: bool_fraction(&ddr_high),
                },
                GenomeFlagFraction {
                    name: "checkpoint_addicted",
                    fraction: bool_fraction(&checkpoint_addicted),
                },
                GenomeFlagFraction {
                    name: "senescent_like",
                    fraction: bool_fraction(&senescent_like),
                },
                GenomeFlagFraction {
                    name: "genomic_instability_risk",
                    fraction: bool_fraction(&genomic_instability_risk),
                },
            ],
        });
    }

    out
}

fn top_clusters(
    clusters: &[ClusterGenomeStats],
    value_fn: impl Fn(&ClusterGenomeStats) -> f32,
) -> Vec<TopClusterStat> {
    let mut ranked = clusters
        .iter()
        .map(|cluster| TopClusterStat {
            cluster_id: cluster.cluster_id.clone(),
            value: value_fn(cluster),
        })
        .collect::<Vec<_>>();

    ranked.sort_by(|a, b| {
        let av = if a.value.is_finite() {
            a.value
        } else {
            f32::NEG_INFINITY
        };
        let bv = if b.value.is_finite() {
            b.value
        } else {
            f32::NEG_INFINITY
        };
        match bv.partial_cmp(&av).unwrap_or(std::cmp::Ordering::Equal) {
            std::cmp::Ordering::Equal => a.cluster_id.cmp(&b.cluster_id),
            other => other,
        }
    });
    ranked.into_iter().take(3).collect()
}

fn distribution(name: &'static str, values: &[f32]) -> GenomeMetricDistribution {
    let finite = values
        .iter()
        .copied()
        .filter(|v| v.is_finite())
        .collect::<Vec<_>>();
    quantile_distribution(name, &finite)
}

fn distribution_idx(
    name: &'static str,
    values: &[f32],
    idxs: &[usize],
) -> GenomeMetricDistribution {
    let mut filtered = Vec::with_capacity(idxs.len());
    for &idx in idxs {
        let value = values[idx];
        if value.is_finite() {
            filtered.push(value);
        }
    }
    quantile_distribution(name, &filtered)
}

fn quantile_distribution(name: &'static str, values: &[f32]) -> GenomeMetricDistribution {
    if values.is_empty() {
        return GenomeMetricDistribution {
            name,
            median: 0.0,
            p10: 0.0,
            p90: 0.0,
        };
    }
    let mut sorted = values.to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

    GenomeMetricDistribution {
        name,
        median: quantile_sorted(&sorted, 0.5),
        p10: quantile_sorted(&sorted, 0.10),
        p90: quantile_sorted(&sorted, 0.90),
    }
}

fn quantile_sorted(sorted: &[f32], q: f32) -> f32 {
    if sorted.is_empty() {
        return 0.0;
    }
    let idx = ((sorted.len() - 1) as f32 * q).ceil() as usize;
    sorted[idx]
}

fn bool_fraction(values: &[bool]) -> f32 {
    if values.is_empty() {
        return 0.0;
    }
    let mut count = 0usize;
    for &value in values {
        if value {
            count += 1;
        }
    }
    count as f32 / values.len() as f32
}

fn missing_fraction(name: &'static str, values: &[f32]) -> GenomeFlagFraction {
    if values.is_empty() {
        return GenomeFlagFraction {
            name,
            fraction: 0.0,
        };
    }
    let mut missing = 0usize;
    for value in values {
        if !value.is_finite() {
            missing += 1;
        }
    }
    GenomeFlagFraction {
        name,
        fraction: missing as f32 / values.len() as f32,
    }
}

fn metric_value(
    metrics: &[GenomeMetricDistribution],
    name: &str,
    select: impl Fn(&GenomeMetricDistribution) -> f32,
) -> f32 {
    for metric in metrics {
        if metric.name == name {
            return select(metric);
        }
    }
    f32::NAN
}

fn flag_value(flags: &[GenomeFlagFraction], name: &str) -> f32 {
    for flag in flags {
        if flag.name == name {
            return flag.fraction;
        }
    }
    f32::NAN
}
