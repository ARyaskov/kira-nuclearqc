use crate::model::axes::clip01;
use crate::model::ddr::{compute_ddr_metrics, rescale_to_0_1};

#[test]
fn test_ddr_metric_formulae() {
    let replication_stress = vec![0.8];
    let checkpoint = vec![0.6];
    let fork_stability = vec![0.4];
    let hr = vec![0.7];
    let nhej = vec![0.2];
    let compaction = vec![0.9];
    let open = vec![0.1];
    let transcription = vec![0.5];

    let out = compute_ddr_metrics(
        &replication_stress,
        &checkpoint,
        &fork_stability,
        &hr,
        &nhej,
        &compaction,
        &open,
        &transcription,
    );

    let fork_instability = clip01(1.0 - 0.4);
    let rss = clip01(0.40 * 0.8 + 0.30 * 0.6 + 0.20 * fork_instability - 0.20 * 0.4);
    let drbi = rescale_to_0_1(0.7 - 0.2);
    let cci = clip01(0.50 * 0.9 - 0.40 * 0.1);
    let trci = clip01(0.35 * 0.8 + 0.35 * 0.5 - 0.25 * 0.4);

    assert!((out.rss[0] - rss).abs() < 1e-6);
    assert!((out.drbi[0] - drbi).abs() < 1e-6);
    assert!((out.cci[0] - cci).abs() < 1e-6);
    assert!((out.trci[0] - trci).abs() < 1e-6);
}

#[test]
fn test_drbi_rescale_boundaries() {
    assert_eq!(rescale_to_0_1(-1.0), 0.0);
    assert_eq!(rescale_to_0_1(0.0), 0.5);
    assert_eq!(rescale_to_0_1(1.0), 1.0);
    assert_eq!(rescale_to_0_1(-2.0), 0.0);
    assert_eq!(rescale_to_0_1(2.0), 1.0);
}

#[test]
fn test_ddr_determinism_bits() {
    let replication_stress = vec![0.1, 0.2];
    let checkpoint = vec![0.2, 0.3];
    let fork_stability = vec![0.4, 0.5];
    let hr = vec![0.6, 0.7];
    let nhej = vec![0.2, 0.1];
    let compaction = vec![0.3, 0.4];
    let open = vec![0.2, 0.1];
    let transcription = vec![0.5, 0.6];

    let a = compute_ddr_metrics(
        &replication_stress,
        &checkpoint,
        &fork_stability,
        &hr,
        &nhej,
        &compaction,
        &open,
        &transcription,
    );
    let b = compute_ddr_metrics(
        &replication_stress,
        &checkpoint,
        &fork_stability,
        &hr,
        &nhej,
        &compaction,
        &open,
        &transcription,
    );

    assert_eq!(a.rss[0].to_bits(), b.rss[0].to_bits());
    assert_eq!(a.drbi[1].to_bits(), b.drbi[1].to_bits());
    assert_eq!(a.cci[0].to_bits(), b.cci[0].to_bits());
    assert_eq!(a.trci[1].to_bits(), b.trci[1].to_bits());
}
