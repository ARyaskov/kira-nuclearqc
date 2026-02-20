use super::*;
use crate::panels::defs::{PanelDef, PanelGroup};
use crate::panels::{Panel, PanelScores, PanelSet};

struct DummyAccessor {
    cols: Vec<Vec<(u32, f32)>>,
    n_genes: usize,
    libsizes: Vec<f32>,
    nnz: Vec<u32>,
}

impl ExprAccessor for DummyAccessor {
    fn n_cells(&self) -> usize {
        self.cols.len()
    }
    fn n_genes(&self) -> usize {
        self.n_genes
    }
    fn for_cell(&self, cell: usize, f: &mut dyn FnMut(u32, f32)) {
        for &(g, v) in &self.cols[cell] {
            f(g, v);
        }
    }
    fn libsize(&self, cell: usize) -> f32 {
        self.libsizes[cell]
    }
    fn nnz(&self, cell: usize) -> u32 {
        self.nnz[cell]
    }
}

fn simple_panel_set() -> PanelSet {
    let panels = vec![
        Panel {
            id: "p1",
            name: "P1",
            group: PanelGroup::Program,
            genes: vec![0, 1],
            missing: Vec::new(),
        },
        Panel {
            id: "p2",
            name: "P2",
            group: PanelGroup::Program,
            genes: vec![2],
            missing: Vec::new(),
        },
        Panel {
            id: "tf1",
            name: "TF",
            group: PanelGroup::Tf,
            genes: vec![0],
            missing: Vec::new(),
        },
        Panel {
            id: "ch1",
            name: "CH",
            group: PanelGroup::Chromatin,
            genes: vec![1],
            missing: Vec::new(),
        },
        Panel {
            id: "stress",
            name: "Stress",
            group: PanelGroup::Stress,
            genes: vec![2],
            missing: Vec::new(),
        },
        Panel {
            id: "dev",
            name: "Dev",
            group: PanelGroup::Developmental,
            genes: vec![1],
            missing: Vec::new(),
        },
    ];
    PanelSet { panels }
}

fn simple_scores() -> PanelScores {
    PanelScores {
        panel_sum: vec![
            vec![3.0, 1.0, 2.0, 1.0, 1.0, 0.5],
            vec![1.0, 1.0, 0.0, 0.0, 0.0, 0.0],
        ],
        panel_detected: vec![vec![2, 1, 1, 1, 1, 1], vec![1, 1, 0, 0, 0, 0]],
        panel_coverage: vec![
            vec![1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
            vec![0.5, 1.0, 0.0, 0.0, 0.0, 0.0],
        ],
    }
}

#[test]
fn test_entropy_uniform() {
    let values = vec![1.0f32, 1.0, 1.0, 1.0];
    let (_h, h_norm) = entropy_norm_from_values(&values);
    assert!((h_norm - 1.0).abs() < 1e-6);
}

#[test]
fn test_tbi_extremes() {
    let accessor = DummyAccessor {
        cols: vec![vec![(0, 10.0)]],
        n_genes: 10,
        libsizes: vec![10.0],
        nnz: vec![1],
    };
    let panel_set = simple_panel_set();
    let panel_scores = simple_scores();
    let mut thresholds = ThresholdProfile::default_v1();
    thresholds.expr_min = 0.0;
    thresholds.frac_rescale_min = 0.0;
    thresholds.frac_rescale_max = 1.0;

    let out = run_stage4(&accessor, &panel_set, &panel_scores, &thresholds);
    assert!(out.axes.tbi[0] >= 0.0 && out.axes.tbi[0] <= 1.0);
}

#[test]
fn test_pds_dominance() {
    let panel_set = simple_panel_set();
    let panel_scores = simple_scores();
    let accessor = DummyAccessor {
        cols: vec![vec![(0, 1.0), (1, 1.0), (2, 1.0)]],
        n_genes: 3,
        libsizes: vec![3.0],
        nnz: vec![3],
    };
    let thresholds = ThresholdProfile::default_v1();
    let out = run_stage4(&accessor, &panel_set, &panel_scores, &thresholds);
    assert!(out.axes.pds[0] > 0.0);
}

#[test]
fn test_rci_guard() {
    let panel_set = simple_panel_set();
    let mut panel_scores = simple_scores();
    panel_scores.panel_sum[0][2] = 0.0;
    panel_scores.panel_sum[0][3] = 0.0;
    let accessor = DummyAccessor {
        cols: vec![vec![(0, 1.0)]],
        n_genes: 1,
        libsizes: vec![1.0],
        nnz: vec![1],
    };
    let mut thresholds = ThresholdProfile::default_v1();
    thresholds.tf_min_sum = 5.0;
    let out = run_stage4(&accessor, &panel_set, &panel_scores, &thresholds);
    assert_eq!(out.axes.rci[0], 0.0);
    assert!(out.flags[0].low_tf_signal);
}

#[test]
fn test_determinism() {
    let panel_set = simple_panel_set();
    let panel_scores = simple_scores();
    let accessor = DummyAccessor {
        cols: vec![vec![(0, 1.0), (1, 2.0), (2, 3.0)]],
        n_genes: 3,
        libsizes: vec![6.0],
        nnz: vec![3],
    };
    let thresholds = ThresholdProfile::default_v1();
    let a = run_stage4(&accessor, &panel_set, &panel_scores, &thresholds);
    let b = run_stage4(&accessor, &panel_set, &panel_scores, &thresholds);

    assert_eq!(a.axes.tbi[0].to_bits(), b.axes.tbi[0].to_bits());
    assert_eq!(a.axes.rci[0].to_bits(), b.axes.rci[0].to_bits());
    assert_eq!(a.axes.pds[0].to_bits(), b.axes.pds[0].to_bits());
    assert_eq!(a.axes.trs[0].to_bits(), b.axes.trs[0].to_bits());
    assert_eq!(a.axes.nsai[0].to_bits(), b.axes.nsai[0].to_bits());
    assert_eq!(a.axes.rss[0].to_bits(), b.axes.rss[0].to_bits());
    assert_eq!(a.axes.drbi[0].to_bits(), b.axes.drbi[0].to_bits());
    assert_eq!(a.axes.cci[0].to_bits(), b.axes.cci[0].to_bits());
    assert_eq!(a.axes.trci[0].to_bits(), b.axes.trci[0].to_bits());
}
