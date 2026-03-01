#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum GenomePanelId {
    ReplicationStress,
    Ddr,
    Hr,
    Nhej,
    SPhase,
    Senescence,
}

impl GenomePanelId {
    pub const fn as_str(self) -> &'static str {
        match self {
            GenomePanelId::ReplicationStress => "replication_stress",
            GenomePanelId::Ddr => "ddr",
            GenomePanelId::Hr => "hr",
            GenomePanelId::Nhej => "nhej",
            GenomePanelId::SPhase => "sphase",
            GenomePanelId::Senescence => "senescence",
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub struct PanelDef {
    pub id: GenomePanelId,
    pub genes: &'static [&'static str],
}

pub const GENOME_STABILITY_PANEL_VERSION: &str = "GENOME_STABILITY_PANEL_V1";
pub const PANEL_MIN_GENES: usize = 3;
pub const TRIM_FRACTION: f32 = 0.10;

pub const REPLICATION_STRESS_GENES: &[&str] = &[
    "ATR", "CHEK1", "CHEK2", "TOPBP1", "CLSPN", "TIMELESS", "TIPIN", "RAD17", "HUS1", "RAD9A",
    "RAD1", "RPA1", "RPA2", "RPA3", "MCM2", "MCM3", "MCM4", "MCM5", "MCM6", "MCM7", "CDC45",
    "GINS1", "GINS2", "GINS3", "GINS4", "PCNA",
];

pub const DDR_GENES: &[&str] = &["TP53", "CDKN1A", "GADD45A", "MDM2", "BAX", "BBC3"];

pub const HR_GENES: &[&str] = &[
    "BRCA1", "BRCA2", "RAD51", "PALB2", "BRIP1", "BARD1", "RAD50", "MRE11", "NBN",
];

pub const NHEJ_GENES: &[&str] = &["XRCC5", "XRCC6", "PRKDC", "LIG4", "XRCC4", "DCLRE1C"];

pub const SPHASE_GENES: &[&str] = &[
    "MKI67", "TOP2A", "CCNB1", "CCNB2", "CDC20", "CDK1", "PCNA", "TYMS", "RRM2",
];

pub const SENESCENCE_GENES: &[&str] = &["CDKN2A", "CDKN1A", "GLB1", "SERPINE1", "IGFBP7"];

pub const PANELS: &[PanelDef] = &[
    PanelDef {
        id: GenomePanelId::ReplicationStress,
        genes: REPLICATION_STRESS_GENES,
    },
    PanelDef {
        id: GenomePanelId::Ddr,
        genes: DDR_GENES,
    },
    PanelDef {
        id: GenomePanelId::Hr,
        genes: HR_GENES,
    },
    PanelDef {
        id: GenomePanelId::Nhej,
        genes: NHEJ_GENES,
    },
    PanelDef {
        id: GenomePanelId::SPhase,
        genes: SPHASE_GENES,
    },
    PanelDef {
        id: GenomePanelId::Senescence,
        genes: SENESCENCE_GENES,
    },
];
