#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PanelGroup {
    Housekeeping,
    Tf,
    Chromatin,
    Stress,
    Developmental,
    Proliferation,
    Program,
    Confounder,
}

#[derive(Debug, Clone, Copy)]
pub struct PanelDef {
    pub id: &'static str,
    pub name: &'static str,
    pub group: PanelGroup,
    pub genes: &'static [&'static str],
}

const HOUSEKEEPING_CORE: &[&str] = &["ACTB", "GAPDH", "RPLP0", "B2M"];
const TF_BASIC: &[&str] = &["POU5F1", "SOX2", "NANOG", "MYC"];
const STRESS_RESPONSE: &[&str] = &["FOS", "JUN", "ATF3", "HSP90AA1"];
const CHROMATIN_CORE: &[&str] = &["SMARCA4", "SMARCB1", "EZH2", "ARID1A"];
const PROLIFERATION_CORE: &[&str] = &["MKI67", "TOP2A", "PCNA", "MCM2"];
const DEVELOPMENTAL_CORE: &[&str] = &["SOX9", "PAX6", "GATA3", "TBX5"];
const IMMUNE_ACTIVATION: &[&str] = &["CD69", "CD83", "HLA-DRA", "HLA-DRB1", "CD74"];
const DIFFERENTIATION_FLUX: &[&str] = &["BCL6", "IRF4", "MYC"];
const CLONAL_ENGAGEMENT: &[&str] = &["HNRNPA1", "SRSF1", "HNRNPC", "RPLP0", "RPL13A"];
const REPLICATION_STRESS_GENES: &[&str] = &[
    "ATR", "CHEK1", "CHEK2", "RPA1", "RPA2", "RPA3", "RAD17", "CLSPN", "TIMELESS", "TIPIN",
];
const DNA_REPAIR_HR: &[&str] = &[
    "BRCA1", "BRCA2", "RAD51", "RAD51B", "RAD51C", "RAD51D", "PALB2", "BARD1", "RAD52",
];
const DNA_REPAIR_NHEJ: &[&str] = &["LIG4", "XRCC4", "XRCC5", "XRCC6", "PRKDC", "NHEJ1", "PNKP"];
const CHROMATIN_COMPACTION: &[&str] = &[
    "CBX1", "CBX3", "CBX5", "SUV39H1", "SUV39H2", "SETDB1", "EHMT2",
];
const CHROMATIN_OPEN_STATE: &[&str] = &[
    "SMARCA4", "SMARCB1", "ARID1A", "ARID1B", "KDM6A", "KAT2B", "EP300",
];
const REPLICATION_FORK_STABILITY: &[&str] = &[
    "MCM2", "MCM3", "MCM4", "MCM5", "MCM6", "MCM7", "CDC45", "GINS1",
];
const CHECKPOINT_ACTIVATION: &[&str] = &["ATM", "ATR", "CHEK1", "CHEK2", "TP53", "CDKN1A"];

const BUILTIN_PANELS: &[PanelDef] = &[
    PanelDef {
        id: "housekeeping_core",
        name: "Housekeeping Core",
        group: PanelGroup::Housekeeping,
        genes: HOUSEKEEPING_CORE,
    },
    PanelDef {
        id: "tf_basic",
        name: "TF Basic",
        group: PanelGroup::Tf,
        genes: TF_BASIC,
    },
    PanelDef {
        id: "chromatin_core",
        name: "Chromatin Core",
        group: PanelGroup::Chromatin,
        genes: CHROMATIN_CORE,
    },
    PanelDef {
        id: "stress_response",
        name: "Stress Response",
        group: PanelGroup::Stress,
        genes: STRESS_RESPONSE,
    },
    PanelDef {
        id: "developmental_core",
        name: "Developmental Core",
        group: PanelGroup::Developmental,
        genes: DEVELOPMENTAL_CORE,
    },
    PanelDef {
        id: "proliferation_core",
        name: "Proliferation Core",
        group: PanelGroup::Proliferation,
        genes: PROLIFERATION_CORE,
    },
    PanelDef {
        id: "immune_activation",
        name: "Immune Activation",
        group: PanelGroup::Program,
        genes: IMMUNE_ACTIVATION,
    },
    PanelDef {
        id: "differentiation_flux",
        name: "Differentiation Flux",
        group: PanelGroup::Program,
        genes: DIFFERENTIATION_FLUX,
    },
    PanelDef {
        id: "clonal_engagement",
        name: "Clonal Engagement",
        group: PanelGroup::Program,
        genes: CLONAL_ENGAGEMENT,
    },
    PanelDef {
        id: "replication_stress_genes",
        name: "Replication Stress",
        group: PanelGroup::Confounder,
        genes: REPLICATION_STRESS_GENES,
    },
    PanelDef {
        id: "dna_repair_hr",
        name: "DNA Repair HR",
        group: PanelGroup::Confounder,
        genes: DNA_REPAIR_HR,
    },
    PanelDef {
        id: "dna_repair_nhej",
        name: "DNA Repair NHEJ",
        group: PanelGroup::Confounder,
        genes: DNA_REPAIR_NHEJ,
    },
    PanelDef {
        id: "chromatin_compaction",
        name: "Chromatin Compaction",
        group: PanelGroup::Confounder,
        genes: CHROMATIN_COMPACTION,
    },
    PanelDef {
        id: "chromatin_open_state",
        name: "Chromatin Open State",
        group: PanelGroup::Confounder,
        genes: CHROMATIN_OPEN_STATE,
    },
    PanelDef {
        id: "replication_fork_stability",
        name: "Replication Fork Stability",
        group: PanelGroup::Confounder,
        genes: REPLICATION_FORK_STABILITY,
    },
    PanelDef {
        id: "checkpoint_activation",
        name: "Checkpoint Activation",
        group: PanelGroup::Confounder,
        genes: CHECKPOINT_ACTIVATION,
    },
];

pub fn builtin_panels() -> &'static [PanelDef] {
    BUILTIN_PANELS
}
