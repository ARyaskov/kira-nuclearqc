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
];

pub fn builtin_panels() -> &'static [PanelDef] {
    BUILTIN_PANELS
}
