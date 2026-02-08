#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Flag {
    LowExprGenes,
    LowPanelCoverage,
    MissingKeyPanels,
    HighProgramDominance,
    HighStressBias,
    LowTfSignal,
    AmbientRnaRisk,
    CellCycleConfounder,
    LowConfidence,
    ModelLimitation,
    BiologicalSilence,
}

pub fn flag_order() -> &'static [Flag] {
    &[
        Flag::LowExprGenes,
        Flag::LowPanelCoverage,
        Flag::MissingKeyPanels,
        Flag::HighProgramDominance,
        Flag::HighStressBias,
        Flag::LowTfSignal,
        Flag::AmbientRnaRisk,
        Flag::CellCycleConfounder,
        Flag::LowConfidence,
        Flag::ModelLimitation,
        Flag::BiologicalSilence,
    ]
}
