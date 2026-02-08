#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NuclearRegime {
    PlasticAdaptive,
    StressAdaptive,
    CommittedState,
    RigidDegenerative,
    TranscriptionallyCollapsed,
    TransientAdaptive,
    Unclassified,
}
