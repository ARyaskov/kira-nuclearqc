#[derive(Debug, Clone)]
pub struct CompositeScores {
    pub nps: Vec<f32>,
    pub ci: Vec<f32>,
    pub rls: Vec<f32>,
    pub confidence: Vec<f32>,
    pub confidence_breakdown: Vec<[f32; 4]>,
}
