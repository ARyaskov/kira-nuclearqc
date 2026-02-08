#[derive(Debug, Clone)]
pub struct ScoreDrivers {
    pub nps: Vec<Vec<(String, f32)>>,
    pub ci: Vec<Vec<(String, f32)>>,
    pub rls: Vec<Vec<(String, f32)>>,
}
