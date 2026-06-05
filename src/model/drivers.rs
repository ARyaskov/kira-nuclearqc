#[derive(Debug, Clone)]
pub struct ScoreDrivers {
    pub nps: Vec<Vec<(&'static str, f32)>>,
    pub ci: Vec<Vec<(&'static str, f32)>>,
    pub rls: Vec<Vec<(&'static str, f32)>>,
}
