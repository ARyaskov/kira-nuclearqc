
use super::*;

#[test]
fn test_parse_args_default_run_mode_standalone() {
    let args = vec![
        "run".to_string(),
        "--input".to_string(),
        "data".to_string(),
        "--out".to_string(),
        "out".to_string(),
    ];
    let parsed = parse_args(&args).unwrap();
    assert_eq!(parsed.run_mode, RunMode::Standalone);
}

#[test]
fn test_parse_args_pipeline_run_mode() {
    let args = vec![
        "run".to_string(),
        "--input".to_string(),
        "data".to_string(),
        "--out".to_string(),
        "out".to_string(),
        "--run-mode".to_string(),
        "pipeline".to_string(),
    ];
    let parsed = parse_args(&args).unwrap();
    assert_eq!(parsed.run_mode, RunMode::Pipeline);
}

#[test]
fn test_resolve_output_dir_pipeline() {
    let out = resolve_output_dir(Path::new("/tmp/out"), RunMode::Pipeline);
    assert_eq!(out, PathBuf::from("/tmp/out/kira-nuclearqc"));
}

#[test]
fn test_resolve_output_dir_standalone() {
    let out = resolve_output_dir(Path::new("/tmp/out"), RunMode::Standalone);
    assert_eq!(out, PathBuf::from("/tmp/out"));
}
