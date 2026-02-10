use super::PanelSet;
use super::defs::{PanelGroup, builtin_panels};
use super::loader::load_panels;
use super::mapping::{build_symbol_map, map_symbol};
use crate::input::{GeneIndex, Species};

fn fake_gene_index(symbols: &[&str]) -> GeneIndex {
    let symbols_by_gene_id: Vec<String> = symbols.iter().map(|s| s.to_string()).collect();
    let gene_id_by_feature = (0..symbols.len()).map(Some).collect();
    GeneIndex {
        gene_id_by_feature,
        symbols_by_gene_id,
    }
}

#[test]
fn test_panel_defs_loaded() {
    let defs = builtin_panels();
    assert!(!defs.is_empty());
    assert_eq!(defs[0].id, "housekeeping_core");
    assert_eq!(defs[0].group, PanelGroup::Housekeeping);
    assert!(
        defs[0]
            .genes
            .iter()
            .all(|g| g.chars().all(|c| !c.is_lowercase()))
    );
}

#[test]
fn test_species_mapping_human_vs_mouse() {
    let gene_index = fake_gene_index(&["ACTB", "H2-K1", "H2-D1", "H2-AB1"]);
    let map = build_symbol_map(&gene_index);

    assert_eq!(map_symbol(Species::Human, "ACTB", &map), Some(0));
    assert_eq!(map_symbol(Species::Human, "HLA-A", &map), None);

    assert_eq!(map_symbol(Species::Mouse, "HLA-A", &map), Some(1));
    assert_eq!(map_symbol(Species::Mouse, "HLA-B", &map), Some(2));
    assert_eq!(map_symbol(Species::Mouse, "HLA-DRB1", &map), Some(3));
}

#[test]
fn test_missing_genes_reported() {
    let gene_index = fake_gene_index(&["ACTB"]);
    let (panels, audits) = load_panels(Species::Human, &gene_index);
    assert!(!panels.panels.is_empty());

    let hk = audits
        .iter()
        .find(|a| a.panel_id == "housekeeping_core")
        .unwrap();
    assert!(hk.panel_size_defined > hk.panel_size_mappable);
    assert!(hk.missing_genes.len() >= 1);
}

#[test]
fn test_panel_set_order_stable() {
    let gene_index = fake_gene_index(&["ACTB", "GAPDH", "RPLP0", "B2M"]);
    let (panels_a, _) = load_panels(Species::Human, &gene_index);
    let (panels_b, _) = load_panels(Species::Human, &gene_index);

    assert_eq!(panels_a.panels.len(), panels_b.panels.len());
    for (a, b) in panels_a.panels.iter().zip(panels_b.panels.iter()) {
        assert_eq!(a.id, b.id);
    }
}
