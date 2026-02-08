use crate::input::{GeneIndex, Species};
use crate::panels::defs::{PanelDef, builtin_panels};
use crate::panels::mapping::{build_symbol_map, map_symbol};
use crate::panels::{Panel, PanelAudit, PanelSet};

pub fn load_panels(species: Species, gene_index: &GeneIndex) -> (PanelSet, Vec<PanelAudit>) {
    let defs = builtin_panels();
    let symbol_map = build_symbol_map(gene_index);

    let mut panels = Vec::with_capacity(defs.len());
    let mut audits = Vec::with_capacity(defs.len());

    for def in defs {
        let (panel, audit) = map_panel(def, species, &symbol_map);
        panels.push(panel);
        audits.push(audit);
    }

    (PanelSet { panels }, audits)
}

fn map_panel(
    def: &PanelDef,
    species: Species,
    symbol_map: &std::collections::BTreeMap<String, u32>,
) -> (Panel, PanelAudit) {
    let mut genes = Vec::new();
    let mut missing = Vec::new();

    for &symbol in def.genes {
        if let Some(gene_id) = map_symbol(species, symbol, symbol_map) {
            genes.push(gene_id);
        } else {
            missing.push(symbol.to_string());
        }
    }

    let audit = PanelAudit {
        panel_id: def.id.to_string(),
        panel_size_defined: def.genes.len(),
        panel_size_mappable: genes.len(),
        missing_genes: missing.clone(),
    };

    let panel = Panel {
        id: def.id,
        name: def.name,
        group: def.group,
        genes,
        missing,
    };

    (panel, audit)
}
