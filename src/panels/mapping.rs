use std::collections::BTreeMap;

use crate::input::{GeneIndex, Species};

pub fn build_symbol_map(gene_index: &GeneIndex) -> BTreeMap<String, u32> {
    let mut map = BTreeMap::new();
    for (gene_id, symbol) in gene_index.symbols_by_gene_id.iter().enumerate() {
        map.insert(symbol.clone(), gene_id as u32);
    }
    map
}

pub fn map_symbol(
    species: Species,
    symbol: &str,
    symbol_map: &BTreeMap<String, u32>,
) -> Option<u32> {
    let sym = normalize_symbol(symbol);
    if let Some(id) = symbol_map.get(&sym) {
        return Some(*id);
    }
    match species {
        Species::Human => None,
        Species::Unknown => None,
        Species::Mouse => {
            if let Some(mapped) = mouse_mapping(&sym) {
                if let Some(id) = symbol_map.get(mapped) {
                    return Some(*id);
                }
            }
            None
        }
    }
}

fn normalize_symbol(s: &str) -> String {
    s.trim().to_ascii_uppercase()
}

fn mouse_mapping(sym: &str) -> Option<&'static str> {
    for (human, mouse) in MOUSE_MAP {
        if *human == sym {
            return Some(*mouse);
        }
    }
    None
}

const MOUSE_MAP: &[(&str, &str)] = &[
    ("HLA-A", "H2-K1"),
    ("HLA-B", "H2-D1"),
    ("HLA-C", "H2-Q7"),
    ("HLA-DRA", "H2-AA"),
    ("HLA-DRB1", "H2-AB1"),
];
