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
    ("ATR", "Atr"),
    ("ATM", "Atm"),
    ("CHEK1", "Chek1"),
    ("CHEK2", "Chek2"),
    ("RPA1", "Rpa1"),
    ("RPA2", "Rpa2"),
    ("RPA3", "Rpa3"),
    ("RAD17", "Rad17"),
    ("CLSPN", "Clspn"),
    ("TIMELESS", "Timeless"),
    ("TIPIN", "Tipin"),
    ("BRCA1", "Brca1"),
    ("BRCA2", "Brca2"),
    ("RAD51", "Rad51"),
    ("RAD51B", "Rad51b"),
    ("RAD51C", "Rad51c"),
    ("RAD51D", "Rad51d"),
    ("PALB2", "Palb2"),
    ("BARD1", "Bard1"),
    ("RAD52", "Rad52"),
    ("LIG4", "Lig4"),
    ("XRCC4", "Xrcc4"),
    ("XRCC5", "Xrcc5"),
    ("XRCC6", "Xrcc6"),
    ("PRKDC", "Prkdc"),
    ("NHEJ1", "Nhej1"),
    ("PNKP", "Pnkp"),
    ("CBX1", "Cbx1"),
    ("CBX3", "Cbx3"),
    ("CBX5", "Cbx5"),
    ("SUV39H1", "Suv39h1"),
    ("SUV39H2", "Suv39h2"),
    ("SETDB1", "Setdb1"),
    ("EHMT2", "Ehmt2"),
    ("ARID1A", "Arid1a"),
    ("ARID1B", "Arid1b"),
    ("KDM6A", "Kdm6a"),
    ("KAT2B", "Kat2b"),
    ("EP300", "Ep300"),
    ("MCM2", "Mcm2"),
    ("MCM3", "Mcm3"),
    ("MCM4", "Mcm4"),
    ("MCM5", "Mcm5"),
    ("MCM6", "Mcm6"),
    ("MCM7", "Mcm7"),
    ("CDC45", "Cdc45"),
    ("GINS1", "Gins1"),
    ("TP53", "Trp53"),
    ("CDKN1A", "Cdkn1a"),
    ("HLA-A", "H2-K1"),
    ("HLA-B", "H2-D1"),
    ("HLA-C", "H2-Q7"),
    ("HLA-DRA", "H2-AA"),
    ("HLA-DRB1", "H2-AB1"),
];
