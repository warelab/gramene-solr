# Universal, species-agnostic PO(int) -> harmonized organ, plus per-clade organ panels.
# Used across all Gramene expression genomes. Organs derive from Plant Ontology ids
# (with label fallbacks), so the same vocabulary applies to any species.

PO_TO_ORGAN = {
 # --- roots ---
 9005: "root", 25: "root", 256: "root", 20127: "root", 252: "root",  # incl primary root, endodermis
 # --- shoot / aerial ---
 9006: "shoot", 20033: "shoot", 1901: "shoot",                       # incl coleoptile, aerial part
 # --- leaf (incl cell types & appendages) ---
 25034: "leaf", 20103: "leaf", 20040: "leaf", 25142: "leaf", 6012: "leaf",
 20104: "leaf", 5645: "leaf", 6070: "leaf", 6023: "leaf", 20038: "leaf",  # petiole
 # --- stem ---
 9047: "stem", 20142: "stem", 3024: "stem",                          # incl stolon
 # --- meristem ---
 229: "meristem", 37: "meristem", 6079: "meristem",                  # floral/shoot apex/shoot meristem
 # --- inflorescence ---
 9049: "inflorescence", 20136: "inflorescence", 20126: "inflorescence",
 1690: "inflorescence", 9051: "inflorescence",
 # --- flower & floral organs ---
 9046: "flower", 9030: "flower", 56: "flower", 6488: "flower", 20003: "flower",
 25074: "flower", 20094: "flower", 9072: "flower", 9064: "flower", 998: "flower",
 # --- anther / pollen / stamen ---
 9066: "anther_pollen", 25281: "anther_pollen", 2: "anther_pollen",
 25121: "anther_pollen", 20099: "anther_pollen", 9029: "anther_pollen",
 # --- seed / grain ---
 9010: "seed", 6814: "seed", 9001: "seed", 7057: "seed",
 # --- fruit (eudicot) ---
 30108: "fruit", 30106: "fruit", 7038: "fruit",                      # berry, silique, ripe fruit
 # --- seed coat ---
 9088: "seed_coat", 20063: "seed_coat",
 # --- cotyledon ---
 20030: "cotyledon", 25471: "cotyledon", 25470: "cotyledon", 6057: "cotyledon", 6058: "cotyledon",
 # --- endosperm ---
 9089: "endosperm", 25196: "endosperm", 5360: "endosperm", 6220: "endosperm", 78: "endosperm",
 # --- embryo ---
 9009: "embryo", 922: "embryo", 1: "embryo", 19018: "embryo", 20108: "embryo",
 # --- tuber (potato) ---
 25522: "tuber", 25047: "tuber", 25055: "tuber", 25052: "tuber", 25081: "tuber",
 # --- pericarp (grass) ---
 9084: "pericarp",
 # --- vasculature ---
 34: "vasculature", 9015: "vasculature",
 # xylem / phloem and their cell types: vascular tissue, not separate organs. Studies that
 # sample wood (e.g. the poplar tissue panel) report these instead of "stem".
 5417: "vasculature", 5352: "vasculature", 274: "vasculature", 25417: "vasculature",
}

# exact-label fallbacks for samples with no PO int_id
LABEL_FALLBACK = {
 "root system": "root", "root (seedling)": "root",
 "nonvascular system": "shoot", "shoot (seedling)": "shoot",
 "vegetative meristem": "meristem", "shoot apex": "meristem",
 "husk": "leaf", "leaf-middle": "leaf", "ligule": "leaf", "leafbud": "leaf",
 "maternal transfer zone": "endosperm", "milk grain": "seed",
 "pod": "fruit", "pod and seed": "seed",
}

# ordered substring fallbacks (after PO + exact-label miss). First match wins, so
# specific patterns precede generic organ words.
SUBSTR_FALLBACK = [
 ("pericarp of fruit", "fruit"), ("fruit", "fruit"),
 ("seed coat", "seed_coat"),
 ("cotyledon", "cotyledon"),
 ("tuber", "tuber"),
 ("caryopsis", "seed"),
 ("axis", "embryo"), ("suspensor", "embryo"), ("embryo", "embryo"),
 ("aleurone", "endosperm"), ("pericarp", "pericarp"), ("endosperm", "endosperm"),
 ("silk", "flower"), ("nucellus", "flower"), ("ovule", "flower"),
 ("panicle", "inflorescence"), ("spikelet", "inflorescence"), ("tassel", "inflorescence"),
 ("ear", "inflorescence"),
 ("internode", "stem"),
 ("pollen", "anther_pollen"), ("anther", "anther_pollen"),
 ("leaf", "leaf"), ("root", "root"), ("stem", "stem"), ("shoot", "shoot"),
]

EXCLUDE = {5052, 9007, 1003, 5421, 5423}              # callus + generic cell types
EXCLUDE_LABELS = {"whole organism", "whole plant", "whole in vitro plant",
                  "epidermis", "parenchyma", "collenchyma"}

# Shared core + per-clade organ panels. Grass panel is unchanged (13 organs).
CLADE_PANELS = {
 "grass":   ["root", "shoot", "leaf", "stem", "meristem", "inflorescence", "flower",
             "anther_pollen", "seed", "endosperm", "embryo", "pericarp", "vasculature"],
 "eudicot": ["root", "shoot", "leaf", "stem", "meristem", "inflorescence", "flower",
             "anther_pollen", "fruit", "pericarp", "seed", "seed_coat", "cotyledon",
             "endosperm", "embryo", "tuber", "vasculature"],
 "bryophyte": [],   # moss: no baseline organ data in Gramene -> stress-only
 "algae": [],       # Chlamydomonas: unicellular -> stress-only
}
# Back-compat: the grass list as the historical ORGANS constant.
ORGANS = CLADE_PANELS["grass"]
# Union of all organs across clades (for reference distributions).
ALL_ORGANS = sorted({o for panel in CLADE_PANELS.values() for o in panel})

CLADE_OF_TAXON = {
 # grasses
 4565: "grass", 4577: "grass", 39947: "grass", 4558: "grass", 15368: "grass",
 200361: "grass", 4555: "grass", 112509: "grass", 39946: "grass",
 # eudicots
 3708: "eudicot", 3847: "eudicot", 109376: "eudicot", 3711: "eudicot", 57577: "eudicot",
 3702: "eudicot", 4113: "eudicot", 4081: "eudicot", 81972: "eudicot", 3641: "eudicot",
 3555: "eudicot", 29760: "eudicot", 3760: "eudicot", 3694: "eudicot",   # populus trichocarpa
 # non-angiosperm
 3218: "bryophyte", 3055: "algae",
}


def organ_for(po_int, label):
    if po_int in EXCLUDE:
        return None
    lab = label.lower() if label else None
    if lab and lab in EXCLUDE_LABELS:
        return None
    if po_int in PO_TO_ORGAN:
        return PO_TO_ORGAN[po_int]
    if lab and lab in LABEL_FALLBACK:
        return LABEL_FALLBACK[lab]
    if lab:
        for sub, organ in SUBSTR_FALLBACK:
            if sub in lab:
                return organ
    return None
