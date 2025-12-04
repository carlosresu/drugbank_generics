"""
Unified Constants for Drug Pipeline

This file consolidates all hardcoded token lists from across the pipeline.
All scripts should import constants from here to ensure consistency.

Created: Nov 28, 2025 (Phase 1, TODO #22)
Sources: constants.py, combos_drugs.py, text_utils_drugs.py, routes_forms_drugs.py,
         generic_normalization.py, scoring.py, and debug/old_files/*.py
"""

from typing import Any, Dict, List, Optional, Set, Tuple

# ============================================================================
# TOKEN CATEGORIES
# ============================================================================

CATEGORY_GENERIC = "generic"
CATEGORY_SALT = "salt"
CATEGORY_DOSE = "dose"
CATEGORY_FORM = "form"
CATEGORY_ROUTE = "route"
CATEGORY_OTHER = "other"

# ============================================================================
# MATCHING RULES
# ============================================================================

GENERIC_MATCH_REQUIRED = True  # No match without generic match
SALT_FLEXIBLE = True           # Different salts of same drug are acceptable
DOSE_FLEXIBLE = True           # Different doses are acceptable (same ATC)

# ============================================================================
# STOPWORDS - Words to ignore during tokenization
# Merged from: NATURAL_STOPWORDS, GENERIC_JUNK_TOKENS, BASE_GENERIC_IGNORE,
#              GENERIC_BAD_SINGLE_TOKENS, COMMON_UNKNOWN_STOPWORDS
# ============================================================================

STOPWORDS: Set[str] = {
    # Natural language stopwords
    "AS", "IN", "FOR", "TO", "WITH", "EQUIV", "EQUIV.", "AND", "OF", "OR",
    "NOT", "THAN", "HAS", "DURING", "THIS", "W/", "W", "PLUS", "APPROX",
    "APPROXIMATELY", "PRE", "FILLED", "PRE-FILLED", "PER", "RATIO",
    
    # Generic junk tokens (formulation words that aren't drug names)
    "SOLUTION", "SOLUTIONS", "SOLN", "IRRIGATION", "IRRIGATING",
    "INJECTION", "INJECTIONS", "INJECTABLE", "INFUSION", "INFUSIONS",
    "DILUENT", "DILUTION", "POWDER", "POWDERS", "MICRONUTRIENT",
    "FORMULA", "FORMULATION", "WATER", "VEHICLE",
    
    # Packaging/container words
    "UNIT", "UNITS", "SYRINGE", "BAG", "BOTTLE", "BOTTLES", "VIAL", "VIALS",
    "AMPULE", "AMPOULE", "AMPUL", "SACHET", "SACHETS", "CAN", "CANS",
    "BOX", "BOXES", "DRUM", "TUBE", "TUBES", "PACK", "CARTRIDGE",
    "GLASS", "PLASTIC", "CONTAINER",
    
    # Form words (should be classified as form, not generic)
    "DROPS", "DROP", "SPRAY", "PATCH", "CAPSULE", "CAPSULES",
    "TABLET", "TABLETS", "CREAM", "OINTMENT", "GEL", "LOTION",
    "SUSPENSION", "SYRUP", "NEBULE", "SOFTGEL", "SOFTGELS",
    
    # Route words
    "ORAL", "INTRAVENOUS", "INTRAMUSCULAR", "SUBCUTANEOUS",
    "OPHTHALMIC", "NASAL", "TOPICAL", "RECTAL", "VAGINAL",
    
    # Descriptor words
    "AGENT", "AGENTS", "FORMING", "GAS", "ELEMENTAL", "EQUIVALENT",
    "PERCENT", "FAT-SOLUBLE", "WATER-SOLUBLE", "TRACE", "ELEMENTS",
    "SOLUBLE", "PARENTERAL",
    
    # Bad single tokens (ambiguous when alone)
    "ACID", "ACIDS", "VITAMIN", "VITAMINS", "ELEMENT", "MIXTURE",
    "PREPARATION", "MC", "NONE", "CONTENT",
    
    # Common noise phrases
    "EYE DROPS", "PRE FILLED", "PRE-FILLED",
    
    # Descriptors (refined - removed drug components like FOLIC, RETINOL)
    "ACETATED", "ADULT", "ANTISEPTIC", "AQUEOUS", "ATTENUATED",
    "BALANCED", "BASED", "BIPHASIC", "BLUE", "CARE", "CARPULE",
    "CELL", "CHEWABLE", "CHICK", "CLEANSER", "COATED", "COMBI",
    "COMPLEX", "CONCENTRATE", "COUNT", "DEGRADED", "DERIVATIVE",
    "DIABETES", "DIALYSATE", "DIALYSIS", "DIBASIC", "DOSE", "DOSES",
    "DRUGS", "DURULES", "EFFERVESCENT", "EMBRYO", "EQUINE", "FEVER",
    "FLUID", "FREE", "GALLON", "HEMODIALYSIS", "HEPATIC", "HEPATITIS",
    "HIGH", "HUMAN", "HYPERTONIC", "INACTIVATED", "INFANTS", "ISOPHANE",
    "JELLY", "JUNIOR", "KCAL", "LACTATED", "LIPID", "LIQUID", "LIVE",
    "LYOPHILIZED", "MAINTENANCE", "MEDICATED", "MEDICINES", "METERED",
    "MILLION", "MODIFIED", "MONOBASIC", "MONODOSE", "MULTI", "MULTIDOSE",
    "MULTIPLE", "NACL", "NEEDLE", "NORMAL", "NUTRITION", "PAINT",
    "PEDIA", "PEDIATRIC", "PETROLEUM", "PLAIN", "POUCH", "PURIFIED",
    "RECOMBINANT", "REGULAR", "REHYDRATION", "RELEASE", "RENAL",
    "REPLACEMENT", "RESPIRATORY", "ROSE", "SALINE", "SALT", "SALTS",
    "SERUM", "SINGLE", "SKIN", "SOFT", "SOLVENT", "SPINAL", "STANDARD",
    "STERILE", "SURGICAL", "WITHOUT", "YELLOW",
}

# Lowercase version for case-insensitive matching
STOPWORDS_LOWER: Set[str] = {s.lower() for s in STOPWORDS}

# ============================================================================
# FORM MODIFIER WORDS - Words that are valid drug names but should be ignored
# when they appear as form/packaging descriptors (after CAPSULE, TABLET, etc.)
# ============================================================================

FORM_MODIFIER_IGNORE: Set[str] = {
    # These are real drugs but commonly appear as form descriptors
    "GELATIN",        # DB11242 - but also describes capsule shells
    "STARCH",         # DB00930 - but also excipient
    "CELLULOSE",      # Excipient
    "LACTOSE",        # Excipient
    
    # Modifiers that appear after form words
    "COATED", "FILM", "ENTERIC", "SUGAR", "HARD", "LIQUID",
    "FILLED", "EXTENDED", "SUSTAINED", "MODIFIED", "DELAYED",
    "IMMEDIATE", "CONTROLLED", "DISPERSIBLE", "CHEWABLE",
    "EFFERVESCENT", "SUBLINGUAL", "BUCCAL", "ORALLY",
    "DISINTEGRATING", "FREEZE", "DRIED", "LYOPHILIZED",
    "DEPOT", "RETARD",
}

# ============================================================================
# SALT TOKENS - Pharmaceutical salt/hydrate suffixes
# Merged from: SALT_TOKENS (3 copies), SPECIAL_SALT_TOKENS, SALT_FORM_SUFFIXES
# ============================================================================

SALT_TOKENS: Set[str] = {
    # Cation salts
    "CALCIUM", "SODIUM", "POTASSIUM", "MAGNESIUM", "ZINC", "AMMONIUM",
    "MEGLUMINE", "ALUMINUM", "ALUMINIUM", "IRON", "FERROUS", "FERRIC",
    "COPPER", "MANGANESE", "SILVER", "LITHIUM", "BARIUM",
    
    # Anion salts / acid forms
    "HYDROCHLORIDE", "HCL", "NITRATE", "NITRITE", "SULFATE", "SULPHATE",
    "PHOSPHATE", "DIPHOSPHATE", "HYDROXIDE", "OXIDE", "DIPROPIONATE",
    "ACETATE", "DIACETATE", "TARTRATE", "BISULFATE",
    "FUMARATE", "OXALATE", "MALEATE", "MALATE", "MESYLATE", "TOSYLATE",
    "BESYLATE", "BESILATE", "BITARTRATE", "SUCCINATE", "CITRATE", "LACTATE",
    "GLUCONATE", "BICARBONATE", "CARBONATE", "BROMIDE", "CHLORIDE",
    "IODIDE", "HYDROIODIDE", "SELENITE", "THIOSULFATE", "FLUORIDE",
    
    # Ester salts
    "DECANOATE", "PALMITATE", "STEARATE", "OLEATE", "PAMOATE", "BENZOATE",
    "VALERATE", "BUTYRATE", "PROPIONATE", "HYDROBROMIDE", "DOCUSATE",
    "HEMISUCCINATE", "FUROATE", "FUSIDATE", "SALICYLATE", "FOLINATE", "ASCORBATE",
    
    # Complex salts (added from drugbank_salt_suffixes.csv)
    "ACETONIDE", "AXETIL", "BENZATHINE", "CLAVULANATE", "DIHYDROCHLORIDE",
    "DINITRATE", "DISODIUM", "GLYCERYL", "MESILATE", "PENTAHYDRATE",
    "POLYMERISATE", "SUCCINYLATED", "SULFONATE", "TRINITRATE",
    "TROMETAMOL", "TROMETHAMOL", "TROMETHAMINE",
    
    # Hydrates
    "DIHYDRATE", "TRIHYDRATE", "MONOHYDRATE", "HYDRATE", "HEMIHYDRATE",
    "ANHYDROUS",
    
    # Release modifiers (sometimes treated as salt-like)
    "SR", "XR", "ER", "CR",
}

SALT_TOKENS_LOWER: Set[str] = {s.lower() for s in SALT_TOKENS}

# ============================================================================
# PURE SALT COMPOUNDS - Should NOT have salt stripped
# These are complete drug names where the "salt" IS the active ingredient
# Merged from: PURE_SALT_COMPOUNDS, COMPOUND_GENERICS, SALT_UNIT_SET
# ============================================================================

PURE_SALT_COMPOUNDS: Set[str] = {
    # Chlorides
    "SODIUM CHLORIDE", "POTASSIUM CHLORIDE", "CALCIUM CHLORIDE",
    "MAGNESIUM CHLORIDE", "ZINC CHLORIDE", "AMMONIUM CHLORIDE",
    
    # Carbonates / Bicarbonates
    "SODIUM BICARBONATE", "POTASSIUM BICARBONATE", "CALCIUM CARBONATE",
    "MAGNESIUM CARBONATE", "SODIUM CARBONATE", "LITHIUM CARBONATE",
    
    # Sulfates
    "MAGNESIUM SULFATE", "MAGNESIUM SULPHATE", "SODIUM SULFATE",
    "POTASSIUM SULFATE", "CALCIUM SULFATE", "FERROUS SULFATE",
    "FERROUS SULPHATE", "FERRIC SULFATE", "ZINC SULFATE", "ZINC SULPHATE",
    "COPPER SULFATE", "MANGANESE SULFATE", "BARIUM SULFATE",
    
    # Phosphates
    "SODIUM PHOSPHATE", "POTASSIUM PHOSPHATE", "CALCIUM PHOSPHATE",
    "MAGNESIUM PHOSPHATE",
    
    # Citrates
    "SODIUM CITRATE", "POTASSIUM CITRATE", "CALCIUM CITRATE",
    "MAGNESIUM CITRATE", "LITHIUM CITRATE",
    
    # Lactates
    "SODIUM LACTATE", "CALCIUM LACTATE",
    
    # Acetates
    "SODIUM ACETATE", "POTASSIUM ACETATE", "CALCIUM ACETATE",
    "MAGNESIUM ACETATE",
    
    # Gluconates
    "SODIUM GLUCONATE", "CALCIUM GLUCONATE", "MAGNESIUM GLUCONATE",
    "POTASSIUM GLUCONATE", "ZINC GLUCONATE", "FERROUS GLUCONATE",
    
    # Hydroxides
    "SODIUM HYDROXIDE", "POTASSIUM HYDROXIDE", "CALCIUM HYDROXIDE",
    "MAGNESIUM HYDROXIDE", "ALUMINUM HYDROXIDE", "ALUMINIUM HYDROXIDE",
    
    # Others
    "SODIUM NITRATE", "POTASSIUM NITRATE", "SILVER NITRATE",
    "SODIUM IODIDE", "POTASSIUM IODIDE", "SODIUM BROMIDE",
    "POTASSIUM BROMIDE", "SODIUM FLUORIDE", "CALCIUM FLUORIDE",
    "POTASSIUM FLUORIDE", "SODIUM SELENITE", "SODIUM THIOSULFATE",
    "FERROUS FUMARATE", "ZINC OXIDE",
}

# ============================================================================
# COMPOUND SALT RECOGNITION - Cation/Anion mapping
# Used to identify related salts (e.g., SODIUM CHLORIDE and POTASSIUM CHLORIDE
# share the CHLORIDE anion)
# ============================================================================

# Common pharmaceutical cations (positively charged ions)
SALT_CATIONS: Set[str] = {
    "SODIUM", "POTASSIUM", "CALCIUM", "MAGNESIUM", "ZINC", "IRON",
    "FERROUS", "FERRIC", "ALUMINUM", "ALUMINIUM", "AMMONIUM",
    "COPPER", "MANGANESE", "SILVER", "LITHIUM", "BARIUM",
}

# Common pharmaceutical anions (negatively charged ions)
SALT_ANIONS: Set[str] = {
    # Halides
    "CHLORIDE", "BROMIDE", "IODIDE", "FLUORIDE",
    # Oxygen-containing
    "SULFATE", "SULPHATE", "PHOSPHATE", "NITRATE", "CARBONATE", "BICARBONATE",
    "ACETATE", "CITRATE", "LACTATE", "GLUCONATE", "FUMARATE", "SUCCINATE",
    "TARTRATE", "MALEATE", "MALATE", "OXIDE", "HYDROXIDE",
    # Other
    "SELENITE", "THIOSULFATE",
}

# Map anion to all its common cation pairs (for identifying related compounds)
ANION_TO_CATIONS: Dict[str, Set[str]] = {
    "CHLORIDE": {"SODIUM", "POTASSIUM", "CALCIUM", "MAGNESIUM", "ZINC", "AMMONIUM"},
    "SULFATE": {"MAGNESIUM", "SODIUM", "POTASSIUM", "CALCIUM", "FERROUS", "ZINC", "COPPER", "MANGANESE"},
    "SULPHATE": {"MAGNESIUM", "SODIUM", "POTASSIUM", "CALCIUM", "FERROUS", "ZINC"},
    "PHOSPHATE": {"SODIUM", "POTASSIUM", "CALCIUM", "MAGNESIUM"},
    "CARBONATE": {"CALCIUM", "MAGNESIUM", "SODIUM"},
    "BICARBONATE": {"SODIUM", "POTASSIUM"},
    "CITRATE": {"SODIUM", "POTASSIUM", "CALCIUM", "MAGNESIUM"},
    "LACTATE": {"SODIUM", "CALCIUM"},
    "ACETATE": {"SODIUM", "POTASSIUM", "CALCIUM", "MAGNESIUM"},
    "GLUCONATE": {"SODIUM", "CALCIUM", "MAGNESIUM", "POTASSIUM", "ZINC", "FERROUS"},
    "HYDROXIDE": {"SODIUM", "POTASSIUM", "CALCIUM", "MAGNESIUM", "ALUMINUM", "ALUMINIUM"},
    "NITRATE": {"SODIUM", "POTASSIUM", "SILVER"},
    "BROMIDE": {"SODIUM", "POTASSIUM"},
    "IODIDE": {"SODIUM", "POTASSIUM"},
    "FLUORIDE": {"SODIUM", "CALCIUM"},
    "FUMARATE": {"FERROUS"},
}


def parse_compound_salt(name: str) -> Tuple[Optional[str], Optional[str]]:
    """
    Parse a compound salt into its cation and anion components.
    
    Examples:
    - "SODIUM CHLORIDE" -> ("SODIUM", "CHLORIDE")
    - "CALCIUM CARBONATE" -> ("CALCIUM", "CARBONATE")
    - "FERROUS SULFATE" -> ("FERROUS", "SULFATE")
    - "PARACETAMOL" -> (None, None)
    
    Returns (cation, anion) or (None, None) if not a compound salt.
    """
    name_upper = name.upper().strip()
    words = name_upper.split()
    
    if len(words) != 2:
        return None, None
    
    first, second = words
    
    # Check if first word is cation and second is anion
    if first in SALT_CATIONS and second in SALT_ANIONS:
        return first, second
    
    return None, None


def get_related_salts(name: str) -> Set[str]:
    """
    Get all related compound salts that share the same anion.
    
    Example:
    - "SODIUM CHLORIDE" -> {"POTASSIUM CHLORIDE", "CALCIUM CHLORIDE", ...}
    """
    cation, anion = parse_compound_salt(name)
    
    if not anion:
        return set()
    
    related = set()
    related_cations = ANION_TO_CATIONS.get(anion, set())
    
    for c in related_cations:
        if c != cation:  # Don't include self
            related.add(f"{c} {anion}")
    
    return related


# ============================================================================
# ELEMENT DRUGS - Elements that can be standalone drugs
# These should be treated as generics when they appear as the main drug
# ============================================================================

ELEMENT_DRUGS: Set[str] = {
    "ZINC", "CALCIUM", "IRON", "MAGNESIUM", "POTASSIUM", "SODIUM",
    "COPPER", "MANGANESE", "SELENIUM", "CHROMIUM", "IODINE",
    "PHOSPHORUS", "FLUORIDE",
}

# ============================================================================
# FORM CANONICALIZATION
# Maps various form spellings/abbreviations to canonical forms
# ============================================================================

FORM_CANON: Dict[str, str] = {
    # Tablets - basic
    "TAB": "TABLET", "TABS": "TABLET", "TABLET": "TABLET", "TABLETS": "TABLET",
    # Tablet variations (normalize to base or keep release info)
    "TABLET, FILM COATED": "TABLET", "TABLET, COATED": "TABLET",
    "TABLET, SUGAR COATED": "TABLET", "TABLET, CHEWABLE": "TABLET",
    "TABLET, ORALLY DISINTEGRATING": "TABLET", "TABLET, EFFERVESCENT": "TABLET",
    "TABLET, SOLUBLE": "TABLET", "TABLET, MULTILAYER": "TABLET",
    "TABLET, FOR SUSPENSION": "TABLET",
    # Extended/modified release tablets (keep release info)
    "TABLET, EXTENDED RELEASE": "TABLET, EXTENDED RELEASE",
    "TABLET, FILM COATED, EXTENDED RELEASE": "TABLET, EXTENDED RELEASE",
    "TABLET, DELAYED RELEASE": "TABLET, DELAYED RELEASE",
    "TABLET, MULTILAYER, EXTENDED RELEASE": "TABLET, EXTENDED RELEASE",
    "TABLET, CHEWABLE, EXTENDED RELEASE": "TABLET, EXTENDED RELEASE",
    
    # Capsules - basic
    "CAP": "CAPSULE", "CAPS": "CAPSULE", "CAPSULE": "CAPSULE", "CAPSULES": "CAPSULE",
    "CAPLET": "CAPSULE",
    # Capsule variations
    "CAPSULE, LIQUID FILLED": "CAPSULE", "CAPSULE, GELATIN COATED": "CAPSULE",
    "CAPSULE, COATED": "CAPSULE", "CAPSULE, COATED PELLETS": "CAPSULE",
    # Extended/modified release capsules
    "CAPSULE, EXTENDED RELEASE": "CAPSULE, EXTENDED RELEASE",
    "CAPSULE, DELAYED RELEASE": "CAPSULE, DELAYED RELEASE",
    "CAPSULE, COATED, EXTENDED RELEASE": "CAPSULE, EXTENDED RELEASE",
    "CAPSULE, DELAYED RELEASE PELLETS": "CAPSULE, DELAYED RELEASE",
    
    # Solutions
    "SOLUTION": "SOLUTION", "SOLN": "SOLUTION", "SOL": "SOLUTION",
    "SOLUTIONS": "SOLUTION", "SOLUTION, CONCENTRATE": "SOLUTION",
    "ORAL SOLUTION": "SOLUTION",
    
    # Suspensions
    "SUSPENSION": "SUSPENSION", "SUSP": "SUSPENSION", "SUSPENSIONS": "SUSPENSION",
    "ORAL SUSPENSION": "SUSPENSION",
    "SUSPENSION, EXTENDED RELEASE": "SUSPENSION, EXTENDED RELEASE",
    
    # Syrups
    "SYR": "SYRUP", "SYRUP": "SYRUP", "SYRUPS": "SYRUP",
    
    # Injections
    "INJ": "INJECTION", "INJECTABLE": "INJECTION", "INJECTION": "INJECTION",
    "INJECTION, SOLUTION": "INJECTION", "INJECTION, EMULSION": "INJECTION",
    "INJECTION, SUSPENSION": "INJECTION",
    "INJECTION, POWDER, FOR SOLUTION": "INJECTION",
    "INJECTION, POWDER, FOR SUSPENSION": "INJECTION",
    "INJECTION, POWDER, LYOPHILIZED, FOR SOLUTION": "INJECTION",
    "INJECTION, POWDER, LYOPHILIZED, FOR SUSPENSION": "INJECTION",
    "INJECTION, SOLUTION, CONCENTRATE": "INJECTION",
    "INJECTABLE, LIPOSOMAL": "INJECTION",
    "INJECTION, SUSPENSION, EXTENDED RELEASE": "INJECTION, EXTENDED RELEASE",
    "INJECTION, POWDER, FOR SUSPENSION, EXTENDED RELEASE": "INJECTION, EXTENDED RELEASE",
    
    # Ampules/Vials
    "AMPULE": "AMPULE", "AMPUL": "AMPULE", "AMP": "AMPULE",
    "AMPOULE": "AMPULE", "AMPS": "AMPULE", "AMPULES": "AMPULE", "AMPULS": "AMPULE",
    "VIAL": "VIAL", "VIALS": "VIAL", "VL": "VIAL",
    
    # Bottles/Containers
    "BOT": "BOTTLE", "BOTT": "BOTTLE", "BOTTLE": "BOTTLE", "BOTTLES": "BOTTLE",
    "BAG": "BAG", "BAGS": "BAG",
    "KIT": "KIT", "KITS": "KIT",
    
    # Film coating abbreviations (formulation markers)
    "FC": "FILM COATED", "EC": "ENTERIC COATED", "SR": "EXTENDED RELEASE", 
    "XR": "EXTENDED RELEASE", "ER": "EXTENDED RELEASE", "DR": "DELAYED RELEASE",
    
    # Other formulation markers (NON-PNF means "not in Philippine National Formulary")
    "NON-PNF": "NON-FORMULARY", "NONPNF": "NON-FORMULARY",
    
    # Topical forms
    "CREAM": "CREAM", "CREAMS": "CREAM", "CREAM, AUGMENTED": "CREAM",
    "LOTION": "LOTION", "LOTIONS": "LOTION", "LOTION, AUGMENTED": "LOTION",
    "GEL": "GEL", "GELS": "GEL", "GEL, METERED": "GEL",
    "OINTMENT": "OINTMENT", "OINTMENTS": "OINTMENT",
    "PASTE": "PASTE", "PASTES": "PASTE",
    "FOAM": "FOAM", "FOAMS": "FOAM",
    "EMULSION": "EMULSION", "EMULSIONS": "EMULSION",
    "SHAMPOO": "SHAMPOO", "SHAMPOOS": "SHAMPOO",
    "SOAP": "SOAP", "SOAPS": "SOAP",
    "WASH": "WASH",
    
    # Powders
    "POWDER": "POWDER", "PWDR": "POWDER", "POWDERS": "POWDER",
    "POWDER, FOR SOLUTION": "POWDER", "POWDER, FOR SUSPENSION": "POWDER",
    "POWDER, METERED": "POWDER, METERED",  # Keep for inhalation
    
    # Granules
    "GRANULE": "GRANULE", "GRANULES": "GRANULE", "GRAN": "GRANULE",
    "GRANULE, FOR SOLUTION": "GRANULE", "GRANULE, FOR SUSPENSION": "GRANULE",
    "GRANULE, EFFERVESCENT": "GRANULE",
    "GRANULE, DELAYED RELEASE": "GRANULE, DELAYED RELEASE",
    
    # Drops
    "DROPS": "DROPS", "DROP": "DROPS", "GTTS": "DROPS",
    "EYE DROPS": "DROPS", "EAR DROPS": "DROPS", "NASAL DROPS": "DROPS",
    "ORAL DROPS": "DROPS",
    "SOLUTION / DROPS": "DROPS", "SUSPENSION / DROPS": "DROPS",
    "SOLUTION, GEL FORMING / DROPS": "DROPS",
    "INSTILL.SOLUTION": "DROPS",
    
    # Sachets
    "SACHET": "SACHET", "SACHETS": "SACHET",
    
    # Sprays
    "SPRAY": "SPRAY", "SPRAYS": "SPRAY", "NASAL SPRAY": "SPRAY",
    "DISPENSER": "SPRAY", "SPRAY, SUSPENSION": "SPRAY",
    "SPRAY, METERED": "SPRAY, METERED",  # Keep for nasal
    
    # Aerosols
    "AEROSOL": "AEROSOL", "AEROSOLS": "AEROSOL",
    "AEROSOL, SPRAY": "AEROSOL", "AEROSOL, FOAM": "AEROSOL", "AEROSOL, POWDER": "AEROSOL",
    "AEROSOL, METERED": "AEROSOL, METERED",  # Keep for inhalation
    
    # Inhalers
    "INHALER": "INHALER", "INHALERS": "INHALER", "INHALANT": "INHALER",
    "MDI": "INHALER", "DPI": "INHALER",
    "METERED DOSE INHALER": "INHALER", "DRY POWDER INHALER": "INHALER",
    "NEBULE": "NEBULE", "NEBULES": "NEBULE", "NEB": "NEBULE", "NEBS": "NEBULE",
    "NEBULIZER SOLUTION": "NEBULE",
    # WHO specific inhalation forms
    "INHAL.AEROSOL": "AEROSOL, METERED", "ORAL AEROSOL": "AEROSOL, METERED",
    "INHAL.POWDER": "POWDER, METERED",
    "INHAL.SOLUTION": "NEBULE",
    
    # Suppositories
    "SUPPOSITORY": "SUPPOSITORY", "SUPP": "SUPPOSITORY", "SUPPOSITORIES": "SUPPOSITORY",
    
    # Patches
    "PATCH": "PATCH", "PATCHES": "PATCH", "PATCH, EXTENDED RELEASE": "PATCH",
    
    # Films
    "FILM": "FILM", "FILMS": "FILM", "FILM, SOLUBLE": "FILM",
    "FILM, EXTENDED RELEASE": "FILM",
    "LAMELLA": "FILM",
    
    # Other forms
    "LOZENGE": "LOZENGE", "LOZENGES": "LOZENGE",
    "MOUTHWASH": "MOUTHWASH",
    "CHEWING GUM": "GUM", "GUM, CHEWING": "GUM", "GUM": "GUM",
    "ELIXIR": "ELIXIR",
    "IMPLANT": "IMPLANT", "IMPLANTS": "IMPLANT", "S.C. IMPLANT": "IMPLANT",
    "OVULE": "OVULE", "OVULES": "OVULE",
    "ENEMA": "ENEMA", "ENEMAS": "ENEMA",
    "BAR": "BAR", "BARS": "BAR",
    "STRIP": "STRIP", "STRIPS": "STRIP",
    "WAFER": "WAFER", "WAFERS": "WAFER",
    "PELLET": "PELLET", "PELLETS": "PELLET",
    "RING": "RING", "RINGS": "RING",
    "INSERT": "INSERT", "INSERTS": "INSERT",
    "PESSARY": "PESSARY",
    "SWAB": "SWAB", "SWABS": "SWAB",
    "CLOTH": "CLOTH", "CLOTHS": "CLOTH",
    "SPONGE": "SPONGE", "SPONGES": "SPONGE",
    "DRESSING": "DRESSING", "DRESSINGS": "DRESSING",
    "STICK": "STICK",
    "LIQUID": "LIQUID",
    "GAS": "GAS",
}

# ============================================================================
# ROUTE CANONICALIZATION
# Maps various route spellings/abbreviations to canonical routes
# ============================================================================

ROUTE_CANON: Dict[str, str] = {
    # Oral
    "PO": "ORAL", "OR": "ORAL", "O": "ORAL", "ORAL": "ORAL",
    "PER OREM": "ORAL", "BY MOUTH": "ORAL", "PER OS": "ORAL",
    
    # Intravenous
    "IV": "INTRAVENOUS", "INTRAVENOUS": "INTRAVENOUS",
    
    # Intramuscular
    "IM": "INTRAMUSCULAR", "INTRAMUSCULAR": "INTRAMUSCULAR",
    
    # Subcutaneous
    "SC": "SUBCUTANEOUS", "SUBCUTANEOUS": "SUBCUTANEOUS",
    "SUBCUT": "SUBCUTANEOUS", "SQ": "SUBCUTANEOUS",
    "SUBDERMAL": "SUBCUTANEOUS",
    
    # Parenteral (generic injectable)
    "P": "PARENTERAL", "PARENTERAL": "PARENTERAL",
    "INJ": "PARENTERAL", "INJECTION": "PARENTERAL",
    
    # Sublingual / Buccal
    "SL": "SUBLINGUAL", "SUBLINGUAL": "SUBLINGUAL",
    "BUCC": "BUCCAL", "BUCCAL": "BUCCAL",
    
    # Topical / Transdermal
    "TOPICAL": "TOPICAL", "CUTANEOUS": "TOPICAL",
    "TD": "TRANSDERMAL", "TRANSDERMAL": "TRANSDERMAL", "DERMAL": "TRANSDERMAL",
    
    # Ophthalmic
    "OPH": "OPHTHALMIC", "EYE": "OPHTHALMIC", "OPHTHALMIC": "OPHTHALMIC",
    "INTRAOCULAR": "OPHTHALMIC",
    
    # Otic
    "OTIC": "OTIC", "EAR": "OTIC",
    
    # Nasal
    "N": "NASAL", "NASAL": "NASAL", "INTRANASAL": "NASAL", "PER NASAL": "NASAL",
    
    # Inhalation
    "INH": "INHALATION", "INHAL": "INHALATION", "INHALATION": "INHALATION",
    "NEB": "INHALATION", "INHALER": "INHALATION",
    "RESPIRATORY (INHALATION)": "INHALATION",
    
    # Rectal
    "R": "RECTAL", "RECTAL": "RECTAL", "PR": "RECTAL", "PER RECTUM": "RECTAL",
    
    # Vaginal
    "V": "VAGINAL", "VAGINAL": "VAGINAL", "PV": "VAGINAL", "PER VAGINAM": "VAGINAL",
    
    # Other
    "INTRATHECAL": "INTRATHECAL",
    "INTRADERMAL": "INTRADERMAL", "ID": "INTRADERMAL",
}

# ============================================================================
# FORM TO ROUTE MAPPING
# Default route inference when only form is known
# ============================================================================

FORM_TO_ROUTE: Dict[str, str] = {
    # Oral forms (high confidence from DrugBank)
    "TABLET": "ORAL", "TAB": "ORAL", "TABS": "ORAL",
    "TABLET, EXTENDED RELEASE": "ORAL", "TABLET, DELAYED RELEASE": "ORAL",
    "CAPSULE": "ORAL", "CAP": "ORAL", "CAPS": "ORAL", "CAPLET": "ORAL",
    "CAPSULE, EXTENDED RELEASE": "ORAL", "CAPSULE, DELAYED RELEASE": "ORAL",
    "SYRUP": "ORAL", "SYRUPS": "ORAL", "SYR": "ORAL",
    "SUSPENSION": "ORAL", "SUSPENSIONS": "ORAL", "SUSP": "ORAL",
    "SUSPENSION, EXTENDED RELEASE": "ORAL",
    "SOLUTION": "ORAL", "SOLUTIONS": "ORAL", "SOLN": "ORAL",
    "SACHET": "ORAL", "GRANULE": "ORAL", "GRANULES": "ORAL",
    "GRANULE, DELAYED RELEASE": "ORAL",
    "POWDER": "ORAL",  # When not metered
    "LOZENGE": "ORAL", "MOUTHWASH": "ORAL", "ELIXIR": "ORAL",
    "GUM": "ORAL", "CHEWING GUM": "ORAL",
    "WAFER": "ORAL", "STRIP": "ORAL",
    "ORAL DROPS": "ORAL",
    
    # Topical forms (high confidence from DrugBank)
    "CREAM": "TOPICAL", "OINTMENT": "TOPICAL", "GEL": "TOPICAL",
    "LOTION": "TOPICAL", "SOAP": "TOPICAL", "SHAMPOO": "TOPICAL",
    "WASH": "TOPICAL", "PASTE": "TOPICAL", "FOAM": "TOPICAL",
    "LIQUID": "TOPICAL", "EMULSION": "TOPICAL",
    "SWAB": "TOPICAL", "CLOTH": "TOPICAL", "STICK": "TOPICAL",
    "SPONGE": "TOPICAL", "DRESSING": "TOPICAL",
    "SPRAY": "TOPICAL", "AEROSOL": "TOPICAL",
    
    # Transdermal
    "PATCH": "TRANSDERMAL",
    
    # Inhalation forms (from DrugBank/WHO)
    "INHALER": "INHALATION", "NEBULE": "INHALATION", "NEB": "INHALATION",
    "MDI": "INHALATION", "DPI": "INHALATION",
    "AEROSOL, METERED": "INHALATION", "POWDER, METERED": "INHALATION",
    "GAS": "INHALATION",
    "METERED DOSE INHALER": "INHALATION", "DRY POWDER INHALER": "INHALATION",
    "INHAL.AEROSOL": "INHALATION", "INHAL.POWDER": "INHALATION",
    "INHAL.SOLUTION": "INHALATION", "ORAL AEROSOL": "INHALATION",
    
    # Injectable forms (from DrugBank)
    "INJECTION": "INTRAVENOUS",  # Most common
    "INJECTION, EXTENDED RELEASE": "INTRAMUSCULAR",
    "AMPULE": "PARENTERAL", "AMP": "PARENTERAL", "AMPUL": "PARENTERAL",
    "AMPOULE": "PARENTERAL", "VIAL": "PARENTERAL", "VL": "PARENTERAL",
    "INJ": "PARENTERAL",
    "BOTTLE": "INTRAVENOUS", "BAG": "INTRAVENOUS",  # IV fluids
    
    # Ophthalmic
    "DROPS": "OPHTHALMIC",  # Most common for drops
    "DROP": "OPHTHALMIC", "EYE DROP": "OPHTHALMIC", "EYE DROPS": "OPHTHALMIC",
    "INSTILL.SOLUTION": "OPHTHALMIC", "LAMELLA": "OPHTHALMIC",
    
    # Otic
    "EAR DROP": "OTIC", "EAR DROPS": "OTIC",
    
    # Nasal
    "SPRAY, METERED": "NASAL", "NASAL SPRAY": "NASAL", "NASAL DROPS": "NASAL",
    
    # Rectal
    "SUPPOSITORY": "RECTAL", "SUPP": "RECTAL", "ENEMA": "RECTAL",
    
    # Vaginal
    "OVULE": "VAGINAL", "OVULES": "VAGINAL",
    "INSERT": "VAGINAL", "RING": "VAGINAL", "PESSARY": "VAGINAL",
    
    # Buccal/Sublingual
    "FILM": "BUCCAL",
    
    # Subcutaneous
    "IMPLANT": "SUBCUTANEOUS", "S.C. IMPLANT": "SUBCUTANEOUS",
}

# ============================================================================
# FORM TO ROUTES MAPPING (PLURAL) - All valid routes for each form
# Derived from DrugBank products dataset (routes with >= 1% of products)
# Routes are ordered by frequency (most common first)
# Use this for exploring all valid matching options
# ============================================================================

FORM_TO_ROUTES: Dict[str, List[str]] = {
    # Basic oral solids - unambiguous
    "TABLET": ["ORAL"],
    "TABLET, CHEWABLE": ["ORAL"],
    "TABLET, COATED": ["ORAL"],
    "TABLET, DELAYED RELEASE": ["ORAL"],
    "TABLET, EFFERVESCENT": ["ORAL"],
    "TABLET, EXTENDED RELEASE": ["ORAL"],
    "TABLET, FILM COATED": ["ORAL"],
    "TABLET, FILM COATED, EXTENDED RELEASE": ["ORAL"],
    "TABLET, FOR SUSPENSION": ["ORAL"],
    "TABLET, MULTILAYER": ["ORAL"],
    "TABLET, MULTILAYER, EXTENDED RELEASE": ["ORAL"],
    "TABLET, SUGAR COATED": ["ORAL"],
    "CAPSULE": ["ORAL"],
    "CAPSULE, COATED": ["ORAL"],
    "CAPSULE, COATED PELLETS": ["ORAL"],
    "CAPSULE, DELAYED RELEASE": ["ORAL"],
    "CAPSULE, DELAYED RELEASE PELLETS": ["ORAL"],
    "CAPSULE, EXTENDED RELEASE": ["ORAL"],
    "CAPSULE, GELATIN COATED": ["ORAL"],
    "CAPSULE, LIQUID FILLED": ["ORAL"],
    "SYRUP": ["ORAL"],
    
    # Oral forms with some alternative routes
    "TABLET, ORALLY DISINTEGRATING": ["ORAL", "SUBLINGUAL"],
    "TABLET, CHEWABLE, EXTENDED RELEASE": ["ORAL", "SUBLINGUAL"],
    "TABLET, SOLUBLE": ["ORAL", "TOPICAL"],
    "LOZENGE": ["ORAL", "TRANSMUCOSAL"],
    "TROCHE": ["SUBLINGUAL", "ORAL"],
    "GUM, CHEWING": ["ORAL", "BUCCAL"],
    "MOUTHWASH": ["ORAL", "DENTAL", "BUCCAL", "TOPICAL"],
    "RINSE": ["ORAL", "DENTAL", "TOPICAL", "BUCCAL", "CUTANEOUS"],
    "PASTILLE": ["ORAL"],
    "ELIXIR": ["ORAL", "TOPICAL"],
    
    # Solutions - highly ambiguous
    "SOLUTION": ["ORAL", "TOPICAL", "INTRAVENOUS", "SUBCUTANEOUS", "HEMODIALYSIS", "OPHTHALMIC", "PERCUTANEOUS", "INTRADERMAL", "INTRAMUSCULAR", "INHALATION"],
    "SOLUTION / DROPS": ["OPHTHALMIC", "ORAL", "OTIC", "TOPICAL"],
    "SOLUTION, CONCENTRATE": ["ORAL", "INTRAVENOUS", "HEMODIALYSIS", "INHALATION", "TOPICAL", "IRRIGATION", "OPHTHALMIC"],
    "SOLUTION, GEL FORMING / DROPS": ["OPHTHALMIC", "TOPICAL", "CUTANEOUS"],
    "SOLUTION, GEL FORMING, EXTENDED RELEASE": ["OPHTHALMIC", "SUBCUTANEOUS", "TOPICAL", "INFILTRATION"],
    
    # Suspensions
    "SUSPENSION": ["ORAL", "INTRAMUSCULAR", "TOPICAL", "OPHTHALMIC", "OTIC", "INHALATION", "RECTAL", "SUBCUTANEOUS"],
    "SUSPENSION / DROPS": ["OPHTHALMIC", "ORAL", "OTIC"],
    "SUSPENSION, EXTENDED RELEASE": ["ORAL", "INTRAMUSCULAR", "SUBCUTANEOUS"],
    
    # Powders - ambiguous
    "POWDER": ["TOPICAL", "ORAL", "INHALATION"],
    "POWDER, FOR SOLUTION": ["ORAL", "INTRAVENOUS", "INTRAMUSCULAR", "SUBCUTANEOUS", "TOPICAL", "NASOGASTRIC"],
    "POWDER, FOR SUSPENSION": ["ORAL", "RECTAL", "SUBCUTANEOUS"],
    "POWDER, METERED": ["INHALATION", "ORAL", "SUBCUTANEOUS"],
    "POWDER, DENTIFRICE": ["ORAL", "DENTAL"],
    
    # Granules
    "GRANULE": ["ORAL", "TOPICAL"],
    "GRANULE, DELAYED RELEASE": ["ORAL", "EXTRACORPOREAL"],
    "GRANULE, EFFERVESCENT": ["ORAL", "TOPICAL"],
    "GRANULE, FOR SOLUTION": ["ORAL", "TOPICAL", "NASAL"],
    "GRANULE, FOR SUSPENSION": ["ORAL"],
    
    # Aerosols - AMBIGUOUS (can be topical, inhalation, nasal, etc.)
    "AEROSOL": ["TOPICAL", "INHALATION", "DENTAL", "ORAL"],
    "AEROSOL, FOAM": ["TOPICAL", "DENTAL", "RECTAL"],
    "AEROSOL, METERED": ["INHALATION", "NASAL", "ORAL", "SUBLINGUAL", "TOPICAL"],
    "AEROSOL, POWDER": ["TOPICAL", "INHALATION", "ORAL"],
    
    # Sprays - AMBIGUOUS
    "SPRAY": ["TOPICAL", "NASAL", "ORAL"],
    "SPRAY, METERED": ["NASAL", "TOPICAL", "INHALATION"],
    "SPRAY, SUSPENSION": ["TOPICAL", "NASAL", "INTRASINAL"],
    
    # Inhalation forms - mostly unambiguous
    "INHALANT": ["INHALATION", "ORAL", "NASAL", "TOPICAL"],
    "INHALER": ["INHALATION"],
    "GAS": ["INHALATION"],
    "NEBULE": ["INHALATION"],
    
    # Topical forms
    "CREAM": ["TOPICAL", "CUTANEOUS"],
    "OINTMENT": ["TOPICAL", "OPHTHALMIC", "RECTAL"],
    "GEL": ["TOPICAL", "DENTAL", "EXTRACORPOREAL", "ORAL"],
    "GEL, METERED": ["TOPICAL", "TRANSDERMAL", "NASAL", "VAGINAL"],
    "GEL, DENTIFRICE": ["DENTAL", "ORAL", "TOPICAL"],
    "LOTION": ["TOPICAL"],
    "LOTION, AUGMENTED": ["TOPICAL"],
    "PASTE": ["TOPICAL", "DENTAL", "ORAL", "CUTANEOUS"],
    "PASTE, DENTIFRICE": ["DENTAL", "ORAL", "TOPICAL"],
    "FOAM": ["TOPICAL", "RECTAL", "VAGINAL"],
    "EMULSION": ["TOPICAL", "INTRAVENOUS", "ORAL", "OPHTHALMIC"],
    "SHAMPOO": ["TOPICAL", "CUTANEOUS"],
    "SOAP": ["TOPICAL", "CUTANEOUS"],
    "OIL": ["TOPICAL", "CUTANEOUS", "ORAL", "OTIC", "TRANSDERMAL"],
    "LINIMENT": ["TOPICAL", "PERCUTANEOUS", "TRANSDERMAL"],
    "TINCTURE": ["TOPICAL", "ORAL", "INHALATION"],
    "LIQUID": ["TOPICAL", "ORAL", "HEMODIALYSIS", "INTRAVENOUS", "OPHTHALMIC", "INTRAMUSCULAR", "SUBCUTANEOUS"],
    "SWAB": ["TOPICAL", "ORAL"],
    "CLOTH": ["TOPICAL", "EXTRACORPOREAL"],
    "SPONGE": ["TOPICAL", "ORAL", "VAGINAL"],
    "DRESSING": ["TOPICAL", "CUTANEOUS", "ORAL"],
    "STICK": ["TOPICAL"],
    "SALVE": ["TOPICAL"],
    "JELLY": ["TOPICAL", "VAGINAL", "NASAL"],
    
    # Transdermal/Patches
    "PATCH": ["TOPICAL", "TRANSDERMAL", "CUTANEOUS", "PERCUTANEOUS"],
    "PATCH, EXTENDED RELEASE": ["TRANSDERMAL", "TOPICAL"],
    "PLASTER": ["TOPICAL", "TRANSDERMAL", "CUTANEOUS"],
    "POULTICE": ["TOPICAL", "TRANSDERMAL", "CUTANEOUS", "PERCUTANEOUS"],
    "FILM": ["BUCCAL", "SUBLINGUAL", "TOPICAL", "ORAL", "DENTAL"],
    "FILM, SOLUBLE": ["BUCCAL", "SUBLINGUAL", "ORAL", "VAGINAL", "TOPICAL"],
    "FILM, EXTENDED RELEASE": ["TRANSDERMAL", "TOPICAL"],
    
    # Drops
    "DROP": ["OPHTHALMIC", "ORAL", "OTIC", "TOPICAL"],
    "DROPS": ["OPHTHALMIC", "ORAL", "OTIC", "TOPICAL"],
    
    # Injections
    "INJECTION": ["INTRAVENOUS", "INTRAMUSCULAR", "SUBCUTANEOUS", "INTRADERMAL", "CUTANEOUS", "INTRATHECAL", "PARENTERAL"],
    "INJECTION, EMULSION": ["INTRAVENOUS", "SUBCUTANEOUS", "INTRAMUSCULAR", "PARENTERAL", "TOPICAL"],
    "INJECTION, LIPID COMPLEX": ["INTRAVENOUS", "EPIDURAL", "INTRATHECAL"],
    "INJECTION, POWDER, FOR SOLUTION": ["INTRAVENOUS", "INTRAMUSCULAR", "SUBCUTANEOUS", "ORAL"],
    "INJECTION, POWDER, FOR SUSPENSION": ["INTRAMUSCULAR", "SUBCUTANEOUS", "INTRAVENOUS"],
    "INJECTION, POWDER, LYOPHILIZED, FOR SOLUTION": ["INTRAVENOUS", "INTRAMUSCULAR", "SUBCUTANEOUS"],
    "INJECTION, POWDER, LYOPHILIZED, FOR SUSPENSION": ["INTRAVENOUS", "SUBCUTANEOUS", "INTRAMUSCULAR"],
    "INJECTION, SOLUTION": ["SUBCUTANEOUS", "INTRADERMAL", "INTRAVENOUS", "INTRAMUSCULAR", "PERCUTANEOUS", "INFILTRATION", "PERINEURAL"],
    "INJECTION, SOLUTION, CONCENTRATE": ["INTRAVENOUS", "INTRAMUSCULAR", "INTRAOCULAR"],
    "INJECTION, SUSPENSION": ["INTRAMUSCULAR", "SUBCUTANEOUS", "INTRA-ARTICULAR", "INTRALESIONAL"],
    "INJECTION, SUSPENSION, EXTENDED RELEASE": ["INTRAMUSCULAR", "SUBCUTANEOUS"],
    "INJECTABLE, LIPOSOMAL": ["INTRAVENOUS"],
    
    # Ampules/Vials - generic parenteral
    "AMPULE": ["PARENTERAL", "INTRAVENOUS", "INTRAMUSCULAR", "SUBCUTANEOUS"],
    "VIAL": ["PARENTERAL", "INTRAVENOUS", "INTRAMUSCULAR", "SUBCUTANEOUS"],
    
    # Containers
    "BOTTLE": ["INTRAVENOUS", "ORAL"],
    "BAG": ["INTRAVENOUS"],
    
    # Rectal forms
    "SUPPOSITORY": ["RECTAL", "VAGINAL", "URETHRAL"],
    "ENEMA": ["RECTAL"],
    
    # Vaginal forms
    "INSERT": ["VAGINAL", "OPHTHALMIC", "CONJUNCTIVAL"],
    "INSERT, EXTENDED RELEASE": ["VAGINAL", "INTRAUTERINE", "PERIODONTAL", "TOPICAL"],
    "RING": ["VAGINAL"],
    "OVULE": ["VAGINAL"],
    
    # Implants
    "IMPLANT": ["SUBCUTANEOUS", "INTRAVITREAL", "INTRA-ARTICULAR", "INTRAOCULAR", "PARENTERAL", "INTRADERMAL"],
    "PELLET, IMPLANTABLE": ["SUBCUTANEOUS", "PARENTERAL"],
    
    # Other
    "KIT": ["ORAL", "TOPICAL", "OPHTHALMIC", "INHALATION", "INTRAVENOUS", "INTRAMUSCULAR", "INFILTRATION", "SUBCUTANEOUS", "EPIDURAL", "INTRA-ARTICULAR"],
    "STRIP": ["ORAL", "OPHTHALMIC", "TOPICAL", "DENTAL"],
    "WAFER": ["ORAL", "INTRACAVITARY", "INTRALESIONAL"],
    "PELLET": ["ORAL", "SUBCUTANEOUS", "OPHTHALMIC", "TOPICAL"],
    "CONCENTRATE": ["ORAL", "INTRADERMAL", "PERCUTANEOUS", "SUBCUTANEOUS", "TOPICAL"],
    "IRRIGANT": ["IRRIGATION", "URETHRAL", "TOPICAL", "INTRAVESICAL", "OPHTHALMIC"],
    "BEAD": ["CUTANEOUS", "TOPICAL", "ORAL"],
}

# ============================================================================
# FORM EQUIVALENCE GROUPS
# Forms within the same group are considered interchangeable for matching
# ============================================================================

FORM_EQUIVALENCE_GROUPS: List[Set[str]] = [
    # Solid oral forms - generally interchangeable
    {"TABLET", "CAPSULE", "CAPLET"},
    # Liquid oral forms
    {"SOLUTION", "SYRUP", "ELIXIR"},
    # Suspensions (separate - different preparation)
    {"SUSPENSION"},
    # Topical semi-solids
    {"CREAM", "OINTMENT", "GEL"},
    # Injectable forms
    {"INJECTION", "AMPULE", "VIAL"},
    # Inhalation forms
    {"INHALER", "AEROSOL", "MDI", "NEBULE", "DPI"},
    # Eye preparations
    {"DROPS"},
]

# Build lookup: form -> set of equivalent forms
FORM_EQUIVALENTS: Dict[str, Set[str]] = {}
for _group in FORM_EQUIVALENCE_GROUPS:
    for _form in _group:
        FORM_EQUIVALENTS[_form] = _group

# ============================================================================
# UNIT/MEASUREMENT TOKENS
# Merged from: UNIT_TOKENS, MEASUREMENT_TOKENS, _PREFIX_UNIT_TOKENS
# ============================================================================

UNIT_TOKENS: Set[str] = {
    # Weight units
    "MG", "G", "MCG", "UG", "KG", "GMS", "GM",
    
    # Volume units
    "ML", "L", "CC",
    
    # Activity units
    "IU", "UNIT", "UNITS", "LSU", "MU",
    
    # Concentration units
    "MEQ", "MEQS", "MOL", "MMOL",
    
    # Percentage
    "PCT", "%",
    
    # Compound units
    "MG/ML", "MCG/ML", "IU/ML", "MG/5ML", "MG/L",
}

UNIT_TOKENS_LOWER: Set[str] = {u.lower() for u in UNIT_TOKENS}

# Weight unit conversion factors (to mg)
WEIGHT_UNIT_FACTORS: Dict[str, float] = {
    "MG": 1.0,
    "G": 1000.0,
    "MCG": 0.001,
    "UG": 0.001,
    "KG": 1_000_000.0,
}

# ============================================================================
# ATC COMBINATION PATTERNS
# ATC codes that indicate combination products
# Merged from: ATC_COMBINATION_PATTERNS, COMBINATION_ATC_PATTERNS
# ============================================================================

ATC_COMBINATION_PATTERNS: List[str] = [
    # Cardiovascular combinations
    "C09DA", "C09DB", "C09DX",  # ARBs + diuretics/CCBs
    "C09BA", "C09BB", "C09BX",  # ACE inhibitors + combos
    "C07FB", "C07BB", "C07CB",  # Beta-blockers + combos
    "C10BA", "C10BX",           # Statins + combos
    
    # Diabetes combinations
    "A10BD",                    # Blood glucose lowering combos
    
    # Pain combinations
    "N02AA55", "N02AA59",       # Opioid combinations
    "N02AJ",                    # Opioid + non-opioid combos
    "N02BE51", "N02BE71",       # Paracetamol combinations
    
    # Antibiotic combinations
    "J01CR", "J01RA",           # Antibiotic combinations
    
    # Respiratory combinations
    "R03AL",                    # Adrenergics + anticholinergics
    "R03AK",                    # Adrenergics + corticosteroids
    "R03DA20", "R03DA55",       # Xanthine combinations
    "R03DB",                    # Xanthines + adrenergics
    
    # Other combinations
    "A02BD",                    # H. pylori eradication combos
    "M05BB",                    # Bisphosphonates combinations
]

# ATC codes ending in these suffixes are typically combinations
COMBINATION_ATC_SUFFIXES: Set[str] = {
    "20", "30", "50", "51", "52", "53", "54", "55", "56", "57", "58", "59"
}

# ============================================================================
# CONNECTIVE WORDS
# Words that connect multiple ingredients in combination products
# ============================================================================

CONNECTIVE_WORDS: Set[str] = {
    "AND", "WITH", "PLUS", "+", "/", "&", "IN",
}

# Tokens that break salt tails (from text_utils_drugs.py)
SALT_TAIL_BREAK_TOKENS: Set[str] = {"+", "/", "&", "AND", "WITH"}

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

def get_canonical_form(form: str) -> str:
    """Get canonical form name, or return uppercase original if not found."""
    return FORM_CANON.get(form.upper(), form.upper())


def get_canonical_route(route: str) -> str:
    """Get canonical route name, or return uppercase original if not found."""
    return ROUTE_CANON.get(route.upper(), route.upper())


def is_stopword(token: str) -> bool:
    """Check if token is a stopword (case-insensitive)."""
    return token.upper() in STOPWORDS or token.lower() in STOPWORDS_LOWER


def is_salt_token(token: str) -> bool:
    """Check if token is a salt/hydrate suffix (case-insensitive)."""
    return token.upper() in SALT_TOKENS


def is_pure_salt_compound(name: str) -> bool:
    """Check if compound name is a pure salt (should not have salt stripped)."""
    return name.upper() in PURE_SALT_COMPOUNDS


def is_element_drug(token: str) -> bool:
    """Check if token is an element that can be a standalone drug."""
    return token.upper() in ELEMENT_DRUGS


def is_unit_token(token: str) -> bool:
    """Check if token is a unit/measurement token."""
    return token.upper() in UNIT_TOKENS or token.lower() in UNIT_TOKENS_LOWER


def is_combination_atc(atc_code: str) -> bool:
    """Check if ATC code indicates a combination product."""
    if not atc_code:
        return False
    atc_upper = atc_code.upper()
    # Check pattern prefixes
    for pattern in ATC_COMBINATION_PATTERNS:
        if atc_upper.startswith(pattern):
            return True
    # Check suffix patterns (last 2 digits)
    if len(atc_upper) >= 2:
        suffix = atc_upper[-2:]
        if suffix in COMBINATION_ATC_SUFFIXES:
            return True
    return False


def forms_are_equivalent(form1: str, form2: str) -> bool:
    """Check if two forms are pharmaceutically equivalent."""
    f1 = get_canonical_form(form1)
    f2 = get_canonical_form(form2)
    if f1 == f2:
        return True
    equiv_set = FORM_EQUIVALENTS.get(f1)
    if equiv_set and f2 in equiv_set:
        return True
    return False


def infer_route_from_form(form: str) -> str | None:
    """Infer the most common route from form. Returns first (most common) route."""
    canonical = get_canonical_form(form)
    # Try FORM_TO_ROUTES first (more complete)
    if canonical in FORM_TO_ROUTES:
        return FORM_TO_ROUTES[canonical][0]
    # Fall back to FORM_TO_ROUTE
    return FORM_TO_ROUTE.get(canonical) or FORM_TO_ROUTE.get(form.upper())


def get_valid_routes_for_form(form: str) -> List[str]:
    """Get ALL valid routes for a form (for exploring matching options).
    
    Returns list of routes ordered by frequency (most common first).
    If form not found, returns empty list.
    """
    canonical = get_canonical_form(form)
    if canonical in FORM_TO_ROUTES:
        return FORM_TO_ROUTES[canonical]
    # Fall back to single-route mapping
    single = FORM_TO_ROUTE.get(canonical) or FORM_TO_ROUTE.get(form.upper())
    return [single] if single else []


def is_valid_form_route_pair(form: str, route: str) -> bool:
    """Check if a form-route combination is valid according to DrugBank data."""
    valid_routes = get_valid_routes_for_form(form)
    if not valid_routes:
        return True  # Unknown form, allow any route
    route_upper = route.upper()
    return route_upper in valid_routes


# ============================================================================
# GARBAGE TOKENS - Tokens to filter out when extracting generic names
# These are not drug names but formulation/packaging/flavor words
# ============================================================================

GARBAGE_TOKENS: Set[str] = {
    # Units
    'MG', 'ML', 'MCG', 'G', 'IU', 'UNIT', 'UNITS',
    # Dosage forms
    'TAB', 'TABLET', 'CAP', 'CAPSULE', 'AMP', 'AMPULE', 'VIAL', 'BOTTLE',
    # Routes
    'ORAL', 'IV', 'IM', 'SC', 'TOPICAL',
    # Marketing descriptors
    'FORTE', 'PLUS', 'EXTRA', 'MAX', 'ULTRA', 'JUNIOR', 'PEDIA', 'ADULT',
    # Flavors
    'ORANGE', 'STRAWBERRY', 'CHERRY', 'GRAPE', 'MINT', 'VANILLA', 'LEMON',
    # Junk tokens
    'PNF', 'NAN', '-', '+', '/', 'AND', 'WITH',
    # Formulation words
    'SOLVENT', 'DILUENT', 'SOLUTION', 'SUSPENSION', 'POWDER',
}

# ============================================================================
# GENERIC SYNONYMS - Drug name synonym mappings
# Bidirectional mappings for matching drugs with different names
# ============================================================================

GENERIC_SYNONYMS: Dict[str, str] = {
    # INN/BAN/USAN variations
    'ACETAMINOPHEN': 'PARACETAMOL',
    'PARACETAMOL': 'ACETAMINOPHEN',
    'SALBUTAMOL': 'ALBUTEROL',
    'ALBUTEROL': 'SALBUTAMOL',
    'ADRENALINE': 'EPINEPHRINE',
    'EPINEPHRINE': 'ADRENALINE',
    'NORADRENALINE': 'NOREPINEPHRINE',
    'NOREPINEPHRINE': 'NORADRENALINE',
    'LIGNOCAINE': 'LIDOCAINE',
    'LIDOCAINE': 'LIGNOCAINE',
    'FRUSEMIDE': 'FUROSEMIDE',
    'FUROSEMIDE': 'FRUSEMIDE',
    'BENZYLPENICILLIN': 'PENICILLIN G',
    'PENICILLIN G': 'BENZYLPENICILLIN',
    
    # DrugBank synonyms found in Annex F
    'CEFALEXIN': 'CEPHALEXIN',
    'CEPHALEXIN': 'CEFALEXIN',
    'CHLORPHENAMINE': 'CHLORPHENIRAMINE',
    'CHLORPHENIRAMINE': 'CHLORPHENAMINE',
    'CICLOSPORIN': 'CYCLOSPORINE',
    'CYCLOSPORINE': 'CICLOSPORIN',
    'DICYCLOVERINE': 'DICYCLOMINE',
    'DICYCLOMINE': 'DICYCLOVERINE',
    'GLIBENCLAMIDE': 'GLYBURIDE',
    'GLYBURIDE': 'GLIBENCLAMIDE',
    'MECLOZINE': 'MECLIZINE',
    'MECLIZINE': 'MECLOZINE',
    'PROXYMETACAINE': 'PROPARACAINE',
    'PROPARACAINE': 'PROXYMETACAINE',
    'THIAMAZOLE': 'METHIMAZOLE',
    'METHIMAZOLE': 'THIAMAZOLE',
    
    # Combination drug synonyms - canonical form is AMOXICILLIN + CLAVULANIC ACID
    'CO-AMOXICLAV': 'AMOXICILLIN + CLAVULANIC ACID',
    'AMOXICILLIN-CLAVULANIC ACID': 'AMOXICILLIN + CLAVULANIC ACID',
    'AMOXICILLIN AND CLAVULANATE POTASSIUM': 'AMOXICILLIN + CLAVULANIC ACID',
    'AUGMENTIN': 'AMOXICILLIN + CLAVULANIC ACID',
    
    # Aluminum/Aluminium spelling
    'ALUMINUM HYDROXIDE': 'ALUMINIUM HYDROXIDE',
    'ALUMINIUM HYDROXIDE': 'ALUMINUM HYDROXIDE',
    'MAGNESIUM HYDROXIDE': 'MAGNESIUM',
}

# ============================================================================
# IV FLUID SYNONYMS - Abbreviations for IV fluids
# ============================================================================

IV_FLUID_SYNONYMS: Dict[str, str] = {
    # NOTE: Main synonyms are in SPELLING_SYNONYMS - this dict is for compatibility
    # D5 = 5% Dextrose, D10 = 10% Dextrose (concentration captured separately)
    'D5': 'GLUCOSE',            # 5% Dextrose
    'D5W': 'GLUCOSE',           # 5% Dextrose in Water
    'D10': 'GLUCOSE',           # 10% Dextrose
    'D10W': 'GLUCOSE',          # 10% Dextrose in Water
    'DEXTROSE': 'GLUCOSE',
    'NSS': 'SODIUM CHLORIDE',   # Normal Saline Solution (0.9%)
    'PNSS': 'SODIUM CHLORIDE',  # Plain Normal Saline Solution
    'NORMAL SALINE': 'SODIUM CHLORIDE',
    'LR': "LACTATED RINGER'S",
    'LRS': "LACTATED RINGER'S",
    'D5LR': 'GLUCOSE + LACTATED RINGER\'S',   # 5% Dextrose + LR
    'D5NSS': 'GLUCOSE + SODIUM CHLORIDE',     # 5% Dextrose + NSS
    'STERILE WATER': 'WATER FOR INJECTION',
}

# ============================================================================
# DRUGBANK COMPONENT SYNONYMS - Map DrugBank component names to parent drug
# These are chemical components that should map to the marketed drug mixture
# ============================================================================

DRUGBANK_COMPONENT_SYNONYMS: Dict[str, str] = {
    # Gentamicin is a mixture of C1, C1a, C2 components
    'GENTAMICIN C2': 'GENTAMICIN',
    'GENTAMICIN C1': 'GENTAMICIN',
    'GENTAMICIN C1A': 'GENTAMICIN',
    'GENTAMICIN': 'GENTAMICIN C2',  # Reverse: search GENTAMICIN finds GENTAMICIN C2
}

# Combined synonym lookup (for convenience)
ALL_DRUG_SYNONYMS: Dict[str, str] = {
    **GENERIC_SYNONYMS,
    **IV_FLUID_SYNONYMS,
    **DRUGBANK_COMPONENT_SYNONYMS,
}

# ============================================================================
# SPELLING SYNONYMS - Not regional variants, just spelling corrections
# These are ONLY for spelling variants that cannot be in the unified dataset
# ============================================================================

SPELLING_SYNONYMS: Dict[str, str] = {
    # Spelling variants (different spellings of same drug)
    "POLYMYXIN": "POLYMIXIN",  # X vs I spelling
    "POLYMYXIN B": "POLYMIXIN B",
    # Common/trade names that are really generics
    "ASPIRIN": "ACETYLSALICYLIC ACID",
    "ASA": "ACETYLSALICYLIC ACID",
    # Combination drug trade names -> canonical generic form
    # CO-AMOXICLAV = Amoxicillin + Clavulanic acid (J01CR02)
    "CO-AMOXICLAV": "AMOXICILLIN + CLAVULANIC ACID",
    "AUGMENTIN": "AMOXICILLIN + CLAVULANIC ACID",
    "AMOXICILLIN AND CLAVULANATE POTASSIUM": "AMOXICILLIN + CLAVULANIC ACID",
    "AMOXICILLIN AND CLAVULANIC ACID": "AMOXICILLIN + CLAVULANIC ACID",
    # Cotrimoxazole = Sulfamethoxazole + Trimethoprim (J01EE01)
    "COTRIMOXAZOLE": "SULFAMETHOXAZOLE + TRIMETHOPRIM",
    "CO-TRIMOXAZOLE": "SULFAMETHOXAZOLE + TRIMETHOPRIM",
    "BACTRIM": "SULFAMETHOXAZOLE + TRIMETHOPRIM",
    # TB drug combinations (first-line treatment for drug-susceptible TB)
    # HRZE = H(Isoniazid) + R(Rifampicin) + Z(Pyrazinamide) + E(Ethambutol)
    # Only use 3+ letter abbreviations to avoid false matches (HR/HE too short)
    "HRZE": "ISONIAZID + RIFAMPICIN + PYRAZINAMIDE + ETHAMBUTOL",
    "HRZ": "ISONIAZID + RIFAMPICIN + PYRAZINAMIDE",
    "HRE": "ISONIAZID + RIFAMPICIN + ETHAMBUTOL",
    "RHZ": "ISONIAZID + RIFAMPICIN + PYRAZINAMIDE",  # Alternate order
    "RHZE": "ISONIAZID + RIFAMPICIN + PYRAZINAMIDE + ETHAMBUTOL",  # Alternate order
    "RHEZ": "ISONIAZID + RIFAMPICIN + PYRAZINAMIDE + ETHAMBUTOL",  # Alternate order
    # Salt form synonyms
    "POTASSIUM CLAVULANATE": "CLAVULANIC ACID",
    "CLAVULANATE": "CLAVULANIC ACID",
    # Vitamin synonyms
    "ALPHA-TOCOPHEROL": "VITAMIN E",
    "TOCOPHEROL": "VITAMIN E",
    # Common abbreviations
    "VIT": "VITAMIN",
    "VIT A": "VITAMIN A",
    "VIT B": "VITAMIN B",
    "VIT C": "VITAMIN C",
    "VIT D": "VITAMIN D",
    "VIT E": "VITAMIN E",
    "VIT K": "VITAMIN K",
    # Regional/INN name synonyms (map to DrugBank name for lookup)
    "HYOSCINE-N-BUTYLBROMIDE": "BUTYLSCOPOLAMINE",
    "HYOSCINE BUTYLBROMIDE": "BUTYLSCOPOLAMINE",
    "HYOSCINE N-BUTYLBROMIDE": "BUTYLSCOPOLAMINE",
    "SCOPOLAMINE BUTYLBROMIDE": "BUTYLSCOPOLAMINE",
    "DIMETICONE": "SIMETHICONE",  # INN name → US name
    # Common IV solution abbreviations
    # Note: D5 = 5% Dextrose, D10 = 10% Dextrose - concentration is in dose field
    # DEXTROSE/D-GLUCOSE = GLUCOSE (same compound, WHO ATC uses "GLUCOSE")
    "PNSS": "SODIUM CHLORIDE",  # Plain Normal Saline Solution
    "NSS": "SODIUM CHLORIDE",   # Normal Saline Solution
    "0.9% SODIUM CHLORIDE": "SODIUM CHLORIDE",  # Normal saline
    "0.45% SODIUM CHLORIDE": "SODIUM CHLORIDE",  # Half-normal saline
    "D5": "DEXTROSE",           # 5% Dextrose (concentration in dose)
    "D5W": "DEXTROSE",          # 5% Dextrose in Water
    "D10": "DEXTROSE",          # 10% Dextrose
    "D10W": "DEXTROSE",         # 10% Dextrose in Water
    "5% DEXTROSE": "DEXTROSE",  # Generic with percentage
    "10% DEXTROSE": "DEXTROSE",
    "50% DEXTROSE": "DEXTROSE",
    "DEXTROSE": "DEXTROSE",     # Keep as DEXTROSE (matched in PNF)
    "D-GLUCOSE": "DEXTROSE",
    "GLUCOSE": "DEXTROSE",      # WHO ATC uses GLUCOSE, PNF uses DEXTROSE
    "LRS": "LACTATED RINGER'S", # Lactated Ringer's Solution
    "D5LR": "DEXTROSE + LACTATED RINGER'S",
    "D5NSS": "DEXTROSE + SODIUM CHLORIDE",
    "5% DEXTROSE IN WATER": "DEXTROSE",
    "10% DEXTROSE IN WATER": "DEXTROSE",
    "5% DEXTROSE IN 0.9% SODIUM CHLORIDE": "DEXTROSE + SODIUM CHLORIDE",
    "5% DEXTROSE IN 0.3% SODIUM CHLORIDE": "DEXTROSE + SODIUM CHLORIDE",
    "5% DEXTROSE IN LACTATED RINGER'S": "DEXTROSE + LACTATED RINGER'S",
    # Brand name mappings
    "REBAMID": "REBAMIPIDE",
    # Alcohol synonyms
    "ETHYL ALCOHOL": "ETHANOL",
    "ALCOHOL ETHYL": "ETHANOL",  # Annex F format: "ALCOHOL, ETHYL"
    "ISOPROPYL ALCOHOL": "ISOPROPANOL",
    "ALCOHOL ISOPROPYL": "ISOPROPANOL",
    # IV Solutions - Philippine formulary names
    "ACETATED RINGER'S SOLUTION": "RINGER'S SOLUTION, ACETATED",
    "ACETATED RINGER'S": "RINGER'S SOLUTION, ACETATED",
    "RINGER'S SOLUTION ACETATED": "RINGER'S SOLUTION, ACETATED",
    "LACTATED RINGER'S SOLUTION": "RINGER'S SOLUTION, LACTATED",
    "LACTATED RINGER'S": "RINGER'S SOLUTION, LACTATED",
    "RINGER'S LACTATE": "RINGER'S SOLUTION, LACTATED",
    # Amino acid solutions (parenteral nutrition)
    "AMINO ACID SOLUTIONS FOR HEPATIC FAILURE": "AMINO ACIDS",
    "AMINO ACID SOLUTIONS FOR INFANTS": "AMINO ACIDS",
    "AMINO ACID SOLUTIONS FOR RENAL CONDITIONS": "AMINO ACIDS",
    "AMINO ACID SOLUTIONS FOR IMMUNONUTRITION": "AMINO ACIDS",
    "AMINO ACID SOLUTIONS": "AMINO ACIDS",
    "AMINO ACID SOLUTION": "AMINO ACIDS",
    # Combination drugs - Annex F format
    "ALUMINUM HYDROXIDE + MAGNESIUM HYDROXIDE": "ALUMINIUM HYDROXIDE + MAGNESIUM HYDROXIDE",
    "AMPICILLIN + SULBACTAM": "SULTAMICILLIN",  # Combination has ATC J01CR04
    "ALENDRONATE + CHOLECALCIFEROL": "ALENDRONATE + VITAMIN D",
    "CHOLECALCIFEROL": "VITAMIN D3",
    "VIT. D3": "VITAMIN D3",
    # Philippine herbal medicine
    "AKAPULKO": "CASSIA ALATA",
    # Vitamin combinations
    "VITAMINS INTRAVENOUS FAT-SOLUBLE": "VITAMINS A, D, E AND K",
    "VITAMINS INTRAVENOUS WATER-SOLUBLE": "VITAMIN B COMPLEX + VITAMIN C",
    "VITAMIN INTRAVENOUS FAT-SOLUBLE": "VITAMINS A, D, E AND K",
    "VITAMIN INTRAVENOUS WATER-SOLUBLE": "VITAMIN B COMPLEX + VITAMIN C",
}

# ============================================================================
# VACCINE CANONICAL NAMES - Normalize vaccine descriptions to canonical form
# Format: "canonical_name" -> list of patterns/aliases
# Output: generic_name = canonical, details = specifics (valency, strains, etc.)
# ============================================================================

VACCINE_CANONICAL: Dict[str, Dict[str, Any]] = {
    # BCG
    "BCG VACCINE": {
        "patterns": ["BCG VACCINE", "BACILLUS CALMETTE-GUERIN", "BACILLUS CALMETTE GUERIN"],
        "acronym": "BCG",
    },
    # Hepatitis
    "HEPATITIS A VACCINE": {
        "patterns": ["HEPATITIS A INACTIVATED VACCINE", "HEPATITIS A VACCINE"],
        "acronym": "HEPA",
    },
    "HEPATITIS B VACCINE": {
        "patterns": ["HEPATITIS B VACCINE", "HEPATITIS B RECOMBINANT"],
        "acronym": "HEPB",
    },
    "HEPATITIS A + B VACCINE": {
        "patterns": ["HEPATITIS A + B VACCINE", "HEPATITIS A AND B VACCINE"],
        "acronym": "HEPA+HEPB",
    },
    # DTP combinations
    "DTP VACCINE": {
        "patterns": ["DIPHTHERIA-TETANUS TOXOIDS AND PERTUSSIS VACCINE", 
                     "DIPHTHERIA-TETANUS TOXOIDS AND ACELLULAR PERTUSSIS VACCINE",
                     "DIPHTHERIA, TETANUS, PERTUSSIS"],
        "acronym": "DTP",
        "aliases": ["DTP", "DTAP"],
    },
    "DT VACCINE": {
        "patterns": ["DIPHTHERIA-TETANUS TOXOIDS"],
        "acronym": "DT",
    },
    "DTP + HIB VACCINE": {
        "patterns": ["DTP + HIB", "DTAP + HIB", "DTP-HIB"],
        "acronym": "DTP-HIB",
    },
    "DTP + HEPATITIS B VACCINE": {
        "patterns": ["DTP + HEPATITIS B VACCINE", "DTAP + HEPATITIS B"],
        "acronym": "DTP-HEPB",
    },
    "DTP + IPV VACCINE": {
        "patterns": ["DTP + INACTIVATED POLIO VACCINE", "DTP + IPV", "DTAP + IPV"],
        "acronym": "DTP-IPV",
    },
    "DTP + IPV + HIB VACCINE": {
        "patterns": ["DTP + IPV + HIB", "DTAP + IPV + HIB"],
        "acronym": "DTP-IPV-HIB",
    },
    # Polio
    "IPV VACCINE": {
        "patterns": ["INACTIVATED POLIOMYELITIS VACCINE", "INACTIVATED POLIO VACCINE", "IPV"],
        "acronym": "IPV",
    },
    "OPV VACCINE": {
        "patterns": ["ORAL POLIO VACCINE", "LIVE ATTENUATED TRIVALENT ORAL POLIO VACCINE", "OPV"],
        "acronym": "OPV",
    },
    # Measles/MMR
    "MEASLES VACCINE": {
        "patterns": ["LIVE ATTENUATED MEASLES VACCINE", "MEASLES VACCINE"],
        "acronym": "MV",
    },
    "MUMPS VACCINE": {
        "patterns": ["LIVE ATTENUATED MUMPS VACCINE", "MUMPS VACCINE"],
        "acronym": "MuV",
    },
    "RUBELLA VACCINE": {
        "patterns": ["LIVE ATTENUATED RUBELLA VACCINE", "RUBELLA VACCINE"],
        "acronym": "RV",
    },
    "MMR VACCINE": {
        "patterns": ["LIVE ATTENUATED MEASLES, MUMPS, AND RUBELLA", "MMR VACCINE", 
                     "MEASLES, MUMPS, AND RUBELLA VACCINE", "MMR"],
        "acronym": "MMR",
    },
    "VARICELLA VACCINE": {
        "patterns": ["LIVE ATTENUATED VARICELLA VACCINE", "VARICELLA VACCINE", "CHICKENPOX VACCINE"],
        "acronym": "VZV",
    },
    # Pneumococcal
    "PNEUMOCOCCAL VACCINE": {
        "patterns": ["PNEUMOCOCCAL CONJUGATE VACCINE", "PNEUMOCOCCAL POLYVALENT VACCINE",
                     "PNEUMOCOCCAL POLYSACCHARIDE"],
        "acronym": "PCV",
    },
    # Meningococcal
    "MENINGOCOCCAL VACCINE": {
        "patterns": ["MENINGOCOCCAL POLYSACCHARIDE", "MENINGOCOCCAL CONJUGATE", "NEISSERIA MENINGITIDIS"],
        "acronym": "MenV",
    },
    # HIB
    "HIB VACCINE": {
        "patterns": ["HEMOPHILUS INFLUENZAE TYPE B", "HIB CONJUGATED VACCINE", "HAEMOPHILUS INFLUENZAE"],
        "acronym": "HIB",
    },
    # Influenza
    "INFLUENZA VACCINE": {
        "patterns": ["INFLUENZA VACCINE", "INFLUENZA POLYVALENT VACCINE", "SPLIT VIRION"],
        "acronym": "FLU",
    },
    # Rotavirus
    "ROTAVIRUS VACCINE": {
        "patterns": ["LIVE ATTENUATED ROTAVIRUS", "ROTAVIRUS VACCINE"],
        "acronym": "RV",
    },
    # Rabies
    "RABIES VACCINE": {
        "patterns": ["RABIES CHICK EMBRYO CELL", "RABIES VACCINE", "PURIFIED INACTIVATED"],
        "acronym": "RabV",
    },
    # Yellow fever
    "YELLOW FEVER VACCINE": {
        "patterns": ["YELLOW FEVER VACCINE"],
        "acronym": "YFV",
    },
    # HPV
    "HPV VACCINE": {
        "patterns": ["HUMAN PAPILLOMAVIRUS", "HPV VACCINE", "QUADRIVALENT", "BIVALENT", "NONAVALENT"],
        "acronym": "HPV",
    },
    # Typhoid
    "TYPHOID VACCINE": {
        "patterns": ["TYPHOID VACCINE", "SALMONELLA TYPHI"],
        "acronym": "TyV",
    },
    # Japanese encephalitis
    "JAPANESE ENCEPHALITIS VACCINE": {
        "patterns": ["JAPANESE ENCEPHALITIS VACCINE", "JE VACCINE"],
        "acronym": "JEV",
    },
    # Pentavalent
    "PENTAVALENT VACCINE": {
        "patterns": ["DIPHTHERIA, TETANUS, PERTUSSIS, HEPATITIS B.*HAEMOPHILUS"],
        "acronym": "PENTA",
    },
}

# Helper function to normalize vaccine names
def normalize_vaccine_name(text: str) -> Tuple[Optional[str], Optional[str]]:
    """
    Normalize a vaccine description to canonical name + details.
    
    Returns: (canonical_name, details) or (None, None) if not a vaccine
    """
    text_upper = text.upper()
    
    # Check if it's a vaccine
    if "VACCINE" not in text_upper and "TOXOID" not in text_upper:
        return None, None
    
    # Try to match against patterns
    import re
    for canonical, info in VACCINE_CANONICAL.items():
        for pattern in info["patterns"]:
            if pattern in text_upper or re.search(pattern, text_upper, re.IGNORECASE):
                # Extract details (valency, strains, etc.)
                details = []
                
                # Valency
                valency_match = re.search(r'(\d+)-?VALENT', text_upper)
                if valency_match:
                    details.append(f"{valency_match.group(1)}-valent")
                
                # Type/strain info in parentheses
                type_match = re.search(r'\(TYPE[S]?\s+([^)]+)\)', text_upper)
                if type_match:
                    details.append(f"Type {type_match.group(1)}")
                
                # Group/serogroup
                group_match = re.search(r'(?:GROUP|SEROGROUP)\s+([A-Z,\s\+]+?)(?:\s|$|\))', text_upper)
                if group_match:
                    details.append(f"Group {group_match.group(1).strip()}")
                
                # Recombinant/Attenuated
                if "RECOMBINANT" in text_upper:
                    details.append("Recombinant")
                if "ATTENUATED" in text_upper and "LIVE" in text_upper:
                    details.append("Live attenuated")
                elif "INACTIVATED" in text_upper:
                    details.append("Inactivated")
                
                # Pediatric/Adult
                if "PEDIATRIC" in text_upper or "JUNIOR" in text_upper:
                    details.append("Pediatric")
                elif "ADULT" in text_upper:
                    details.append("Adult")
                
                detail_str = "; ".join(details) if details else None
                return canonical, detail_str
    
    # Fallback: return generic "VACCINE" if contains vaccine
    if "VACCINE" in text_upper:
        return "VACCINE", text_upper.replace("VACCINE", "").strip()
    
    return None, None


# ============================================================================
# VACCINE ACRONYM LOOKUP - Bidirectional mapping for vaccine abbreviations
# WHO/CDC standard abbreviations that enable matching between:
#   - Acronym: "DTP" ↔ Components: "DIPHTHERIA + TETANUS + PERTUSSIS"
# Used to match Annex F (often uses acronyms) with ESOA (often spells out components)
# ============================================================================

# Acronym → List of component antigens (sorted alphabetically for matching)
VACCINE_ACRONYM_TO_COMPONENTS: Dict[str, List[str]] = {
    # Single antigens
    "BCG": ["BACILLUS CALMETTE-GUERIN"],
    "D": ["DIPHTHERIA"],
    "T": ["TETANUS"],
    "P": ["PERTUSSIS"],
    "AP": ["ACELLULAR PERTUSSIS"],
    "WP": ["WHOLE-CELL PERTUSSIS"],
    "HIB": ["HAEMOPHILUS INFLUENZAE TYPE B"],
    "HEPB": ["HEPATITIS B"],
    "HEPA": ["HEPATITIS A"],
    "IPV": ["INACTIVATED POLIO", "INACTIVATED POLIOVIRUS", "INACTIVATED POLIOMYELITIS"],
    "OPV": ["ORAL POLIO", "ORAL POLIOVIRUS", "LIVE ATTENUATED POLIO"],
    "MV": ["MEASLES"],
    "MR": ["MEASLES", "RUBELLA"],
    "MMR": ["MEASLES", "MUMPS", "RUBELLA"],
    "MMRV": ["MEASLES", "MUMPS", "RUBELLA", "VARICELLA"],
    "VAR": ["VARICELLA"],
    "VZV": ["VARICELLA", "VARICELLA-ZOSTER"],
    "RV": ["ROTAVIRUS"],
    "PCV": ["PNEUMOCOCCAL CONJUGATE"],
    "PPSV": ["PNEUMOCOCCAL POLYSACCHARIDE"],
    "FLU": ["INFLUENZA"],
    "IIV": ["INACTIVATED INFLUENZA"],
    "LAIV": ["LIVE ATTENUATED INFLUENZA"],
    "HPV": ["HUMAN PAPILLOMAVIRUS"],
    "YF": ["YELLOW FEVER"],
    "JE": ["JAPANESE ENCEPHALITIS"],
    "RAB": ["RABIES"],
    "TYP": ["TYPHOID"],
    "MEN": ["MENINGOCOCCAL"],
    
    # Two-component combinations
    "DT": ["DIPHTHERIA", "TETANUS"],
    "TD": ["TETANUS", "DIPHTHERIA"],  # Adult formulation
    "DP": ["DIPHTHERIA", "PERTUSSIS"],
    "TP": ["TETANUS", "PERTUSSIS"],
    
    # Three-component combinations (DTP family)
    "DTP": ["DIPHTHERIA", "TETANUS", "PERTUSSIS"],
    "DTAP": ["DIPHTHERIA", "TETANUS", "ACELLULAR PERTUSSIS"],
    "DTWP": ["DIPHTHERIA", "TETANUS", "WHOLE-CELL PERTUSSIS"],
    "TDAP": ["TETANUS", "DIPHTHERIA", "ACELLULAR PERTUSSIS"],  # Adult formulation
    
    # Four-component combinations
    "DTP-HIB": ["DIPHTHERIA", "TETANUS", "PERTUSSIS", "HAEMOPHILUS INFLUENZAE TYPE B"],
    "DTP-HEPB": ["DIPHTHERIA", "TETANUS", "PERTUSSIS", "HEPATITIS B"],
    "DTP-IPV": ["DIPHTHERIA", "TETANUS", "PERTUSSIS", "INACTIVATED POLIO"],
    "DTAP-HIB": ["DIPHTHERIA", "TETANUS", "ACELLULAR PERTUSSIS", "HAEMOPHILUS INFLUENZAE TYPE B"],
    "DTAP-HEPB": ["DIPHTHERIA", "TETANUS", "ACELLULAR PERTUSSIS", "HEPATITIS B"],
    "DTAP-IPV": ["DIPHTHERIA", "TETANUS", "ACELLULAR PERTUSSIS", "INACTIVATED POLIO"],
    
    # Five-component combinations (Pentavalent)
    "PENTA": ["DIPHTHERIA", "TETANUS", "PERTUSSIS", "HEPATITIS B", "HAEMOPHILUS INFLUENZAE TYPE B"],
    "DTP-HEPB-HIB": ["DIPHTHERIA", "TETANUS", "PERTUSSIS", "HEPATITIS B", "HAEMOPHILUS INFLUENZAE TYPE B"],
    "DTP-IPV-HIB": ["DIPHTHERIA", "TETANUS", "PERTUSSIS", "INACTIVATED POLIO", "HAEMOPHILUS INFLUENZAE TYPE B"],
    "DTAP-HEPB-HIB": ["DIPHTHERIA", "TETANUS", "ACELLULAR PERTUSSIS", "HEPATITIS B", "HAEMOPHILUS INFLUENZAE TYPE B"],
    "DTAP-IPV-HIB": ["DIPHTHERIA", "TETANUS", "ACELLULAR PERTUSSIS", "INACTIVATED POLIO", "HAEMOPHILUS INFLUENZAE TYPE B"],
    
    # Six-component combinations (Hexavalent)
    "HEXA": ["DIPHTHERIA", "TETANUS", "PERTUSSIS", "HEPATITIS B", "HAEMOPHILUS INFLUENZAE TYPE B", "INACTIVATED POLIO"],
    "DTP-HEPB-HIB-IPV": ["DIPHTHERIA", "TETANUS", "PERTUSSIS", "HEPATITIS B", "HAEMOPHILUS INFLUENZAE TYPE B", "INACTIVATED POLIO"],
    "DTAP-HEPB-HIB-IPV": ["DIPHTHERIA", "TETANUS", "ACELLULAR PERTUSSIS", "HEPATITIS B", "HAEMOPHILUS INFLUENZAE TYPE B", "INACTIVATED POLIO"],
    
    # Hepatitis combinations
    "HEPA-HEPB": ["HEPATITIS A", "HEPATITIS B"],
    "TWINRIX": ["HEPATITIS A", "HEPATITIS B"],  # Brand name used as generic
    
    # Meningococcal serogroups
    "MENACWY": ["MENINGOCOCCAL A", "MENINGOCOCCAL C", "MENINGOCOCCAL W", "MENINGOCOCCAL Y"],
    "MENB": ["MENINGOCOCCAL B"],
    "MENABCWY": ["MENINGOCOCCAL A", "MENINGOCOCCAL B", "MENINGOCOCCAL C", "MENINGOCOCCAL W", "MENINGOCOCCAL Y"],
    
    # Pneumococcal valency
    "PCV7": ["PNEUMOCOCCAL CONJUGATE 7-VALENT"],
    "PCV10": ["PNEUMOCOCCAL CONJUGATE 10-VALENT"],
    "PCV13": ["PNEUMOCOCCAL CONJUGATE 13-VALENT"],
    "PCV15": ["PNEUMOCOCCAL CONJUGATE 15-VALENT"],
    "PCV20": ["PNEUMOCOCCAL CONJUGATE 20-VALENT"],
    "PPSV23": ["PNEUMOCOCCAL POLYSACCHARIDE 23-VALENT"],
}

# Component keyword → Normalized component name (for text parsing)
VACCINE_COMPONENT_KEYWORDS: Dict[str, str] = {
    # Diphtheria
    "DIPHTHERIA": "DIPHTHERIA",
    "DIPHTERIA": "DIPHTHERIA",  # Common misspelling
    "CORYNEBACTERIUM DIPHTHERIAE": "DIPHTHERIA",
    
    # Tetanus
    "TETANUS": "TETANUS",
    "CLOSTRIDIUM TETANI": "TETANUS",
    
    # Pertussis
    "PERTUSSIS": "PERTUSSIS",
    "WHOOPING COUGH": "PERTUSSIS",
    "BORDETELLA PERTUSSIS": "PERTUSSIS",
    "ACELLULAR PERTUSSIS": "ACELLULAR PERTUSSIS",
    "WHOLE-CELL PERTUSSIS": "WHOLE-CELL PERTUSSIS",
    "WHOLE CELL PERTUSSIS": "WHOLE-CELL PERTUSSIS",
    
    # Haemophilus influenzae type B
    "HAEMOPHILUS INFLUENZAE TYPE B": "HAEMOPHILUS INFLUENZAE TYPE B",
    "HAEMOPHILUS INFLUENZAE B": "HAEMOPHILUS INFLUENZAE TYPE B",
    "H. INFLUENZAE TYPE B": "HAEMOPHILUS INFLUENZAE TYPE B",
    "HIB": "HAEMOPHILUS INFLUENZAE TYPE B",
    
    # Hepatitis
    "HEPATITIS A": "HEPATITIS A",
    "HEPATITIS B": "HEPATITIS B",
    "HEP A": "HEPATITIS A",
    "HEP B": "HEPATITIS B",
    
    # Polio
    "POLIO": "POLIO",
    "POLIOVIRUS": "POLIO",
    "POLIOMYELITIS": "POLIO",
    "INACTIVATED POLIO": "INACTIVATED POLIO",
    "INACTIVATED POLIOVIRUS": "INACTIVATED POLIO",
    "INACTIVATED POLIOMYELITIS": "INACTIVATED POLIO",
    "ORAL POLIO": "ORAL POLIO",
    "LIVE ATTENUATED POLIO": "ORAL POLIO",
    
    # MMR components
    "MEASLES": "MEASLES",
    "RUBEOLA": "MEASLES",
    "MUMPS": "MUMPS",
    "RUBELLA": "RUBELLA",
    "GERMAN MEASLES": "RUBELLA",
    
    # Varicella
    "VARICELLA": "VARICELLA",
    "CHICKENPOX": "VARICELLA",
    "VARICELLA-ZOSTER": "VARICELLA",
    "VARICELLA ZOSTER": "VARICELLA",
    
    # Others
    "ROTAVIRUS": "ROTAVIRUS",
    "INFLUENZA": "INFLUENZA",
    "FLU": "INFLUENZA",
    "PNEUMOCOCCAL": "PNEUMOCOCCAL",
    "MENINGOCOCCAL": "MENINGOCOCCAL",
    "HUMAN PAPILLOMAVIRUS": "HUMAN PAPILLOMAVIRUS",
    "HPV": "HUMAN PAPILLOMAVIRUS",
    "YELLOW FEVER": "YELLOW FEVER",
    "JAPANESE ENCEPHALITIS": "JAPANESE ENCEPHALITIS",
    "RABIES": "RABIES",
    "TYPHOID": "TYPHOID",
    "TUBERCULOSIS": "TUBERCULOSIS",
    "BCG": "TUBERCULOSIS",
    "BACILLUS CALMETTE-GUERIN": "TUBERCULOSIS",
    "BACILLUS CALMETTE GUERIN": "TUBERCULOSIS",
}

def _build_components_to_acronym() -> Dict[str, str]:
    """Build reverse mapping from sorted component key to acronym."""
    result = {}
    for acronym, components in VACCINE_ACRONYM_TO_COMPONENTS.items():
        # Normalize and sort components for order-independent matching
        normalized = sorted([c.upper() for c in components])
        key = " + ".join(normalized)
        # Prefer shorter acronyms (DTP over DTP-HIB-HEPB)
        if key not in result or len(acronym) < len(result[key]):
            result[key] = acronym
    return result

# Reverse mapping: sorted component key → acronym
VACCINE_COMPONENTS_TO_ACRONYM: Dict[str, str] = _build_components_to_acronym()


def normalize_vaccine_components(text: str) -> List[str]:
    """
    Extract and normalize vaccine components from text.
    
    Returns a sorted list of normalized component names.
    """
    text_upper = text.upper()
    components = []
    
    # Check for each component keyword
    for keyword, normalized in sorted(VACCINE_COMPONENT_KEYWORDS.items(), key=lambda x: -len(x[0])):
        if keyword in text_upper:
            if normalized not in components:
                components.append(normalized)
            # Remove the keyword to avoid double-matching
            text_upper = text_upper.replace(keyword, " ")
    
    return sorted(components)


def expand_vaccine_acronym(acronym: str) -> Optional[List[str]]:
    """
    Expand a vaccine acronym to its component antigens.
    
    Args:
        acronym: Vaccine abbreviation (e.g., "DTP", "PENTA")
    
    Returns:
        List of component names, or None if not a recognized acronym
    """
    acronym_upper = acronym.upper().strip()
    # Remove common suffixes
    for suffix in [" VACCINE", "-VACCINE", "VACCINE"]:
        if acronym_upper.endswith(suffix):
            acronym_upper = acronym_upper[:-len(suffix)].strip()
    
    return VACCINE_ACRONYM_TO_COMPONENTS.get(acronym_upper)


def get_vaccine_acronym(components: List[str]) -> Optional[str]:
    """
    Get the standard acronym for a list of vaccine components.
    
    Args:
        components: List of component names (e.g., ["DIPHTHERIA", "TETANUS", "PERTUSSIS"])
    
    Returns:
        Standard acronym (e.g., "DTP"), or None if no match
    """
    if not components:
        return None
    
    # Normalize and sort for order-independent matching
    normalized = sorted([c.upper().strip() for c in components])
    key = " + ".join(normalized)
    
    return VACCINE_COMPONENTS_TO_ACRONYM.get(key)


def match_vaccine_text(text: str) -> Tuple[Optional[str], Optional[List[str]]]:
    """
    Match vaccine text to acronym and components.
    
    Works bidirectionally:
    - If text contains an acronym, expands it to components
    - If text contains component names, finds the matching acronym
    
    Args:
        text: Vaccine description text
    
    Returns:
        (acronym, components) tuple, or (None, None) if no match
    """
    text_upper = text.upper()
    
    # Check if text starts with or contains a known acronym
    for acronym in sorted(VACCINE_ACRONYM_TO_COMPONENTS.keys(), key=len, reverse=True):
        # Check for acronym as standalone word
        import re
        if re.search(rf'\b{re.escape(acronym)}\b', text_upper):
            return acronym, VACCINE_ACRONYM_TO_COMPONENTS[acronym]
    
    # Extract components from text
    components = normalize_vaccine_components(text)
    if components:
        acronym = get_vaccine_acronym(components)
        return acronym, components
    
    return None, None


# ============================================================================
# REGIONAL CANONICAL NAMES - Map US names to PH/WHO names
# Philippines uses WHO naming conventions (e.g., PARACETAMOL not ACETAMINOPHEN)
# Applied when generating output to prefer regional names
# ============================================================================

REGIONAL_CANONICAL: Dict[str, str] = {
    # US name → PH/WHO name
    "ACETAMINOPHEN": "PARACETAMOL",
    "ALBUTEROL": "SALBUTAMOL",
    "EPINEPHRINE": "ADRENALINE",
    "NOREPINEPHRINE": "NORADRENALINE",
    "MEPERIDINE": "PETHIDINE",
}

# Reverse mapping: PH/WHO name → US name (for lookups)
REGIONAL_TO_US: Dict[str, str] = {v: k for k, v in REGIONAL_CANONICAL.items()}

def get_regional_canonical(name: str) -> str:
    """Get the regional (PH/WHO) canonical name for a drug."""
    return REGIONAL_CANONICAL.get(name.upper(), name.upper())

def get_us_canonical(name: str) -> str:
    """Get the US canonical name for a drug (for database lookups)."""
    return REGIONAL_TO_US.get(name.upper(), name.upper())


# ============================================================================
# MULTIWORD GENERICS - Drug names that contain spaces
# These should be preserved as single tokens during tokenization
# ============================================================================

MULTIWORD_GENERICS: Set[str] = {
    # Acids
    "TRANEXAMIC ACID", "FOLIC ACID", "ASCORBIC ACID", "VALPROIC ACID",
    "ACETYLSALICYLIC ACID", "HYALURONIC ACID", "RETINOIC ACID",
    "CLAVULANIC ACID", "FUSIDIC ACID", "NALIDIXIC ACID", "MEFENAMIC ACID",
    "URSODEOXYCHOLIC ACID", "CHENODEOXYCHOLIC ACID", "AMINO ACID",
    # Nitrates
    "ISOSORBIDE MONONITRATE", "ISOSORBIDE DINITRATE",
    "GLYCERYL TRINITRATE",
    # Salts (the compound IS the drug)
    "SODIUM CHLORIDE", "POTASSIUM CHLORIDE", "CALCIUM CHLORIDE",
    "MAGNESIUM SULFATE", "FERROUS SULFATE", "ZINC SULFATE",
    "SODIUM BICARBONATE", "CALCIUM CARBONATE", "MAGNESIUM HYDROXIDE",
    "ALUMINUM HYDROXIDE", "FERROUS FUMARATE", "CALCIUM GLUCONATE",
    # IV solution bases
    "LACTATED RINGER'S", "LACTATED RINGER'S SOLUTION",
    "ACETATED RINGER'S", "ACETATED RINGER'S SOLUTION",
    "RINGER'S SOLUTION", "WATER FOR INJECTION",
    # Insulins
    "INSULIN GLARGINE", "INSULIN LISPRO", "INSULIN ASPART",
    "INSULIN DETEMIR", "INSULIN DEGLUDEC", "INSULIN GLULISINE",
    "INSULIN HUMAN", "INSULIN REGULAR",
    # Biologics
    "HUMAN ALBUMIN", "ALBUMIN HUMAN",
    "TETANUS ANTITOXIN", "TETANUS IMMUNOGLOBULIN",
    "HEPATITIS B IMMUNOGLOBULIN",
    # Others
    "VITAMIN A", "VITAMIN B", "VITAMIN B1", "VITAMIN B2", "VITAMIN B6",
    "VITAMIN B12", "VITAMIN C", "VITAMIN D", "VITAMIN D3", "VITAMIN E",
    "VITAMIN K", "VITAMIN K1",
    "FERRIC CARBOXYMALTOSE", "IRON SUCROSE", "IRON DEXTRAN",
    "SODIUM HYALURONATE", "CHONDROITIN SULFATE",
    # Vaccines - canonical names (must match VACCINE_CANONICAL keys)
    "BCG VACCINE", "HEPATITIS A VACCINE", "HEPATITIS B VACCINE",
    "HEPATITIS A + B VACCINE", "DTP VACCINE", "DT VACCINE",
    "DTP + HIB VACCINE", "DTP + HEPATITIS B VACCINE",
    "DTP + IPV VACCINE", "DTP + IPV + HIB VACCINE",
    "IPV VACCINE", "OPV VACCINE", "MEASLES VACCINE", "MUMPS VACCINE",
    "RUBELLA VACCINE", "MMR VACCINE", "VARICELLA VACCINE",
    "PNEUMOCOCCAL VACCINE", "MENINGOCOCCAL VACCINE", "HIB VACCINE",
    "INFLUENZA VACCINE", "ROTAVIRUS VACCINE", "RABIES VACCINE",
    "YELLOW FEVER VACCINE", "HPV VACCINE", "TYPHOID VACCINE",
    "JAPANESE ENCEPHALITIS VACCINE", "PENTAVALENT VACCINE",
}


# ============================================================================
# TEXT NORMALIZATION UTILITIES
# These are included here so submodules only need to import unified_constants.py
# ============================================================================

import re as _re
import unicodedata as _unicodedata

# Pre-compiled patterns for normalize_text
_GM_TOKEN_RX = _re.compile(r"(?<![a-z])gms?(?![a-z])")

# Lowercase form words for matching (sorted by length, longest first)
_FORM_WORDS_LOWER = sorted(
    {k.lower() for k in FORM_TO_ROUTE.keys()},
    key=len, reverse=True
)

def normalize_text(s: str) -> str:
    """
    Produce the canonical normalized text used for matching and parsing.
    
    This is the standard text normalization function used across the pipeline.
    Converts to lowercase, normalizes unicode, expands abbreviations, etc.
    """
    if not isinstance(s, str):
        return ""
    s = _unicodedata.normalize("NFKD", s)
    s = "".join(c for c in s if not _unicodedata.combining(c))
    s = s.lower()
    s = _re.sub(r"\biv\b", "intravenous", s)
    s = _re.sub(r"[^\w%/+\.\- ]+", " ", s)
    s = s.replace("microgram", "mcg").replace("μg", "mcg").replace("µg", "mcg")
    s = _re.sub(r"(?<![a-z])cc(?![a-z])", "ml", s).replace("milli litre", "ml").replace("milliliter", "ml")
    s = _GM_TOKEN_RX.sub("g", s)
    s = s.replace("milligram", "mg")
    s = s.replace("polymixin", "polymyxin")
    s = s.replace("hydrochlorde", "hydrochloride")
    s = _re.sub(r"\s+", " ", s).strip()
    return s

def parse_form_from_text(s_norm: str) -> str | None:
    """
    Extract a recognized dosage form keyword from normalized text.
    
    Args:
        s_norm: Normalized text (should be output of normalize_text())
    
    Returns:
        The first matching form keyword, or None if not found.
    """
    for fw in _FORM_WORDS_LOWER:
        if _re.search(rf"\b{_re.escape(fw)}\b", s_norm):
            return fw
    return None


# ============================================================================
# CANONICAL GENERICS - DrugBank-preferred names with DrugBank IDs
# These are combination drugs or special products that need explicit DrugBank IDs
# Only use when DrugBank doesn't already have the combination
# ============================================================================

CANONICAL_GENERICS: List[Dict[str, Any]] = [
    # Combination antibiotics
    {"drugbank_id": "DB00766", "generic_name": "AMOXICILLIN + CLAVULANIC ACID", "source": "canonical"},
    {"drugbank_id": None, "generic_name": "COTRIMOXAZOLE", "source": "canonical"},
    {"drugbank_id": None, "generic_name": "SULFAMETHOXAZOLE + TRIMETHOPRIM", "source": "canonical"},
    {"drugbank_id": None, "generic_name": "AMPICILLIN + SULBACTAM", "source": "canonical"},
    {"drugbank_id": None, "generic_name": "SULTAMICILLIN", "source": "canonical"},
    # TB drug combinations (first-line anti-tuberculosis)
    {"drugbank_id": None, "generic_name": "ISONIAZID + RIFAMPICIN + PYRAZINAMIDE + ETHAMBUTOL", "source": "canonical"},  # HRZE
    {"drugbank_id": None, "generic_name": "ISONIAZID + RIFAMPICIN + PYRAZINAMIDE", "source": "canonical"},  # HRZ
    {"drugbank_id": None, "generic_name": "ISONIAZID + RIFAMPICIN + ETHAMBUTOL", "source": "canonical"},  # HRE
    {"drugbank_id": None, "generic_name": "ISONIAZID + RIFAMPICIN", "source": "canonical"},
    {"drugbank_id": None, "generic_name": "ISONIAZID + ETHAMBUTOL", "source": "canonical"},
    # IV Solutions
    {"drugbank_id": None, "generic_name": "ACETATED RINGER'S SOLUTION", "source": "canonical"},
    {"drugbank_id": None, "generic_name": "RINGER'S SOLUTION, ACETATED", "source": "canonical"},
    {"drugbank_id": None, "generic_name": "LACTATED RINGER'S SOLUTION", "source": "canonical"},
    {"drugbank_id": None, "generic_name": "RINGER'S SOLUTION, LACTATED", "source": "canonical"},
    # Parenteral nutrition
    {"drugbank_id": None, "generic_name": "AMINO ACIDS", "source": "canonical"},
    {"drugbank_id": None, "generic_name": "AMINO ACID SOLUTIONS FOR HEPATIC FAILURE", "source": "canonical"},
    {"drugbank_id": None, "generic_name": "AMINO ACID SOLUTIONS FOR INFANTS", "source": "canonical"},
    {"drugbank_id": None, "generic_name": "AMINO ACID SOLUTIONS FOR RENAL CONDITIONS", "source": "canonical"},
    {"drugbank_id": None, "generic_name": "AMINO ACID SOLUTIONS FOR IMMUNONUTRITION/IMMUNOENHANCEMENT", "source": "canonical"},
    # Antacid combinations - DB06723 is Maalox (aluminum hydroxide + magnesium hydroxide)
    {"drugbank_id": "DB06723", "generic_name": "ALUMINIUM HYDROXIDE + MAGNESIUM HYDROXIDE", "source": "canonical"},
    {"drugbank_id": "DB06723", "generic_name": "ALUMINUM HYDROXIDE + MAGNESIUM HYDROXIDE", "source": "canonical"},
    # Vitamin combinations
    {"drugbank_id": None, "generic_name": "ALENDRONATE + CHOLECALCIFEROL", "source": "canonical"},
    {"drugbank_id": None, "generic_name": "VITAMINS INTRAVENOUS, FAT-SOLUBLE", "source": "canonical"},
    {"drugbank_id": None, "generic_name": "VITAMINS INTRAVENOUS, WATER-SOLUBLE", "source": "canonical"},
    # Philippine herbal
    {"drugbank_id": None, "generic_name": "AKAPULKO", "source": "canonical"},
    {"drugbank_id": None, "generic_name": "CASSIA ALATA", "source": "canonical"},
    # Topical combinations - DB00443 is Enstilar (calcipotriol + betamethasone)
    {"drugbank_id": None, "generic_name": "BENZOIC ACID + SALICYLIC ACID", "source": "canonical"},
    {"drugbank_id": None, "generic_name": "BACITRACIN + NEOMYCIN + POLYMYXIN B", "source": "canonical"},
    {"drugbank_id": "DB00443", "generic_name": "CALCIPOTRIOL + BETAMETHASONE", "source": "canonical"},
    # Respiratory combinations - DB01222 is Budesonide/Formoterol
    {"drugbank_id": "DB01222", "generic_name": "BUDESONIDE + FORMOTEROL", "source": "canonical"},
    {"drugbank_id": None, "generic_name": "SALMETEROL + FLUTICASONE", "source": "canonical"},
    # Calcium + Vitamin D combinations
    {"drugbank_id": None, "generic_name": "CALCIUM CARBONATE + CHOLECALCIFEROL", "source": "canonical"},
    {"drugbank_id": None, "generic_name": "CALCIUM + VITAMIN D", "source": "canonical"},
    # Dialysis solutions
    {"drugbank_id": None, "generic_name": "DIALYSATE", "source": "canonical"},
    {"drugbank_id": None, "generic_name": "DIALYSATE (ACETATE BASED)", "source": "canonical"},
    {"drugbank_id": None, "generic_name": "DIALYSATE (BICARBONATE BASED)", "source": "canonical"},
    {"drugbank_id": None, "generic_name": "DOUSOL", "source": "canonical"},
    # Insulin - DB00030 is insulin human
    {"drugbank_id": "DB00030", "generic_name": "BIPHASIC ISOPHANE HUMAN INSULIN", "source": "canonical"},
    {"drugbank_id": "DB00030", "generic_name": "INSULIN HUMAN", "source": "canonical"},
    # Vaccines (CDC standard abbreviations) - using DrugBank mixture IDs where available
    {"drugbank_id": None, "generic_name": "DIPHTHERIA-TETANUS TOXOIDS", "source": "canonical"},
    {"drugbank_id": None, "generic_name": "DT", "source": "canonical"},  # Diphtheria-Tetanus (pediatric)
    {"drugbank_id": None, "generic_name": "TD", "source": "canonical"},  # Tetanus-Diphtheria (adult)
    {"drugbank_id": None, "generic_name": "DTP", "source": "canonical"},  # Diphtheria-Tetanus-Pertussis
    {"drugbank_id": None, "generic_name": "DTAP", "source": "canonical"},  # Diphtheria-Tetanus-acellular Pertussis
    {"drugbank_id": "DB10584", "generic_name": "DTP + HIB", "source": "canonical"},  # Infanrix-Hexa
    {"drugbank_id": "DB10584", "generic_name": "DTAP + HIB", "source": "canonical"},
    {"drugbank_id": "DB10797", "generic_name": "DTP + HEPATITIS B VACCINE", "source": "canonical"},  # Pediarix
    {"drugbank_id": "DB11627", "generic_name": "HEPATITIS A + B VACCINE", "source": "canonical"},  # Twinrix
    {"drugbank_id": None, "generic_name": "HEPATITIS B VACCINE", "source": "canonical"},
    {"drugbank_id": None, "generic_name": "HIB VACCINE", "source": "canonical"},
    {"drugbank_id": None, "generic_name": "PNEUMOCOCCAL VACCINE", "source": "canonical"},
    {"drugbank_id": None, "generic_name": "MMR", "source": "canonical"},  # Measles-Mumps-Rubella
    {"drugbank_id": None, "generic_name": "MMRV", "source": "canonical"},  # Measles-Mumps-Rubella-Varicella
    {"drugbank_id": None, "generic_name": "IPV", "source": "canonical"},  # Inactivated Polio Vaccine
    {"drugbank_id": None, "generic_name": "OPV", "source": "canonical"},  # Oral Polio Vaccine
    {"drugbank_id": None, "generic_name": "BCG", "source": "canonical"},  # Bacillus Calmette-Guérin
    # Special non-drug category
    {"drugbank_id": None, "generic_name": "DRUGS AND MEDICINES NOT NEEDED DURING THIS PARTICULAR EPISODE OF CARE", "source": "canonical"},
    # Hormone combinations
    {"drugbank_id": None, "generic_name": "CONJUGATED EQUINE ESTROGEN + MEDROXYPROGESTERONE ACETATE", "source": "canonical"},
    # Cephalosporins (already in DrugBank but adding for robustness)
    {"drugbank_id": "DB00493", "generic_name": "CEFOTAXIME", "source": "canonical"},
    {"drugbank_id": "DB01212", "generic_name": "CEFTRIAXONE", "source": "canonical"},
    # Special compounds
    {"drugbank_id": None, "generic_name": "DIMERCAPTOPROPANE SULFONATE", "source": "canonical"},
    {"drugbank_id": None, "generic_name": "DMPS", "source": "canonical"},
]

# ATC codes for canonical generics (only for items not in DrugBank ATC table)
CANONICAL_ATC_MAPPINGS: List[Dict[str, Any]] = [
    # Amoxicillin + Clavulanic acid -> J01CR02
    {"drugbank_id": "DB00766", "generic_name": "AMOXICILLIN + CLAVULANIC ACID", "atc_code": "J01CR02"},
    {"drugbank_id": "DB00766", "generic_name": "CO-AMOXICLAV", "atc_code": "J01CR02"},
    # Sulfamethoxazole + Trimethoprim -> J01EE01
    {"drugbank_id": None, "generic_name": "COTRIMOXAZOLE", "atc_code": "J01EE01"},
    {"drugbank_id": None, "generic_name": "SULFAMETHOXAZOLE + TRIMETHOPRIM", "atc_code": "J01EE01"},
    # Ampicillin + Sulbactam
    {"drugbank_id": None, "generic_name": "AMPICILLIN + SULBACTAM", "atc_code": "J01CR01"},
    {"drugbank_id": None, "generic_name": "SULTAMICILLIN", "atc_code": "J01CR04"},
    # TB drug combinations (WHO ATC J04AM - combinations of anti-mycobacterials)
    # Only use 3+ letter abbreviations to avoid false matches
    {"drugbank_id": None, "generic_name": "ISONIAZID + RIFAMPICIN + PYRAZINAMIDE + ETHAMBUTOL", "atc_code": "J04AM06"},  # HRZE
    {"drugbank_id": None, "generic_name": "HRZE", "atc_code": "J04AM06"},
    {"drugbank_id": None, "generic_name": "ISONIAZID + RIFAMPICIN + PYRAZINAMIDE", "atc_code": "J04AM05"},  # HRZ
    {"drugbank_id": None, "generic_name": "HRZ", "atc_code": "J04AM05"},
    {"drugbank_id": None, "generic_name": "ISONIAZID + RIFAMPICIN + ETHAMBUTOL", "atc_code": "J04AM06"},  # HRE (3-drug, same ATC as HRZE)
    {"drugbank_id": None, "generic_name": "HRE", "atc_code": "J04AM06"},
    {"drugbank_id": None, "generic_name": "ISONIAZID + RIFAMPICIN", "atc_code": "J04AM02"},
    {"drugbank_id": None, "generic_name": "ISONIAZID + ETHAMBUTOL", "atc_code": "J04AM03"},
    # IV Solutions - Ringer's
    {"drugbank_id": None, "generic_name": "ACETATED RINGER'S SOLUTION", "atc_code": "B05BB01"},
    {"drugbank_id": None, "generic_name": "RINGER'S SOLUTION, ACETATED", "atc_code": "B05BB01"},
    {"drugbank_id": None, "generic_name": "LACTATED RINGER'S SOLUTION", "atc_code": "B05BB01"},
    {"drugbank_id": None, "generic_name": "RINGER'S SOLUTION, LACTATED", "atc_code": "B05BB01"},
    # Parenteral nutrition - Amino acids
    {"drugbank_id": None, "generic_name": "AMINO ACIDS", "atc_code": "B05BA01"},
    {"drugbank_id": None, "generic_name": "AMINO ACID SOLUTIONS FOR HEPATIC FAILURE", "atc_code": "B05BA01"},
    {"drugbank_id": None, "generic_name": "AMINO ACID SOLUTIONS FOR INFANTS", "atc_code": "B05BA01"},
    {"drugbank_id": None, "generic_name": "AMINO ACID SOLUTIONS FOR RENAL CONDITIONS", "atc_code": "B05BA01"},
    {"drugbank_id": None, "generic_name": "AMINO ACID SOLUTIONS FOR IMMUNONUTRITION/IMMUNOENHANCEMENT", "atc_code": "B05BA01"},
    # Antacids
    {"drugbank_id": "DB06723", "generic_name": "ALUMINIUM HYDROXIDE + MAGNESIUM HYDROXIDE", "atc_code": "A02AD01"},
    {"drugbank_id": "DB06723", "generic_name": "ALUMINUM HYDROXIDE + MAGNESIUM HYDROXIDE", "atc_code": "A02AD01"},
    # Bone drugs
    {"drugbank_id": None, "generic_name": "ALENDRONATE + CHOLECALCIFEROL", "atc_code": "M05BB03"},
    # Vitamins IV
    {"drugbank_id": None, "generic_name": "VITAMINS INTRAVENOUS, FAT-SOLUBLE", "atc_code": "B05XC"},
    {"drugbank_id": None, "generic_name": "VITAMINS INTRAVENOUS, WATER-SOLUBLE", "atc_code": "B05XC"},
    # Philippine herbal
    {"drugbank_id": None, "generic_name": "AKAPULKO", "atc_code": "D01AE"},
    {"drugbank_id": None, "generic_name": "CASSIA ALATA", "atc_code": "D01AE"},
    # Topical combinations
    {"drugbank_id": None, "generic_name": "BENZOIC ACID + SALICYLIC ACID", "atc_code": "D01AE20"},
    {"drugbank_id": None, "generic_name": "BACITRACIN + NEOMYCIN + POLYMYXIN B", "atc_code": "D06AX"},
    {"drugbank_id": "DB00443", "generic_name": "CALCIPOTRIOL + BETAMETHASONE", "atc_code": "D05AX52"},
    # Respiratory combinations
    {"drugbank_id": "DB01222", "generic_name": "BUDESONIDE + FORMOTEROL", "atc_code": "R03AK07"},
    {"drugbank_id": None, "generic_name": "SALMETEROL + FLUTICASONE", "atc_code": "R03AK06"},
    # Calcium + Vitamin D
    {"drugbank_id": None, "generic_name": "CALCIUM CARBONATE + CHOLECALCIFEROL", "atc_code": "A12AX"},
    {"drugbank_id": None, "generic_name": "CALCIUM + VITAMIN D", "atc_code": "A12AX"},
    # Dialysis solutions
    {"drugbank_id": None, "generic_name": "DIALYSATE", "atc_code": "B05DB"},
    {"drugbank_id": None, "generic_name": "DIALYSATE (ACETATE BASED)", "atc_code": "B05DB"},
    {"drugbank_id": None, "generic_name": "DIALYSATE (BICARBONATE BASED)", "atc_code": "B05DB"},
    {"drugbank_id": None, "generic_name": "DOUSOL", "atc_code": "B05DB"},
    # Insulin
    {"drugbank_id": "DB00030", "generic_name": "BIPHASIC ISOPHANE HUMAN INSULIN", "atc_code": "A10AD01"},
    {"drugbank_id": "DB00030", "generic_name": "INSULIN HUMAN", "atc_code": "A10AB01"},
    # Vaccines
    {"drugbank_id": None, "generic_name": "DIPHTHERIA-TETANUS TOXOIDS", "atc_code": "J07AM51"},
    {"drugbank_id": None, "generic_name": "DTP", "atc_code": "J07CA01"},
    {"drugbank_id": None, "generic_name": "DTAP", "atc_code": "J07CA01"},
    {"drugbank_id": "DB10584", "generic_name": "DTP + HIB", "atc_code": "J07CA06"},
    {"drugbank_id": "DB10584", "generic_name": "DTAP + HIB", "atc_code": "J07CA06"},
    {"drugbank_id": "DB10797", "generic_name": "DTP + HEPATITIS B VACCINE", "atc_code": "J07CA05"},
    {"drugbank_id": "DB11627", "generic_name": "HEPATITIS A + B VACCINE", "atc_code": "J07BC20"},
    {"drugbank_id": None, "generic_name": "HEPATITIS B VACCINE", "atc_code": "J07BC01"},
    {"drugbank_id": None, "generic_name": "HIB VACCINE", "atc_code": "J07AG01"},
    {"drugbank_id": None, "generic_name": "PNEUMOCOCCAL VACCINE", "atc_code": "J07AL01"},
    {"drugbank_id": None, "generic_name": "MMR", "atc_code": "J07BD52"},
    {"drugbank_id": None, "generic_name": "MMRV", "atc_code": "J07BD54"},
    {"drugbank_id": None, "generic_name": "IPV", "atc_code": "J07BF03"},
    {"drugbank_id": None, "generic_name": "OPV", "atc_code": "J07BF02"},
    {"drugbank_id": None, "generic_name": "BCG", "atc_code": "J07AN01"},
    # Hormone combinations
    {"drugbank_id": None, "generic_name": "CONJUGATED EQUINE ESTROGEN + MEDROXYPROGESTERONE ACETATE", "atc_code": "G03FB"},
    # Cephalosporins
    {"drugbank_id": "DB00493", "generic_name": "CEFOTAXIME", "atc_code": "J01DD01"},
    {"drugbank_id": "DB01212", "generic_name": "CEFTRIAXONE", "atc_code": "J01DD04"},
    # Special compounds
    {"drugbank_id": None, "generic_name": "DIMERCAPTOPROPANE SULFONATE", "atc_code": "V03AB"},
    {"drugbank_id": None, "generic_name": "DMPS", "atc_code": "V03AB"},
]


# ============================================================================
# EXPORTS
# ============================================================================

__all__ = [
    # Categories
    "CATEGORY_GENERIC", "CATEGORY_SALT", "CATEGORY_DOSE",
    "CATEGORY_FORM", "CATEGORY_ROUTE", "CATEGORY_OTHER",
    
    # Rules
    "GENERIC_MATCH_REQUIRED", "SALT_FLEXIBLE", "DOSE_FLEXIBLE",
    
    # Token sets
    "STOPWORDS", "STOPWORDS_LOWER",
    "SALT_TOKENS", "SALT_TOKENS_LOWER",
    "PURE_SALT_COMPOUNDS",
    "SALT_CATIONS", "SALT_ANIONS", "ANION_TO_CATIONS",
    "ELEMENT_DRUGS",
    "UNIT_TOKENS", "UNIT_TOKENS_LOWER",
    "CONNECTIVE_WORDS", "SALT_TAIL_BREAK_TOKENS",
    
    # Mappings
    "FORM_CANON", "ROUTE_CANON", "FORM_TO_ROUTE", "FORM_TO_ROUTES",
    "FORM_EQUIVALENCE_GROUPS", "FORM_EQUIVALENTS",
    "WEIGHT_UNIT_FACTORS",
    
    # ATC patterns
    "ATC_COMBINATION_PATTERNS", "COMBINATION_ATC_SUFFIXES",
    
    # Synonyms, multiword generics, regional names
    "GARBAGE_TOKENS",
    "GENERIC_SYNONYMS", "IV_FLUID_SYNONYMS", "DRUGBANK_COMPONENT_SYNONYMS",
    "ALL_DRUG_SYNONYMS",
    "SPELLING_SYNONYMS", "MULTIWORD_GENERICS",
    "REGIONAL_CANONICAL", "REGIONAL_TO_US",
    "get_regional_canonical", "get_us_canonical",
    
    # Vaccine normalization
    "VACCINE_CANONICAL", "normalize_vaccine_name",
    # Vaccine acronym bidirectional lookup (WHO/CDC standard)
    "VACCINE_ACRONYM_TO_COMPONENTS", "VACCINE_COMPONENT_KEYWORDS",
    "VACCINE_COMPONENTS_TO_ACRONYM",
    "normalize_vaccine_components", "expand_vaccine_acronym",
    "get_vaccine_acronym", "match_vaccine_text",
    
    # Helper functions
    "get_canonical_form", "get_canonical_route",
    "is_stopword", "is_salt_token", "is_pure_salt_compound",
    "is_element_drug", "is_unit_token", "is_combination_atc",
    "forms_are_equivalent", "infer_route_from_form",
    "get_valid_routes_for_form", "is_valid_form_route_pair",
    "parse_compound_salt", "get_related_salts",
    
    # Text utilities (for submodules)
    "normalize_text", "parse_form_from_text",
    
    # Canonical generics and ATC mappings
    "CANONICAL_GENERICS", "CANONICAL_ATC_MAPPINGS",
]
