"""
Unified Constants for Drug Pipeline

This file consolidates all hardcoded token lists from across the pipeline.
All scripts should import constants from here to ensure consistency.

Created: Nov 28, 2025 (Phase 1, TODO #22)
Sources: constants.py, combos_drugs.py, text_utils_drugs.py, routes_forms_drugs.py,
         generic_normalization.py, scoring.py, and debug/old_files/*.py
"""

from typing import Dict, List, Optional, Set, Tuple

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
    # Tablets
    "TAB": "TABLET", "TABS": "TABLET", "TABLET": "TABLET", "TABLETS": "TABLET",
    
    # Capsules
    "CAP": "CAPSULE", "CAPS": "CAPSULE", "CAPSULE": "CAPSULE", "CAPSULES": "CAPSULE",
    "CAPLET": "CAPSULE",
    
    # Bottles (container, not form - but often used as form)
    "BOT": "BOTTLE", "BOTT": "BOTTLE", "BOTTLE": "BOTTLE", "BOTTLES": "BOTTLE",
    
    # Vials
    "VIAL": "VIAL", "VIALS": "VIAL", "VL": "VIAL",
    
    # Injections
    "INJ": "INJECTION", "INJECTABLE": "INJECTION", "INJECTION": "INJECTION",
    
    # Syrups
    "SYR": "SYRUP", "SYRUP": "SYRUP", "SYRUPS": "SYRUP",
    
    # Topical forms
    "LOTION": "LOTION", "CREAM": "CREAM", "GEL": "GEL", "OINTMENT": "OINTMENT",
    "PASTE": "PASTE", "FOAM": "FOAM", "EMULSION": "EMULSION", "SHAMPOO": "SHAMPOO",
    "SOAP": "SOAP", "WASH": "WASH",
    
    # Solutions
    "SOLUTION": "SOLUTION", "SOLN": "SOLUTION", "SOL": "SOLUTION",
    "SOLUTIONS": "SOLUTION",
    
    # Suspensions
    "SUSPENSION": "SUSPENSION", "SUSP": "SUSPENSION", "SUSPENSIONS": "SUSPENSION",
    
    # Drops
    "DROPS": "DROPS", "DROP": "DROPS", "GTTS": "DROPS",
    "EYE DROPS": "DROPS", "EAR DROPS": "DROPS", "NASAL DROPS": "DROPS",
    
    # Powders
    "POWDER": "POWDER", "PWDR": "POWDER", "POWDERS": "POWDER",
    
    # Granules
    "GRANULES": "GRANULES", "GRANULE": "GRANULES", "GRAN": "GRANULES",
    
    # Sachets
    "SACHET": "SACHET", "SACHETS": "SACHET",
    
    # Inhalers
    "INHALER": "INHALER", "INH": "INHALER", "INHALERS": "INHALER",
    "NEBULE": "NEBULE", "NEBULES": "NEBULE", "NEB": "NEBULE", "NEBS": "NEBULE",
    "DPI": "DPI", "MDI": "MDI", "AEROSOL": "AEROSOL",
    "METERED DOSE INHALER": "MDI", "DRY POWDER INHALER": "DPI",
    "INHAL.AEROSOL": "AEROSOL", "INHAL.POWDER": "DPI",
    "INHAL.SOLUTION": "NEBULE", "ORAL AEROSOL": "AEROSOL",
    
    # Sprays
    "SPRAY": "SPRAY", "SPRAYS": "SPRAY", "NASAL SPRAY": "SPRAY",
    "DISPENSER": "SPRAY",
    
    # Suppositories
    "SUPPOSITORY": "SUPPOSITORY", "SUPP": "SUPPOSITORY", "SUPPOSITORIES": "SUPPOSITORY",
    
    # Patches
    "PATCH": "PATCH", "PATCHES": "PATCH",
    
    # Films
    "FILM": "FILM", "FILMS": "FILM",
    
    # Ampules
    "AMPULE": "AMPULE", "AMPUL": "AMPULE", "AMP": "AMPULE",
    "AMPOULE": "AMPULE", "AMPS": "AMPULE", "AMPULES": "AMPULE",
    
    # Other forms
    "LOZENGE": "LOZENGE", "LOZENGES": "LOZENGE",
    "MOUTHWASH": "MOUTHWASH",
    "CHEWING GUM": "CHEWING GUM",
    "ELIXIR": "ELIXIR",
    "IMPLANT": "IMPLANT", "IMPLANTS": "IMPLANT",
    "OVULE": "OVULE", "OVULES": "OVULE",
    "ENEMA": "ENEMA", "ENEMAS": "ENEMA",
    "BAR": "BAR", "BARS": "BAR",
    "STRIP": "STRIP", "STRIPS": "STRIP",
    "WAFER": "WAFER", "WAFERS": "WAFER",
    "PELLET": "PELLET", "PELLETS": "PELLET",
    "RING": "RING", "RINGS": "RING",
    "INSERT": "INSERT", "INSERTS": "INSERT",
    "SWAB": "SWAB", "SWABS": "SWAB",
    "CLOTH": "CLOTH", "CLOTHS": "CLOTH",
    "SPONGE": "SPONGE", "SPONGES": "SPONGE",
    "DRESSING": "DRESSING", "DRESSINGS": "DRESSING",
    "LAMELLA": "LAMELLA",
    "INSTILL.SOLUTION": "DROPS",
    "S.C. IMPLANT": "IMPLANT",
}

# ============================================================================
# ROUTE CANONICALIZATION
# Maps various route spellings/abbreviations to canonical routes
# ============================================================================

ROUTE_CANON: Dict[str, str] = {
    # Oral
    "PO": "ORAL", "OR": "ORAL", "ORAL": "ORAL",
    "PER OREM": "ORAL", "BY MOUTH": "ORAL", "PER OS": "ORAL",
    
    # Intravenous
    "IV": "INTRAVENOUS", "INTRAVENOUS": "INTRAVENOUS",
    
    # Intramuscular
    "IM": "INTRAMUSCULAR", "INTRAMUSCULAR": "INTRAMUSCULAR",
    
    # Subcutaneous
    "SC": "SUBCUTANEOUS", "SUBCUTANEOUS": "SUBCUTANEOUS",
    "SUBCUT": "SUBCUTANEOUS", "SQ": "SUBCUTANEOUS",
    "SUBDERMAL": "SUBCUTANEOUS",
    
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
    "NASAL": "NASAL", "INTRANASAL": "NASAL", "PER NASAL": "NASAL",
    
    # Inhalation
    "INH": "INHALATION", "INHALATION": "INHALATION",
    "NEB": "INHALATION", "INHALER": "INHALATION",
    
    # Rectal
    "RECTAL": "RECTAL", "PR": "RECTAL", "PER RECTUM": "RECTAL",
    
    # Vaginal
    "VAGINAL": "VAGINAL", "PV": "VAGINAL", "PER VAGINAM": "VAGINAL",
    
    # Other
    "INTRATHECAL": "INTRATHECAL",
    "INTRADERMAL": "INTRADERMAL", "ID": "INTRADERMAL",
}

# ============================================================================
# FORM TO ROUTE MAPPING
# Default route inference when only form is known
# ============================================================================

FORM_TO_ROUTE: Dict[str, str] = {
    # Oral forms
    "TABLET": "ORAL", "TAB": "ORAL", "TABS": "ORAL",
    "CAPSULE": "ORAL", "CAP": "ORAL", "CAPS": "ORAL", "CAPLET": "ORAL",
    "SYRUP": "ORAL", "SYRUPS": "ORAL", "SYR": "ORAL",
    "SUSPENSION": "ORAL", "SUSPENSIONS": "ORAL", "SUSP": "ORAL",
    "SOLUTION": "ORAL", "SOLUTIONS": "ORAL", "SOLN": "ORAL",
    "SACHET": "ORAL", "GRANULE": "ORAL", "GRANULES": "ORAL",
    "LOZENGE": "ORAL", "MOUTHWASH": "ORAL", "ELIXIR": "ORAL",
    "DROPS": "ORAL", "ORAL DROPS": "ORAL",
    "CHEWING GUM": "ORAL",
    
    # Topical forms
    "CREAM": "TOPICAL", "OINTMENT": "TOPICAL", "GEL": "TOPICAL",
    "LOTION": "TOPICAL", "SOAP": "TOPICAL", "SHAMPOO": "TOPICAL",
    "WASH": "TOPICAL", "PASTE": "TOPICAL", "FOAM": "TOPICAL",
    
    # Transdermal
    "PATCH": "TRANSDERMAL",
    
    # Inhalation forms
    "INHALER": "INHALATION", "NEBULE": "INHALATION", "NEB": "INHALATION",
    "MDI": "INHALATION", "DPI": "INHALATION", "AEROSOL": "INHALATION",
    "METERED DOSE INHALER": "INHALATION", "DRY POWDER INHALER": "INHALATION",
    "INHAL.AEROSOL": "INHALATION", "INHAL.POWDER": "INHALATION",
    "INHAL.SOLUTION": "INHALATION", "ORAL AEROSOL": "INHALATION",
    
    # Injectable forms (default to IV, may need context)
    "AMPULE": "INTRAVENOUS", "AMP": "INTRAVENOUS", "AMPUL": "INTRAVENOUS",
    "AMPOULE": "INTRAVENOUS", "VIAL": "INTRAVENOUS", "VL": "INTRAVENOUS",
    "INJ": "INTRAVENOUS", "INJECTION": "INTRAVENOUS",
    
    # Ophthalmic
    "DROP": "OPHTHALMIC", "EYE DROP": "OPHTHALMIC", "EYE DROPS": "OPHTHALMIC",
    "INSTILL.SOLUTION": "OPHTHALMIC", "LAMELLA": "OPHTHALMIC",
    
    # Otic
    "EAR DROP": "OTIC", "EAR DROPS": "OTIC",
    
    # Nasal
    "SPRAY": "NASAL", "NASAL SPRAY": "NASAL", "NASAL DROPS": "NASAL",
    
    # Rectal
    "SUPPOSITORY": "RECTAL", "SUPP": "RECTAL", "ENEMA": "RECTAL",
    
    # Vaginal
    "OVULE": "VAGINAL", "OVULES": "VAGINAL",
    
    # Subcutaneous
    "IMPLANT": "SUBCUTANEOUS", "S.C. IMPLANT": "SUBCUTANEOUS",
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
    """Infer route from form if possible."""
    canonical = get_canonical_form(form)
    return FORM_TO_ROUTE.get(canonical) or FORM_TO_ROUTE.get(form.lower())


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
    # Combination drug trade names -> generic components
    "CO-AMOXICLAV": "AMOXICILLIN + CLAVULANIC ACID",
    "COTRIMOXAZOLE": "SULFAMETHOXAZOLE + TRIMETHOPRIM",
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
}

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
    "FORM_CANON", "ROUTE_CANON", "FORM_TO_ROUTE",
    "FORM_EQUIVALENCE_GROUPS", "FORM_EQUIVALENTS",
    "WEIGHT_UNIT_FACTORS",
    
    # ATC patterns
    "ATC_COMBINATION_PATTERNS", "COMBINATION_ATC_SUFFIXES",
    
    # Synonyms, multiword generics, regional names
    "SPELLING_SYNONYMS", "MULTIWORD_GENERICS",
    "REGIONAL_CANONICAL", "REGIONAL_TO_US",
    "get_regional_canonical", "get_us_canonical",
    
    # Helper functions
    "get_canonical_form", "get_canonical_route",
    "is_stopword", "is_salt_token", "is_pure_salt_compound",
    "is_element_drug", "is_unit_token", "is_combination_atc",
    "forms_are_equivalent", "infer_route_from_form",
    "parse_compound_salt", "get_related_salts",
    
    # Text utilities (for submodules)
    "normalize_text", "parse_form_from_text",
]
