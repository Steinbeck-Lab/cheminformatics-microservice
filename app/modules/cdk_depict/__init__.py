"""CDK Enhancement Modules for Chemical Depiction.

This package provides advanced CDK-based functionality for molecular depiction:

Core Functionality:
- Hydrogen Display: Control hydrogen atom visibility (minimal, explicit, stereo, smart)
- CXSMILES: Extended SMILES parsing with highlighting and coordinates
- Abbreviations: Chemical group and reagent abbreviation system
- Dative Bonds: Coordinate bond perception and display
- Multicenter Bonds: η-complex and π-bonding support

Advanced Features:
- Annotations: Atom/bond numbering, mapping, CIP labels
- Aromatic Display: Circle-in-ring (donut) representation
- Style Presets: Color schemes for different backgrounds
- Advanced Controls: Zoom, stroke ratio, SVG units, structure flipping
- Reaction Arrows: Different arrow types for reactions

Specialized:
- Radical Perception: Unpaired electron detection
- MDL HILITE: Highlighting from MDL/SDF files
"""

# Hydrogen Display
from app.modules.cdk_depict.hydrogen_display import (
    set_hydrogen_display,
    setHydrogenDisplay,
    HydrogenDisplayType,
)

# CXSMILES Parser
from app.modules.cdk_depict.cxsmiles_parser import (
    extract_cxsmiles_highlighting,
)

# Abbreviations
from app.modules.cdk_depict.abbreviations import (
    ChemicalAbbreviations,
    AbbreviationMode,
    AbbreviationOptions,
)

# Dative Bonds
from app.modules.cdk_depict.dative_bonds import (
    DativeBondPerception,
    DativeBondMode,
)

# Multicenter Bonds
from app.modules.cdk_depict.multicenter_bonds import (
    MulticenterBonds,
    MulticenterStyle,
)

# Annotations
from app.modules.cdk_depict.annotations import (
    AnnotationSystem,
    AnnotationMode,
    MAPPING_COLORS,
)

# Aromatic Display
from app.modules.cdk_depict.aromatic_display import (
    AromaticDisplaySystem,
)

# Style Presets
from app.modules.cdk_depict.style_presets import (
    StylePresetSystem,
    StylePreset,
    StyleConfiguration,
    STYLE_CONFIGURATIONS,
)

# Advanced Controls
from app.modules.cdk_depict.advanced_controls import (
    AdvancedControls,
    SVGUnits,
)

# Reaction Depiction
from app.modules.cdk_depict.reaction_depiction import (
    ReactionArrowSystem,
    ReactionArrowType,
    ARROW_ABBREV_MAP,
)

# Radical Perception
from app.modules.cdk_depict.radical_perception import (
    RadicalPerception,
)

# MDL HILITE
from app.modules.cdk_depict.mdl_hilite import (
    MDLHiliteParser,
)


__all__ = [
    # Hydrogen Display
    "set_hydrogen_display",
    "setHydrogenDisplay",
    "HydrogenDisplayType",
    # CXSMILES
    "CXSMILESParser",
    "CXSMILESData",
    "CXSMILESError",
    "StereoGroupType",
    "extract_cxsmiles_highlighting",
    # Abbreviations
    "ChemicalAbbreviations",
    "AbbreviationMode",
    "AbbreviationOptions",
    # Dative Bonds
    "DativeBondPerception",
    "DativeBondMode",
    # Multicenter Bonds
    "MulticenterBonds",
    "MulticenterStyle",
    # Annotations
    "AnnotationSystem",
    "AnnotationMode",
    "MAPPING_COLORS",
    # Aromatic Display
    "AromaticDisplaySystem",
    # Style Presets
    "StylePresetSystem",
    "StylePreset",
    "StyleConfiguration",
    "STYLE_CONFIGURATIONS",
    # Advanced Controls
    "AdvancedControls",
    "SVGUnits",
    # Reaction Depiction
    "ReactionArrowSystem",
    "ReactionArrowType",
    "ARROW_ABBREV_MAP",
    # Radical Perception
    "RadicalPerception",
    # MDL HILITE
    "MDLHiliteParser",
]
