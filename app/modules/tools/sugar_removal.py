from __future__ import annotations
from enum import Enum

import app.modules.toolkits.cdk_wrapper as cdk

SCOB_CLASS = cdk.JClass(cdk.cdk_base + ".silent.SilentChemObjectBuilder")

SDU_PATH_BASE = cdk.cdk_base + ".tools"

SDU_CLASS_NAME = "SugarDetectionUtility"

SMI_FLAVOR_PATH_AND_CLASS_NAME = cdk.cdk_base + ".smiles.SmiFlavor"

SMILES_GENERATOR_PATH_AND_CLASS_NAME = cdk.cdk_base + ".smiles.SmilesGenerator"


class preservation_modes_enum(Enum):
    """
    Enum for sugar removal preservation modes. This specifies under what circumstances to preserve or
    discard structures that get disconnected from the central core in the sugar removal process.
    """

    ALL = 1  # Preserve all disconnected structures (note: this might lead to no circular sugar moieties being detected, depending on the other settings)
    HEAVY_ATOM_COUNT = (
        2  # Remove disconnected structures that do not have enough heavy atoms
    )
    MOLECULAR_WEIGHT = 3  # Remove disconnected structures that do not have a sufficient molecular weight


def get_sugar_info(
    molecule: any,
    gly_bond: bool = False,
    oxygen_atoms: bool = True,
    oxygen_atoms_threshold: float = 0.5,
    linear_sugars_in_rings: bool = False,
    linear_sugars_min_size: int = 4,
    linear_sugars_max_size: int = 7,
    linear_acidic_sugars: bool = False,
    spiro_sugars: bool = False,
    keto_sugars: bool = False,
) -> tuple:
    """
    Analyzes a molecule represented by a CDK IAtomContainer object to determine whether it contains sugars.

    This function utilizes the Sugar Removal/Detection Utility to check for the presence of circular or linear sugars.

    Args:
        molecule (IAtomContainer): CDK molecule object.
        glyBond (bool): Whether to consider only circular sugars with glycosidic bonds in the analysis. Default is False.
        oxygenAtoms (bool): Whether to consider only circular sugars with a sufficient number of exocyclic oxygen atoms
                            in the analysis (see oxygenAtomsThreshold). Default is True.
        oxygenAtomsThreshold (float): A number giving the minimum attached exocyclic oxygen atoms to atom number in the
                                      ring ratio a circular sugar needs to have to be considered in the analysis. Default
                                      is 0.5 (a 6-membered ring needs at least 3 attached exocyclic oxygen atoms). Must
                                      be positive!
        linearSugarsInRings (bool): Whether to consider linear sugars in rings. Default is False.
        linearSugarsMinSize (int): Minimum size of linear sugars to consider. Default is 4. Must be positive and higher
                                   than or equal to 1 and also smaller than the linear sugar candidate maximum size.
        linearSugarsMaxSize (int): Maximum size of linear sugars to consider. Default is 7. Must be positive and higher
                                   than or equal to 1 and also higher than the linear sugar candidate minimum size.
        linearAcidicSugars (bool): Whether to consider linear acidic sugars. Default is False.
        spiroSugars (bool): Whether spiro rings (rings that share one atom with another cycle) should be included in the
                            circular sugar detection. Default is False.
        ketoSugars (bool): Whether circular sugars with keto groups should be detected. Default is False.

    Returns:
        tuple: A tuple containing two boolean values indicating whether the molecule has linear sugars
               and whether the molecule has circular sugars. If no sugars are found, both values will be False.

    Raises:
        ValueError: If one of the numeric parameters is not valid.
    """
    _sugar_detection_utility = cdk.JClass(SDU_PATH_BASE + "." + SDU_CLASS_NAME)(
        SCOB_CLASS.getInstance(),
    )

    _sugar_detection_utility.setDetectCircularSugarsOnlyWithOGlycosidicBondSetting(
        gly_bond
    )
    _sugar_detection_utility.setDetectCircularSugarsOnlyWithEnoughExocyclicOxygenAtomsSetting(
        oxygen_atoms
    )
    _is_valid = _sugar_detection_utility.setExocyclicOxygenAtomsToAtomsInRingRatioThresholdSetting(
        oxygen_atoms_threshold
    )
    if not _is_valid:
        raise ValueError("oxygenAtomsThreshold must be a positive number.")
    _sugar_detection_utility.setDetectLinearSugarsInRingsSetting(linear_sugars_in_rings)
    if (
        linear_sugars_min_size < 0
        or linear_sugars_max_size < 1
        or linear_sugars_min_size >= linear_sugars_max_size
    ):
        raise ValueError(
            "linearSugarsMinSize and linearSugarsMaxSize must be positive integers, "
            "with linearSugarsMinSize being smaller than linearSugarsMaxSize."
        )
    _sugar_detection_utility.setLinearSugarCandidateMinSizeSetting(
        linear_sugars_min_size
    )
    _sugar_detection_utility.setLinearSugarCandidateMaxSizeSetting(
        linear_sugars_max_size
    )
    _sugar_detection_utility.setDetectLinearAcidicSugarsSetting(linear_acidic_sugars)
    _sugar_detection_utility.setDetectSpiroRingsAsCircularSugarsSetting(spiro_sugars)
    _sugar_detection_utility.setDetectCircularSugarsWithKetoGroupsSetting(keto_sugars)

    _has_linear_sugars = _sugar_detection_utility.hasLinearSugars(molecule)
    _has_circular_sugars = _sugar_detection_utility.hasCircularSugars(molecule)
    return _has_linear_sugars, _has_circular_sugars


def remove_linear_sugars(
    molecule: any,
    only_terminal: bool = True,
    preservation_mode: preservation_modes_enum = preservation_modes_enum.HEAVY_ATOM_COUNT,
    preservation_threshold: int = 5,
    linear_sugars_in_rings: bool = False,
    linear_sugars_min_size: int = 4,
    linear_sugars_max_size: int = 7,
    linear_acidic_sugars: bool = False,
    mark_attach_points: bool = False,
) -> str:
    """
    Removes linear sugars from a given CDK IAtomContainer object using the Sugar Removal/Detection Utility.

    Args:
        molecule (IAtomContainer): CDK molecule object.
        only_terminal (bool): Whether to remove only terminal linear sugars. Default is True.
        preservation_mode (preservation_modes_enum): Mode to determine which disconnected structures to preserve.
                                                     Default is preservation_modes_enum.HEAVY_ATOM_COUNT.
        preservation_threshold (int): Threshold value for the selected preservation mode. Default is 5 (heavy atoms).
        linear_sugars_in_rings (bool): Whether to consider linear sugars in rings. Default is False.
        linear_sugars_min_size (int): Minimum size of linear sugars to consider. Default is 4. Must be positive and
                                      higher than or equal to 1 and also smaller than the linear sugar candidate
                                      maximum size.
        linear_sugars_max_size (int): Maximum size of linear sugars to consider. Default is 7. Must be positive and
                                      higher than or equal to 1 and also higher than the linear sugar candidate
                                      minimum size.
        linear_acidic_sugars (bool): Whether to consider linear acidic sugars. Default is False.
        mark_attach_points (bool): Whether to mark the attachment points of removed sugars with a dummy atom. Default
                                   is False.

    Returns:
        str: The SMILES string with linear sugars removed, or a message indicating that no linear sugars were
        found ("No Linear sugar found").

    Raises:
        ValueError: If one of the numeric parameters is not valid.
        Exception: If an error occurs during the output SMILES generation process.
    """
    _smi_flavor = cdk.JClass(SMI_FLAVOR_PATH_AND_CLASS_NAME)
    _smiles_generator = cdk.JClass(SMILES_GENERATOR_PATH_AND_CLASS_NAME)(
        _smi_flavor.Absolute,
    )
    _sugar_detection_utility = cdk.JClass(SDU_PATH_BASE + "." + SDU_CLASS_NAME)(
        SCOB_CLASS.getInstance(),
    )

    _sugar_detection_utility.setRemoveOnlyTerminalSugarsSetting(only_terminal)
    _sru_preservation_modes_enum = cdk.JClass(
        SDU_PATH_BASE + "." + "SugarRemovalUtility$PreservationMode"
    )
    if preservation_mode == preservation_modes_enum.ALL:
        _sugar_detection_utility.setPreservationModeSetting(
            _sru_preservation_modes_enum.ALL
        )
    elif preservation_mode == preservation_modes_enum.HEAVY_ATOM_COUNT:
        _sugar_detection_utility.setPreservationModeSetting(
            _sru_preservation_modes_enum.HEAVY_ATOM_COUNT
        )
    elif preservation_mode == preservation_modes_enum.MOLECULAR_WEIGHT:
        _sugar_detection_utility.setPreservationModeSetting(
            _sru_preservation_modes_enum.MOLECULAR_WEIGHT
        )
    else:
        raise ValueError("Invalid preservation_mode specified.")
    if preservation_threshold < 0:
        raise ValueError("preservation_threshold must be a non-negative integer.")
    _sugar_detection_utility.setPreservationModeThresholdSetting(preservation_threshold)
    _sugar_detection_utility.setDetectLinearSugarsInRingsSetting(linear_sugars_in_rings)
    if (
        linear_sugars_min_size < 0
        or linear_sugars_max_size < 1
        or linear_sugars_min_size >= linear_sugars_max_size
    ):
        raise ValueError(
            "linearSugarsMinSize and linearSugarsMaxSize must be positive integers, "
            "with linearSugarsMinSize being smaller than linearSugarsMaxSize."
        )
    _sugar_detection_utility.setLinearSugarCandidateMinSizeSetting(
        linear_sugars_min_size
    )
    _sugar_detection_utility.setLinearSugarCandidateMaxSizeSetting(
        linear_sugars_max_size
    )
    _sugar_detection_utility.setDetectLinearAcidicSugarsSetting(linear_acidic_sugars)

    _has_linear_sugars = _sugar_detection_utility.hasLinearSugars(molecule)

    if _has_linear_sugars:
        _aglycone_and_sugars = _sugar_detection_utility.copyAndExtractAglyconeAndSugars(
            molecule,
            False,
            True,
            mark_attach_points,
        )
        if not _aglycone_and_sugars.get(0).isEmpty():
            try:
                _smiles = _smiles_generator.create(_aglycone_and_sugars.get(0))
                return str(_smiles)
            except Exception as e:
                raise Exception(f"{str(e)}")
        else:
            return ""
    else:
        return "No Linear sugar found"


def remove_circular_sugars(
    molecule: any,
    gly_bond: bool = False,
    only_terminal: bool = True,
    preservation_mode: preservation_modes_enum = preservation_modes_enum.HEAVY_ATOM_COUNT,
    preservation_threshold: int = 5,
    oxygen_atoms: bool = True,
    oxygen_atoms_threshold: float = 0.5,
    spiro_sugars: bool = False,
    keto_sugars: bool = False,
    mark_attach_points: bool = False,
) -> str:
    """
    Removes circular sugars from a given CDK IAtomContainer object using the Sugar Removal/Detection Utility.

    Args:
        molecule (IAtomContainer): CDK molecule object.
        gly_bond (bool): Whether to consider only circular sugars with glycosidic bonds in the removal process. Default is False.
        only_terminal (bool): Whether to remove only terminal circular sugars. Default is True.
        preservation_mode (preservation_modes_enum): Mode to determine which disconnected structures to preserve. Default is
                                                     preservation_modes_enum.HEAVY_ATOM_COUNT.
        preservation_threshold (int): Threshold value for the selected preservation mode. Default is 5 (heavy atoms).
        oxygen_atoms (bool): Whether to consider only circular sugars with a sufficient number of exocyclic oxygen atoms in the
                             removal process (see oxygen_atoms_threshold). Default is True.
        oxygen_atoms_threshold (float): A number giving the minimum attached exocyclic oxygen atoms to atom number in the ring
                                        ratio a circular sugar needs to have to be considered in the removal process. Default is
                                        0.5 (a 6-membered ring needs at least 3 attached exocyclic oxygen atoms). Must be positive!
        spiro_sugars (bool): Whether spiro rings (rings that share one atom with another cycle) should be included in the circular
                             sugar detection. Default is False.
        keto_sugars (bool): Whether circular sugars with keto groups should be detected. Default is False.
        mark_attach_points (bool): Whether to mark the attachment points of removed sugars with a dummy atom. Default is False.

    Returns:
        str: The SMILES string with circular sugars removed, or a message indicating that no circular sugars were
        found ("No Circular sugars found").

    Raises:
        ValueError: If one of the numeric parameters is not valid.
        Exception: If an error occurs during the output SMILES generation process.
    """
    _smi_flavor = cdk.JClass(SMI_FLAVOR_PATH_AND_CLASS_NAME)
    _smiles_generator = cdk.JClass(SMILES_GENERATOR_PATH_AND_CLASS_NAME)(
        _smi_flavor.Absolute,
    )
    _sugar_detection_utility = cdk.JClass(SDU_PATH_BASE + "." + SDU_CLASS_NAME)(
        SCOB_CLASS.getInstance(),
    )

    _sugar_detection_utility.setDetectCircularSugarsOnlyWithOGlycosidicBondSetting(
        gly_bond
    )
    _sugar_detection_utility.setRemoveOnlyTerminalSugarsSetting(only_terminal)
    _sru_preservation_modes_enum = cdk.JClass(
        SDU_PATH_BASE + "." + "SugarRemovalUtility$PreservationMode"
    )
    if preservation_mode == preservation_modes_enum.ALL:
        _sugar_detection_utility.setPreservationModeSetting(
            _sru_preservation_modes_enum.ALL
        )
    elif preservation_mode == preservation_modes_enum.HEAVY_ATOM_COUNT:
        _sugar_detection_utility.setPreservationModeSetting(
            _sru_preservation_modes_enum.HEAVY_ATOM_COUNT
        )
    elif preservation_mode == preservation_modes_enum.MOLECULAR_WEIGHT:
        _sugar_detection_utility.setPreservationModeSetting(
            _sru_preservation_modes_enum.MOLECULAR_WEIGHT
        )
    else:
        raise ValueError("Invalid preservation_mode specified.")
    if preservation_threshold < 0:
        raise ValueError("preservation_threshold must be a non-negative integer.")
    _sugar_detection_utility.setPreservationModeThresholdSetting(preservation_threshold)
    _sugar_detection_utility.setDetectCircularSugarsOnlyWithEnoughExocyclicOxygenAtomsSetting(
        oxygen_atoms
    )
    _is_valid = _sugar_detection_utility.setExocyclicOxygenAtomsToAtomsInRingRatioThresholdSetting(
        oxygen_atoms_threshold
    )
    if not _is_valid:
        raise ValueError("oxygenAtomsThreshold must be a positive number.")
    _sugar_detection_utility.setDetectSpiroRingsAsCircularSugarsSetting(spiro_sugars)
    _sugar_detection_utility.setDetectCircularSugarsWithKetoGroupsSetting(keto_sugars)

    _has_circular_sugars = _sugar_detection_utility.hasCircularSugars(molecule)

    if _has_circular_sugars:
        _aglycone_and_sugars = _sugar_detection_utility.copyAndExtractAglyconeAndSugars(
            molecule,
            True,
            False,
            mark_attach_points,
        )
        if not _aglycone_and_sugars.get(0).isEmpty():
            try:
                _smiles = _smiles_generator.create(_aglycone_and_sugars.get(0))
                return str(_smiles)
            except Exception as e:
                raise Exception(f"{str(e)}")
        else:
            return ""
    else:
        return "No Circular sugars found"


def remove_linear_and_circular_sugars(
    molecule: any,
    gly_bond: bool = False,
    only_terminal: bool = True,
    preservation_mode: preservation_modes_enum = preservation_modes_enum.HEAVY_ATOM_COUNT,
    preservation_threshold: int = 5,
    oxygen_atoms: bool = True,
    oxygen_atoms_threshold: float = 0.5,
    linear_sugars_in_rings: bool = False,
    linear_sugars_min_size: int = 4,
    linear_sugars_max_size: int = 7,
    linear_acidic_sugars: bool = False,
    spiro_sugars: bool = False,
    keto_sugars: bool = False,
    mark_attach_points: bool = False,
) -> str:
    """
    Removes circular and linear sugars from a given CDK IAtomContainer object using the Sugar Removal/Detection Utility.

    Args:
        molecule (IAtomContainer): CDK molecule object.
        gly_bond (bool): Whether to consider only circular sugars with glycosidic bonds in the removal process. Default is False.
        only_terminal (bool): Whether to remove only terminal sugars. Default is True.
        preservation_mode (preservation_modes_enum): Mode to determine which disconnected structures to preserve. Default is
                                                     preservation_modes_enum.HEAVY_ATOM_COUNT.
        preservation_threshold (int): Threshold value for the selected preservation mode. Default is 5 (heavy atoms).
        oxygen_atoms (bool): Whether to consider only circular sugars with a sufficient number of exocyclic oxygen atoms in the
                             removal process (see oxygen_atoms_threshold). Default is True.
        oxygen_atoms_threshold (float): A number giving the minimum attached exocyclic oxygen atoms to atom number in the ring
                                        ratio a circular sugar needs to have to be considered in the removal process. Default is
                                        0.5 (a 6-membered ring needs at least 3 attached exocyclic oxygen atoms). Must be positive!
        linear_sugars_in_rings (bool): Whether to consider linear sugars in rings. Default is False.
        linear_sugars_min_size (int): Minimum size of linear sugars to consider. Default is 4. Must be positive and higher than or
                                      equal to 1 and also smaller than the linear sugar candidate maximum size.
        linear_sugars_max_size (int): Maximum size of linear sugars to consider. Default is 7. Must be positive and higher than or
                                      equal to 1 and also higher than the linear sugar candidate minimum size.
        linear_acidic_sugars (bool): Whether to consider linear acidic sugars. Default is False.
        spiro_sugars (bool): Whether spiro rings (rings that share one atom with another cycle) should be included in the circular
                             sugar detection. Default is False.
        keto_sugars (bool): Whether circular sugars with keto groups should be detected. Default is False.
        mark_attach_points (bool): Whether to mark the attachment points of removed sugars with a dummy atom. Default is False.

    Returns:
        str: The SMILES string with circular and linear sugars removed, or a message indicating that no sugars were
        found ("No Linear or Circular sugars found").

    Raises:
        ValueError: If one of the numeric parameters is not valid.
        Exception: If an error occurs during the output SMILES generation process.
    """
    _smi_flavor = cdk.JClass(SMI_FLAVOR_PATH_AND_CLASS_NAME)
    _smiles_generator = cdk.JClass(SMILES_GENERATOR_PATH_AND_CLASS_NAME)(
        _smi_flavor.Absolute,
    )
    _sugar_detection_utility = cdk.JClass(SDU_PATH_BASE + "." + SDU_CLASS_NAME)(
        SCOB_CLASS.getInstance(),
    )

    _sugar_detection_utility.setDetectCircularSugarsOnlyWithOGlycosidicBondSetting(
        gly_bond
    )
    _sugar_detection_utility.setRemoveOnlyTerminalSugarsSetting(only_terminal)
    _sru_preservation_modes_enum = cdk.JClass(
        SDU_PATH_BASE + "." + "SugarRemovalUtility$PreservationMode"
    )
    if preservation_mode == preservation_modes_enum.ALL:
        _sugar_detection_utility.setPreservationModeSetting(
            _sru_preservation_modes_enum.ALL
        )
    elif preservation_mode == preservation_modes_enum.HEAVY_ATOM_COUNT:
        _sugar_detection_utility.setPreservationModeSetting(
            _sru_preservation_modes_enum.HEAVY_ATOM_COUNT
        )
    elif preservation_mode == preservation_modes_enum.MOLECULAR_WEIGHT:
        _sugar_detection_utility.setPreservationModeSetting(
            _sru_preservation_modes_enum.MOLECULAR_WEIGHT
        )
    else:
        raise ValueError("Invalid preservation_mode specified.")
    if preservation_threshold < 0:
        raise ValueError("preservation_threshold must be a non-negative integer.")
    _sugar_detection_utility.setPreservationModeThresholdSetting(preservation_threshold)
    _sugar_detection_utility.setDetectCircularSugarsOnlyWithEnoughExocyclicOxygenAtomsSetting(
        oxygen_atoms
    )
    _is_valid = _sugar_detection_utility.setExocyclicOxygenAtomsToAtomsInRingRatioThresholdSetting(
        oxygen_atoms_threshold
    )
    if not _is_valid:
        raise ValueError("oxygenAtomsThreshold must be a positive number.")
    _sugar_detection_utility.setDetectLinearSugarsInRingsSetting(linear_sugars_in_rings)
    if (
        linear_sugars_min_size < 0
        or linear_sugars_max_size < 1
        or linear_sugars_min_size >= linear_sugars_max_size
    ):
        raise ValueError(
            "linearSugarsMinSize and linearSugarsMaxSize must be positive integers, with "
            "linearSugarsMinSize being smaller than linearSugarsMaxSize."
        )
    _sugar_detection_utility.setLinearSugarCandidateMinSizeSetting(
        linear_sugars_min_size
    )
    _sugar_detection_utility.setLinearSugarCandidateMaxSizeSetting(
        linear_sugars_max_size
    )
    _sugar_detection_utility.setDetectLinearAcidicSugarsSetting(linear_acidic_sugars)
    _sugar_detection_utility.setDetectSpiroRingsAsCircularSugarsSetting(spiro_sugars)
    _sugar_detection_utility.setDetectCircularSugarsWithKetoGroupsSetting(keto_sugars)

    _has_circular_or_linear_sugars = _sugar_detection_utility.hasCircularOrLinearSugars(
        molecule,
    )

    if _has_circular_or_linear_sugars:
        _aglycone_and_sugars = _sugar_detection_utility.copyAndExtractAglyconeAndSugars(
            molecule,
            True,
            True,
            mark_attach_points,
        )
        if not _aglycone_and_sugars.get(0).isEmpty():
            try:
                _smiles = _smiles_generator.create(_aglycone_and_sugars.get(0))
                return str(_smiles)
            except Exception as e:
                raise Exception(f"{str(e)}")
        else:
            return ""
    else:
        return "No Linear or Circular sugars found"


def extract_aglycone_and_sugars(
    molecule: any,
    extract_circular_sugars: bool = True,
    extract_linear_sugars: bool = False,
    gly_bond: bool = False,
    only_terminal: bool = True,
    preservation_mode: preservation_modes_enum = preservation_modes_enum.HEAVY_ATOM_COUNT,
    preservation_threshold: int = 5,
    oxygen_atoms: bool = True,
    oxygen_atoms_threshold: float = 0.5,
    linear_sugars_in_rings: bool = False,
    linear_sugars_min_size: int = 4,
    linear_sugars_max_size: int = 7,
    linear_acidic_sugars: bool = False,
    spiro_sugars: bool = False,
    keto_sugars: bool = False,
    mark_attach_points: bool = False,
    post_process_sugars: bool = False,
    limit_post_process_by_size: bool = False,
) -> tuple:
    """
    Extracts aglycone and sugar components from a given CDK IAtomContainer object using the Sugar Removal/Detection Utility.
    Depending on the parameters, it can extract circular sugars, linear sugars, or both.

    Args:
        molecule (IAtomContainer): CDK molecule object.
        extract_circular_sugars (bool): Whether to extract circular sugars. Default is True.
        extract_linear_sugars (bool): Whether to extract linear sugars. Default is False.
        gly_bond (bool): Whether to consider only circular sugars with glycosidic bonds in the removal process. Default is False.
        only_terminal (bool): Whether to remove only terminal sugars. Default is True.
        preservation_mode (preservation_modes_enum): Mode to determine which disconnected structures to preserve. Default is
                                                     preservation_modes_enum.HEAVY_ATOM_COUNT.
        preservation_threshold (int): Threshold value for the selected preservation mode. Default is 5 (heavy atoms).
        oxygen_atoms (bool): Whether to consider only circular sugars with a sufficient number of exocyclic oxygen atoms in the
                             removal process (see oxygen_atoms_threshold). Default is True.
        oxygen_atoms_threshold (float): A number giving the minimum attached exocyclic oxygen atoms to atom number in the ring
                                        ratio a circular sugar needs to have to be considered in the removal process. Default is
                                        0.5 (a 6-membered ring needs at least 3 attached exocyclic oxygen atoms). Must be positive!
        linear_sugars_in_rings (bool): Whether to consider linear sugars in rings. Default is False.
        linear_sugars_min_size (int): Minimum size of linear sugars to consider. Default is 4. Must be positive and higher than or
                                      equal to 1 and also smaller than the linear sugar candidate maximum size.
        linear_sugars_max_size (int): Maximum size of linear sugars to consider. Default is 7. Must be positive and higher than or
                                      equal to 1 and also higher than the linear sugar candidate minimum size.
        linear_acidic_sugars (bool): Whether to consider linear acidic sugars. Default is False.
        spiro_sugars (bool): Whether spiro rings (rings that share one atom with another cycle) should be included in the circular
                             sugar detection. Default is False.
        keto_sugars (bool): Whether circular sugars with keto groups should be detected. Default is False.
        mark_attach_points (bool): Whether to mark the attachment points of removed sugars with a dummy atom. Default is False.
        post_process_sugars (bool): Whether the extracted sugar moieties should be post-processed, i.e. bond splitting (O-glycosidic,
                                    ether, ester, peroxide) to separate the individual sugars, before being output. Default is False.
        limit_post_process_by_size (bool): Whether the post-processing of extracted sugar moieties should be limited to structures
                                           bigger than a defined size (see preservation mode (threshold)) to preserve smaller
                                           modifications. Default is False.

    Returns:
        tuple: The SMILES strings of the generated aglycone (position 0) and sugars (position(s) 1..n).
               If no sugars were found, the tuple only contains the aglycone SMILES string.

    Raises:
        ValueError: If one of the numeric parameters is not valid.
        Exception: If an error occurs during the output SMILES generation process.
    """
    _smi_flavor = cdk.JClass(SMI_FLAVOR_PATH_AND_CLASS_NAME)
    _smiles_generator = cdk.JClass(SMILES_GENERATOR_PATH_AND_CLASS_NAME)(
        _smi_flavor.Absolute,
    )
    _sugar_detection_utility = cdk.JClass(SDU_PATH_BASE + "." + SDU_CLASS_NAME)(
        SCOB_CLASS.getInstance(),
    )

    _sugar_detection_utility.setDetectCircularSugarsOnlyWithOGlycosidicBondSetting(
        gly_bond
    )
    _sugar_detection_utility.setRemoveOnlyTerminalSugarsSetting(only_terminal)
    _sru_preservation_modes_enum = cdk.JClass(
        SDU_PATH_BASE + "." + "SugarRemovalUtility$PreservationMode"
    )
    if preservation_mode == preservation_modes_enum.ALL:
        _sugar_detection_utility.setPreservationModeSetting(
            _sru_preservation_modes_enum.ALL
        )
    elif preservation_mode == preservation_modes_enum.HEAVY_ATOM_COUNT:
        _sugar_detection_utility.setPreservationModeSetting(
            _sru_preservation_modes_enum.HEAVY_ATOM_COUNT
        )
    elif preservation_mode == preservation_modes_enum.MOLECULAR_WEIGHT:
        _sugar_detection_utility.setPreservationModeSetting(
            _sru_preservation_modes_enum.MOLECULAR_WEIGHT
        )
    else:
        raise ValueError("Invalid preservation_mode specified.")
    if preservation_threshold < 0:
        raise ValueError("preservation_threshold must be a non-negative integer.")
    _sugar_detection_utility.setPreservationModeThresholdSetting(preservation_threshold)
    _sugar_detection_utility.setDetectCircularSugarsOnlyWithEnoughExocyclicOxygenAtomsSetting(
        oxygen_atoms
    )
    _is_valid = _sugar_detection_utility.setExocyclicOxygenAtomsToAtomsInRingRatioThresholdSetting(
        oxygen_atoms_threshold
    )
    if not _is_valid:
        raise ValueError("oxygenAtomsThreshold must be a positive number.")
    _sugar_detection_utility.setDetectLinearSugarsInRingsSetting(linear_sugars_in_rings)
    if (
        linear_sugars_min_size < 0
        or linear_sugars_max_size < 1
        or linear_sugars_min_size >= linear_sugars_max_size
    ):
        raise ValueError(
            "linearSugarsMinSize and linearSugarsMaxSize must be positive integers, with linearSugarsMinSize being smaller than linearSugarsMaxSize."
        )
    _sugar_detection_utility.setLinearSugarCandidateMinSizeSetting(
        linear_sugars_min_size
    )
    _sugar_detection_utility.setLinearSugarCandidateMaxSizeSetting(
        linear_sugars_max_size
    )
    _sugar_detection_utility.setDetectLinearAcidicSugarsSetting(linear_acidic_sugars)
    _sugar_detection_utility.setDetectSpiroRingsAsCircularSugarsSetting(spiro_sugars)
    _sugar_detection_utility.setDetectCircularSugarsWithKetoGroupsSetting(keto_sugars)

    _aglycone_and_sugars = _sugar_detection_utility.copyAndExtractAglyconeAndSugars(
        molecule,
        extract_circular_sugars,
        extract_linear_sugars,
        mark_attach_points,
        post_process_sugars,
        limit_post_process_by_size,
    )
    _result = []
    for _structure in _aglycone_and_sugars:
        if _structure.isEmpty():
            _result.append("")
            continue
        try:
            _smiles = _smiles_generator.create(_structure)
            _result.append(str(_smiles))
        except Exception as e:
            raise Exception(f"{str(e)}")
    return tuple(_result)


def get_aglycone_and_sugar_indices(
    molecule: any,
    extract_circular_sugars: bool = True,
    extract_linear_sugars: bool = False,
    gly_bond: bool = False,
    only_terminal: bool = True,
    preservation_mode: preservation_modes_enum = preservation_modes_enum.HEAVY_ATOM_COUNT,
    preservation_threshold: int = 5,
    oxygen_atoms: bool = True,
    oxygen_atoms_threshold: float = 0.5,
    linear_sugars_in_rings: bool = False,
    linear_sugars_min_size: int = 4,
    linear_sugars_max_size: int = 7,
    linear_acidic_sugars: bool = False,
    spiro_sugars: bool = False,
    keto_sugars: bool = False,
    mark_attach_points: bool = False,
    post_process_sugars: bool = False,
    limit_post_process_by_size: bool = False,
) -> list[list[int]]:
    """
    Extracts aglycone and sugar components from a given CDK IAtomContainer object using the Sugar Removal/Detection Utility
    and returns their atom indices.
    Depending on the parameters, it can extract circular sugars, linear sugars, or both.

    Args:
        molecule (IAtomContainer): CDK molecule object.
        extract_circular_sugars (bool): Whether to extract circular sugars. Default is True.
        extract_linear_sugars (bool): Whether to extract linear sugars. Default is False.
        gly_bond (bool): Whether to consider only circular sugars with glycosidic bonds in the removal process. Default is False.
        only_terminal (bool): Whether to remove only terminal sugars. Default is True.
        preservation_mode (preservation_modes_enum): Mode to determine which disconnected structures to preserve. Default is
                                                     preservation_modes_enum.HEAVY_ATOM_COUNT.
        preservation_threshold (int): Threshold value for the selected preservation mode. Default is 5 (heavy atoms).
        oxygen_atoms (bool): Whether to consider only circular sugars with a sufficient number of exocyclic oxygen atoms in the
                             removal process (see oxygen_atoms_threshold). Default is True.
        oxygen_atoms_threshold (float): A number giving the minimum attached exocyclic oxygen atoms to atom number in the ring
                                        ratio a circular sugar needs to have to be considered in the removal process. Default is
                                        0.5 (a 6-membered ring needs at least 3 attached exocyclic oxygen atoms). Must be positive!
        linear_sugars_in_rings (bool): Whether to consider linear sugars in rings. Default is False.
        linear_sugars_min_size (int): Minimum size of linear sugars to consider. Default is 4. Must be positive and higher than or
                                      equal to 1 and also smaller than the linear sugar candidate maximum size.
        linear_sugars_max_size (int): Maximum size of linear sugars to consider. Default is 7. Must be positive and higher than or
                                      equal to 1 and also higher than the linear sugar candidate minimum size.
        linear_acidic_sugars (bool): Whether to consider linear acidic sugars. Default is False.
        spiro_sugars (bool): Whether spiro rings (rings that share one atom with another cycle) should be included in the circular
                             sugar detection. Default is False.
        keto_sugars (bool): Whether circular sugars with keto groups should be detected. Default is False.
        mark_attach_points (bool): Whether to mark the attachment points of removed sugars with a dummy atom. Default is False.
        post_process_sugars (bool): Whether the extracted sugar moieties should be post-processed, i.e. bond splitting (O-glycosidic,
                                    ether, ester, peroxide) to separate the individual sugars, before being output. Default is False.
        limit_post_process_by_size (bool): Whether the post-processing of extracted sugar moieties should be limited to structures
                                           bigger than a defined size (see preservation mode (threshold)) to preserve smaller
                                           modifications. Default is False.

    Returns:
        list: The atom indices of the generated aglycone (position 0) and sugars (position(s) 1..n).
               If no sugars were found, the list only contains the aglycone atom indices.

    Raises:
        ValueError: If one of the numeric parameters is not valid.
        Exception: If an error occurs during the output generation process.
    """
    _sugar_detection_utility = cdk.JClass(SDU_PATH_BASE + "." + SDU_CLASS_NAME)(
        SCOB_CLASS.getInstance(),
    )
    _input_atom_to_aglycone_atom_map = cdk.JClass("java.util.HashMap")(
        int((molecule.getAtomCount() / 0.75) + 2), 0.75
    )
    _input_atom_to_sugar_atom_map = cdk.JClass("java.util.HashMap")(
        int((molecule.getAtomCount() / 0.75) + 2), 0.75
    )

    _sugar_detection_utility.setDetectCircularSugarsOnlyWithOGlycosidicBondSetting(
        gly_bond
    )
    _sugar_detection_utility.setRemoveOnlyTerminalSugarsSetting(only_terminal)
    _sru_preservation_modes_enum = cdk.JClass(
        SDU_PATH_BASE + "." + "SugarRemovalUtility$PreservationMode"
    )
    if preservation_mode == preservation_modes_enum.ALL:
        _sugar_detection_utility.setPreservationModeSetting(
            _sru_preservation_modes_enum.ALL
        )
    elif preservation_mode == preservation_modes_enum.HEAVY_ATOM_COUNT:
        _sugar_detection_utility.setPreservationModeSetting(
            _sru_preservation_modes_enum.HEAVY_ATOM_COUNT
        )
    elif preservation_mode == preservation_modes_enum.MOLECULAR_WEIGHT:
        _sugar_detection_utility.setPreservationModeSetting(
            _sru_preservation_modes_enum.MOLECULAR_WEIGHT
        )
    else:
        raise ValueError("Invalid preservation_mode specified.")
    if preservation_threshold < 0:
        raise ValueError("preservation_threshold must be a non-negative integer.")
    _sugar_detection_utility.setPreservationModeThresholdSetting(preservation_threshold)
    _sugar_detection_utility.setDetectCircularSugarsOnlyWithEnoughExocyclicOxygenAtomsSetting(
        oxygen_atoms
    )
    _is_valid = _sugar_detection_utility.setExocyclicOxygenAtomsToAtomsInRingRatioThresholdSetting(
        oxygen_atoms_threshold
    )
    if not _is_valid:
        raise ValueError("oxygenAtomsThreshold must be a positive number.")
    _sugar_detection_utility.setDetectLinearSugarsInRingsSetting(linear_sugars_in_rings)
    if (
        linear_sugars_min_size < 0
        or linear_sugars_max_size < 1
        or linear_sugars_min_size >= linear_sugars_max_size
    ):
        raise ValueError(
            "linearSugarsMinSize and linearSugarsMaxSize must be positive integers, with linearSugarsMinSize being smaller than linearSugarsMaxSize."
        )
    _sugar_detection_utility.setLinearSugarCandidateMinSizeSetting(
        linear_sugars_min_size
    )
    _sugar_detection_utility.setLinearSugarCandidateMaxSizeSetting(
        linear_sugars_max_size
    )
    _sugar_detection_utility.setDetectLinearAcidicSugarsSetting(linear_acidic_sugars)
    _sugar_detection_utility.setDetectSpiroRingsAsCircularSugarsSetting(spiro_sugars)
    _sugar_detection_utility.setDetectCircularSugarsWithKetoGroupsSetting(keto_sugars)

    _aglycone_and_sugars = _sugar_detection_utility.copyAndExtractAglyconeAndSugars(
        molecule,
        extract_circular_sugars,
        extract_linear_sugars,
        mark_attach_points,
        post_process_sugars,
        limit_post_process_by_size,
        _input_atom_to_aglycone_atom_map,
        cdk.JClass("java.util.HashMap")(int((molecule.getBondCount() / 0.75) + 2), 0.75),
        _input_atom_to_sugar_atom_map,
        cdk.JClass("java.util.HashMap")(int((molecule.getBondCount() / 0.75) + 2), 0.75),
    )
    _result = []
    _aglycone_atom_indices = _sugar_detection_utility.getAtomIndicesOfGroup(molecule, _aglycone_and_sugars.get(0), _input_atom_to_aglycone_atom_map)
    _result.append([int(x) for x in _aglycone_atom_indices])
    for i in range(1, _aglycone_and_sugars.size()):
        try:
            _sugar_atom_indices = _sugar_detection_utility.getAtomIndicesOfGroup(molecule, _aglycone_and_sugars.get(i), _input_atom_to_sugar_atom_map)
            _result.append([int(x) for x in _sugar_atom_indices])
        except:
            _result.append([])
    return _result
