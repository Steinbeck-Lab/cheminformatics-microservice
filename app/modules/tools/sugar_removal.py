from __future__ import annotations

import app.modules.toolkits.cdk_wrapper as cdk

SCOB_CLASS = cdk.JClass(cdk.cdk_base + ".silent.SilentChemObjectBuilder")

SDU_PATH_BASE = cdk.cdk_base + ".tools"

SDU_CLASS_NAME = "SugarDetectionUtility"

SMI_FLAVOR_PATH_AND_CLASS_NAME = cdk.cdk_base + ".smiles.SmiFlavor"

SMILES_GENERATOR_PATH_AND_CLASS_NAME = cdk.cdk_base + ".smiles.SmilesGenerator"

def get_sugar_info(molecule: any, 
                   gly_bond: bool = False, 
                   oxygen_atoms : bool = True,
                   oxygen_atoms_threshold : float = 0.5,
                   linear_sugars_in_rings : bool = False,
                   linear_sugars_min_size : int = 4,
                   linear_sugars_max_size : int = 7,
                   linear_acidic_sugars : bool = False,
                   spiro_sugars : bool = False,
                   keto_sugars : bool = False
                   ) -> tuple:
    """Analyzes a molecule represented by a CDK IAtomContainer object to determine whether it contains sugars.

    This function utilizes the Sugar Removal/Detection Utility to check for the presence of circular or linear sugars.

    Args:
        molecule (IAtomContainer): CDK molecule object.
        glyBond (bool): Whether to consider only circular sugars with glycosidic bonds in the analysis. Default is False.
        oxygenAtoms (bool): Whether to consider only circular sugars with a sufficient number of exocyclic oxygen atoms in the analysis (see oxygenAtomsThreshold). Default is True.
        oxygenAtomsThreshold (float): A number giving the minimum attached exocyclic oxygen atoms to atom number in the ring ratio a circular sugar needs to have to be considered in the analysis. Default is 0.5 (a 6-membered ring needs at least 3 attached exocyclic oxygen atoms). Must be positive!
        linearSugarsInRings (bool): Whether to consider linear sugars in rings. Default is False.
        linearSugarsMinSize (int): Minimum size of linear sugars to consider. Default is 4. Must be positive and higher than or equal to 1 and also smaller than the linear sugar candidate maximum size.
        linearSugarsMaxSize (int): Maximum size of linear sugars to consider. Default is 7. Must be positive and higher than or equal to 1 and also higher than the linear sugar candidate minimum size.
        linearAcidicSugars (bool): Whether to consider linear acidic sugars. Default is False.
        spiroSugars (bool): Whether spiro rings (rings that share one atom with another cycle) should be included in the circular sugar detection. Default is False.
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

    _sugar_detection_utility.setDetectCircularSugarsOnlyWithOGlycosidicBondSetting(gly_bond)
    _sugar_detection_utility.setDetectCircularSugarsOnlyWithEnoughExocyclicOxygenAtomsSetting(oxygen_atoms)
    _is_valid = _sugar_detection_utility.setExocyclicOxygenAtomsToAtomsInRingRatioThresholdSetting(oxygen_atoms_threshold)
    if not _is_valid:
        raise ValueError("oxygenAtomsThreshold must be a positive number.")
    _sugar_detection_utility.setDetectLinearSugarsInRingsSetting(linear_sugars_in_rings)
    if linear_sugars_min_size < 0 or linear_sugars_max_size < 1 or linear_sugars_min_size >= linear_sugars_max_size:
        raise ValueError("linearSugarsMinSize and linearSugarsMaxSize must be positive integers, with linearSugarsMinSize being smaller than linearSugarsMaxSize.")
    _sugar_detection_utility.setLinearSugarCandidateMinSizeSetting(linear_sugars_min_size)
    _sugar_detection_utility.setLinearSugarCandidateMaxSizeSetting(linear_sugars_max_size)
    _sugar_detection_utility.setDetectLinearAcidicSugarsSetting(linear_acidic_sugars)
    _sugar_detection_utility.setDetectSpiroRingsAsCircularSugarsSetting(spiro_sugars)
    _sugar_detection_utility.setDetectCircularSugarsWithKetoGroupsSetting(keto_sugars)

    _has_linear_sugars = _sugar_detection_utility.hasLinearSugars(molecule)
    _has_circular_sugars = _sugar_detection_utility.hasCircularSugars(molecule)
    return _has_linear_sugars, _has_circular_sugars


def remove_linear_sugar(molecule: any) -> str:
    """Removes linear sugars from a given CDK IAtomContainer object using the Sugar Removal/Detection Utility.

    Args:
        molecule (IAtomContainer): CDK molecule object.

    Returns:
        str: The SMILES string with linear sugars removed, or a message indicating that no linear sugars were 
        found ("No Linear sugar found").
    Raises:
        Exception: If an error occurs during the output SMILES generation process.
    """
    _smi_flavor = cdk.JClass(SMI_FLAVOR_PATH_AND_CLASS_NAME)
    _smiles_generator = cdk.JClass(SMILES_GENERATOR_PATH_AND_CLASS_NAME)(
        _smi_flavor.Absolute,
    )
    _sugar_detection_utility = cdk.JClass(SDU_PATH_BASE + "." + SDU_CLASS_NAME)(
        SCOB_CLASS.getInstance(),
    )

    _has_linear_sugars = _sugar_detection_utility.hasLinearSugars(molecule)

    if _has_linear_sugars:
        _aglycone_and_sugars = _sugar_detection_utility.copyAndExtractAglyconeAndSugars(
            molecule,
            False,
            True,
        )
        if not _aglycone_and_sugars.isEmpty():
            try:
                _smiles = _smiles_generator.create(_aglycone_and_sugars.get(0))
                return str(_smiles)
            except Exception as e:
                raise Exception(f"{str(e)}")
        else:
            return ""
    else:
        return "No Linear sugar found"


def remove_circular_sugar(molecule: any) -> str:
    """Removes circular sugars from a given CDK IAtomContainer object using the Sugar Removal/Detection Utility.

    Args:
        molecule (IAtomContainer): CDK molecule object.

    Returns:
        str: The SMILES string with circular sugars removed, or a message indicating that no circular sugars were 
        found ("No Circular sugars found").
    Raises:
        Exception: If an error occurs during the output SMILES generation process.
    """
    _smi_flavor = cdk.JClass(SMI_FLAVOR_PATH_AND_CLASS_NAME)
    _smiles_generator = cdk.JClass(SMILES_GENERATOR_PATH_AND_CLASS_NAME)(
        _smi_flavor.Absolute,
    )
    _sugar_detection_utility = cdk.JClass(SDU_PATH_BASE + "." + SDU_CLASS_NAME)(
        SCOB_CLASS.getInstance(),
    )

    _has_circular_sugars = _sugar_detection_utility.hasCircularSugars(molecule)

    if _has_circular_sugars:
        _aglycone_and_sugars = _sugar_detection_utility.copyAndExtractAglyconeAndSugars(
            molecule,
            True,
            False,
        )
        if not _aglycone_and_sugars.isEmpty():
            try:
                _smiles = _smiles_generator.create(_aglycone_and_sugars.get(0))
                return str(_smiles)
            except Exception as e:
                raise Exception(f"{str(e)}")
        else:
            return ""
    else:
        return "No Circular sugars found"


def remove_linear_and_circular_sugar(molecule: any):
    """Removes circular and linearsugars from a given CDK IAtomContainer object using the Sugar Removal/Detection Utility.

    Args:
        molecule (IAtomContainer): CDK molecule object.

    Returns:
        str: The SMILES string with circular and linear sugars removed, or a message indicating that no sugars were 
        found ("No Linear or Circular sugars found").
    Raises:
        Exception: If an error occurs during the output SMILES generation process.
    """
    _smi_flavor = cdk.JClass(SMI_FLAVOR_PATH_AND_CLASS_NAME)
    _smiles_generator = cdk.JClass(SMILES_GENERATOR_PATH_AND_CLASS_NAME)(
        _smi_flavor.Absolute,
    )
    _sugar_detection_utility = cdk.JClass(SDU_PATH_BASE + "." + SDU_CLASS_NAME)(
        SCOB_CLASS.getInstance(),
    )

    _has_circular_or_linear_sugars = _sugar_detection_utility.hasCircularOrLinearSugars(
        molecule,
    )

    if _has_circular_or_linear_sugars:
        _aglycone_and_sugars = _sugar_detection_utility.copyAndExtractAglyconeAndSugars(
            molecule,
            True,
            True,
        )
        if not _aglycone_and_sugars.isEmpty():
            try:
                _smiles = _smiles_generator.create(_aglycone_and_sugars.get(0))
                return str(_smiles)
            except Exception as e:
                raise Exception(f"{str(e)}")
        else:
            return ""
    else:
        return "No Linear or Circular sugars found"
