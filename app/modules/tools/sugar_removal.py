from __future__ import annotations

import app.modules.toolkits.cdk_wrapper as cdk

SCOB_CLASS = cdk.JClass(cdk.cdk_base + ".silent.SilentChemObjectBuilder")

SDU_PATH_BASE = cdk.cdk_base + ".tools"

SDU_CLASS_NAME = "SugarDetectionUtility"

SMI_FLAVOR_PATH_AND_CLASS_NAME = cdk.cdk_base + ".smiles.SmiFlavor"

SMILES_GENERATOR_PATH_AND_CLASS_NAME = cdk.cdk_base + ".smiles.SmilesGenerator"

def get_sugar_info(molecule: any) -> tuple:
    """Analyzes a molecule represented by a CDK IAtomContainer object to determine whether it contains sugars.

    This function utilizes the Sugar Removal/Detection Utility to check for the presence of circular or linear sugars.

    Args:
        molecule (IAtomContainer): CDK molecule object.

    Returns:
        tuple: A tuple containing two boolean values indicating whether the molecule has linear sugars
               and whether the molecule has circular sugars. If no sugars are found, both values will be False.
    """
    _sugar_detection_utility = cdk.JClass(SDU_PATH_BASE + "." + SDU_CLASS_NAME)(
        SCOB_CLASS.getInstance(),
    )
    
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
