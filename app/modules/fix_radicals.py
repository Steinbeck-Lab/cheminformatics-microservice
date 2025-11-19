from __future__ import annotations

from jpype import JClass
from app.modules.cdk_depict.radical_perception import perceive_radicals
from app.modules.toolkits.cdk_wrapper import get_canonical_SMILES
from app.modules.toolkits.cdk_wrapper import cdk_base


def fixradicals(molecule: any) -> dict:
    """Fix radicals in the given molecule using CDK.

    IMPORTANT: This function expects ISingleElectron objects to already exist
    on the molecule. You must call perceive_radicals() BEFORE calling this function.

    This function handles radicals (single electrons) in molecules by:
    1. Identifying atoms with single electrons (C, N, O only)
    2. Resetting their properties (hybridization, valency, formal neighbour count)
    3. Incrementing implicit hydrogen count
    4. Removing the single electron
    5. Re-perceiving atom types and adding implicit hydrogens

    Args:
        molecule (IAtomContainer): molecule given by the user.
    Returns:
        dict: A dictionary containing:
            - fixed_smiles (str): The canonical SMILES with fixed radicals
            - radicals_detected (int): Number of radicals detected
            - radicals_fixed (int): Number of radicals successfully fixed

    Raises:
        ValueError: If molecule is None or empty.
    """

    # Null check
    if molecule is None:
        raise ValueError("Given molecule is None.")

    # Get CDK classes
    AtomContainerManipulator = JClass(
        cdk_base + ".tools.manipulator.AtomContainerManipulator",
    )
    CDKHydrogenAdder = JClass(cdk_base + ".tools.CDKHydrogenAdder")
    IElement = JClass(cdk_base + ".interfaces.IElement")
    Integer = JClass("java.lang.Integer")
    perceive_radicals(molecule)

    # Track statistics
    radicals_detected = molecule.getSingleElectronCount()
    radicals_fixed = 0

    # Check if molecule has single electrons
    if radicals_detected > 0:

        # Fix properties of atoms that are radicals
        encountered_supported_radical = False
        i = 0

        # Use while loop because removeSingleElectron shifts indices
        while i < molecule.getSingleElectronCount():
            single_electron = molecule.getSingleElectron(i)
            atom = single_electron.getAtom()
            atomic_number = atom.getAtomicNumber()

            # Only handle C, N, O for now
            if atomic_number not in [IElement.C, IElement.N, IElement.O]:
                i += 1
                continue

            encountered_supported_radical = True
            radicals_fixed += 1

            # Setting to null now, will be re-detected correctly by
            # AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms below
            atom.setHybridization(None)
            atom.setValency(None)
            atom.setFormalNeighbourCount(None)

            # Increment implicit hydrogen count
            h_count = atom.getImplicitHydrogenCount()
            if h_count is None:
                h_count = 0
            atom.setImplicitHydrogenCount(Integer(h_count + 1))

            # Remove the single electron
            # Important: Don't increment i since removal shifts indices (mimics i-- in Java)
            molecule.removeSingleElectron(i)

        if encountered_supported_radical:
            # Needs to be redone now to set the correct atom types
            AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(molecule)
            # The newly added line to complete valences if fixing the radicals did not do the trick
            hydrogen_adder = CDKHydrogenAdder.getInstance(molecule.getBuilder())
            hydrogen_adder.addImplicitHydrogens(molecule)

    fixed_SMILES = get_canonical_SMILES(molecule)

    return {
        "fixed_smiles": fixed_SMILES,
        "radicals_detected": radicals_detected,
        "radicals_fixed": radicals_fixed,
    }
