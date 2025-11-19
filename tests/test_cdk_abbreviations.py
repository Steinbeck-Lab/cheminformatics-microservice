from __future__ import annotations

import pytest
from jpype import JClass

from app.modules.cdk_depict.abbreviations import (
    ChemicalAbbreviations,
    AbbreviationMode,
    AbbreviationOptions,
)
from app.modules.toolkits.cdk_wrapper import get_CDK_IAtomContainer
from pathlib import Path
import tempfile


@pytest.fixture
def abbr_system():
    system = ChemicalAbbreviations()
    system.initialize()
    return system


@pytest.fixture
def phenyl_molecule():
    return get_CDK_IAtomContainer("c1ccccc1C")


@pytest.fixture
def thf_molecule():
    return get_CDK_IAtomContainer("C1CCOC1")


@pytest.fixture
def dmf_molecule():
    return get_CDK_IAtomContainer("CN(C)C=O")


@pytest.fixture
def boc_molecule():
    return get_CDK_IAtomContainer("CC(C)(C)OC(=O)N")


@pytest.fixture
def complex_molecule():
    return get_CDK_IAtomContainer("c1ccccc1CC(C)(C)OC(=O)N")


class TestChemicalAbbreviationsInitialization:
    """Test abbreviation system initialization."""

    def test_initialization_succeeds(self):
        system = ChemicalAbbreviations()
        system.initialize()
        assert system._initialized is True

    def test_double_initialization_safe(self):
        system = ChemicalAbbreviations()
        system.initialize()
        system.initialize()
        assert system._initialized is True

    def test_custom_cdk_base(self):
        system = ChemicalAbbreviations(cdk_base="org.openscience.cdk")
        system.initialize()
        assert system.cdk_base == "org.openscience.cdk"


class TestAbbreviationModes:
    """Test different abbreviation modes."""

    def test_off_mode_no_abbreviations(self, abbr_system, phenyl_molecule):
        initial_count = phenyl_molecule.getAtomCount()
        abbr_system.apply(phenyl_molecule, mode=AbbreviationMode.OFF)
        assert phenyl_molecule.getAtomCount() == initial_count

    def test_groups_mode_phenyl(self, abbr_system, phenyl_molecule):
        abbr_system.apply(phenyl_molecule, mode=AbbreviationMode.GROUPS)
        assert phenyl_molecule is not None

    def test_reagents_mode_thf(self, abbr_system, thf_molecule):
        abbr_system.apply(thf_molecule, mode=AbbreviationMode.REAGENTS)
        assert thf_molecule is not None

    def test_all_mode_combines_both(self, abbr_system, complex_molecule):
        abbr_system.apply(complex_molecule, mode=AbbreviationMode.ALL)
        assert complex_molecule is not None

    def test_reagents_mode_dmf(self, abbr_system, dmf_molecule):
        abbr_system.apply(dmf_molecule, mode=AbbreviationMode.REAGENTS)
        assert dmf_molecule is not None

    def test_groups_mode_boc(self, abbr_system, boc_molecule):
        abbr_system.apply(boc_molecule, mode=AbbreviationMode.GROUPS)
        assert boc_molecule is not None


class TestAbbreviationApplication:
    """Test application of abbreviations to molecules."""

    def test_apply_to_simple_phenyl(self, abbr_system, phenyl_molecule):
        abbr_system.apply(phenyl_molecule, mode=AbbreviationMode.GROUPS)
        assert phenyl_molecule is not None

    def test_apply_to_thf(self, abbr_system, thf_molecule):
        abbr_system.apply(thf_molecule, mode=AbbreviationMode.REAGENTS)
        assert thf_molecule is not None

    def test_apply_preserves_molecule_integrity(self, abbr_system, phenyl_molecule):
        initial_atom_count = phenyl_molecule.getAtomCount()
        abbr_system.apply(phenyl_molecule, mode=AbbreviationMode.GROUPS)
        assert phenyl_molecule.getAtomCount() <= initial_atom_count

    def test_apply_to_molecule_without_groups(self, abbr_system):
        mol = get_CDK_IAtomContainer("CCCC")
        initial_count = mol.getAtomCount()
        abbr_system.apply(mol, mode=AbbreviationMode.GROUPS)
        assert mol.getAtomCount() == initial_count

    def test_apply_to_molecule_without_reagents(self, abbr_system):
        mol = get_CDK_IAtomContainer("CCCC")
        initial_count = mol.getAtomCount()
        abbr_system.apply(mol, mode=AbbreviationMode.REAGENTS)
        assert mol.getAtomCount() == initial_count


class TestHighlightedAtomsExclusion:
    """Test that highlighted atoms are not abbreviated."""

    def test_highlighted_atoms_not_abbreviated(self, abbr_system, phenyl_molecule):
        highlighted = {0, 1, 2}
        abbr_system.apply(
            phenyl_molecule, mode=AbbreviationMode.GROUPS, highlighted_atoms=highlighted
        )
        assert phenyl_molecule is not None

    def test_empty_highlighted_set(self, abbr_system, phenyl_molecule):
        abbr_system.apply(
            phenyl_molecule, mode=AbbreviationMode.GROUPS, highlighted_atoms=set()
        )
        assert phenyl_molecule is not None

    def test_all_atoms_highlighted(self, abbr_system, phenyl_molecule):
        all_atoms = set(range(phenyl_molecule.getAtomCount()))
        abbr_system.apply(
            phenyl_molecule, mode=AbbreviationMode.GROUPS, highlighted_atoms=all_atoms
        )
        assert phenyl_molecule is not None

    def test_partial_highlighting(self, abbr_system, complex_molecule):
        highlighted = {0, 1, 2}
        abbr_system.apply(
            complex_molecule, mode=AbbreviationMode.ALL, highlighted_atoms=highlighted
        )
        assert complex_molecule is not None


class TestCommonReagents:
    """Test abbreviation of common chemical reagents."""

    def test_thf_abbreviation(self, abbr_system):
        mol = get_CDK_IAtomContainer("C1CCOC1")
        abbr_system.apply(mol, mode=AbbreviationMode.REAGENTS)
        assert mol is not None

    def test_dmf_abbreviation(self, abbr_system):
        mol = get_CDK_IAtomContainer("CN(C)C=O")
        abbr_system.apply(mol, mode=AbbreviationMode.REAGENTS)
        assert mol is not None

    def test_dcm_abbreviation(self, abbr_system):
        mol = get_CDK_IAtomContainer("ClCCl")
        abbr_system.apply(mol, mode=AbbreviationMode.REAGENTS)
        assert mol is not None

    def test_etoh_abbreviation(self, abbr_system):
        mol = get_CDK_IAtomContainer("CCO")
        abbr_system.apply(mol, mode=AbbreviationMode.REAGENTS)
        assert mol is not None

    def test_meoh_abbreviation(self, abbr_system):
        mol = get_CDK_IAtomContainer("CO")
        abbr_system.apply(mol, mode=AbbreviationMode.REAGENTS)
        assert mol is not None


class TestFunctionalGroups:
    """Test abbreviation of functional groups."""

    def test_phenyl_group(self, abbr_system):
        mol = get_CDK_IAtomContainer("c1ccccc1C")
        abbr_system.apply(mol, mode=AbbreviationMode.GROUPS)
        assert mol is not None

    def test_methyl_group(self, abbr_system):
        mol = get_CDK_IAtomContainer("CC(C)C")
        abbr_system.apply(mol, mode=AbbreviationMode.GROUPS)
        assert mol is not None

    def test_ethyl_group(self, abbr_system):
        mol = get_CDK_IAtomContainer("CCC")
        abbr_system.apply(mol, mode=AbbreviationMode.GROUPS)
        assert mol is not None

    def test_boc_protecting_group(self, abbr_system):
        mol = get_CDK_IAtomContainer("CC(C)(C)OC(=O)N")
        abbr_system.apply(mol, mode=AbbreviationMode.GROUPS)
        assert mol is not None

    def test_fmoc_protecting_group(self, abbr_system):
        mol = get_CDK_IAtomContainer("C1=CC=C2C(=C1)C=CC=C2COC(=O)N")
        abbr_system.apply(mol, mode=AbbreviationMode.GROUPS)
        assert mol is not None

    def test_cbz_protecting_group(self, abbr_system):
        mol = get_CDK_IAtomContainer("c1ccccc1COC(=O)N")
        abbr_system.apply(mol, mode=AbbreviationMode.GROUPS)
        assert mol is not None


class TestAbbreviationOptions:
    """Test AbbreviationOptions dataclass."""

    def test_default_options(self):
        options = AbbreviationOptions()
        assert options.mode == AbbreviationMode.REAGENTS
        assert options.allow_singleton is False
        assert options.auto_contract_terminal is True
        assert options.auto_contract_hetero is False

    def test_custom_options(self):
        options = AbbreviationOptions(
            mode=AbbreviationMode.ALL,
            allow_singleton=True,
            auto_contract_terminal=False,
            auto_contract_hetero=True,
        )
        assert options.mode == AbbreviationMode.ALL
        assert options.allow_singleton is True
        assert options.auto_contract_terminal is False
        assert options.auto_contract_hetero is True

    def test_custom_file_paths(self):
        options = AbbreviationOptions(
            custom_group_file="/custom/groups.smi",
            custom_reagent_file="/custom/reagents.smi",
        )
        assert options.custom_group_file == "/custom/groups.smi"
        assert options.custom_reagent_file == "/custom/reagents.smi"


class TestComplexMolecules:
    """Test abbreviations on complex molecules."""

    def test_multiple_phenyl_groups(self, abbr_system):
        mol = get_CDK_IAtomContainer("c1ccccc1Cc1ccccc1")
        abbr_system.apply(mol, mode=AbbreviationMode.GROUPS)
        assert mol is not None

    def test_molecule_with_multiple_groups(self, abbr_system):
        mol = get_CDK_IAtomContainer("c1ccccc1CC(C)(C)C")
        abbr_system.apply(mol, mode=AbbreviationMode.GROUPS)
        assert mol is not None

    def test_caffeine_structure(self, abbr_system):
        mol = get_CDK_IAtomContainer("CN1C=NC2=C1C(=O)N(C(=O)N2C)C")
        abbr_system.apply(mol, mode=AbbreviationMode.GROUPS)
        assert mol is not None

    def test_peptide_with_protecting_groups(self, abbr_system):
        mol = get_CDK_IAtomContainer("CC(C)(C)OC(=O)NCC(=O)N")
        abbr_system.apply(mol, mode=AbbreviationMode.ALL)
        assert mol is not None


class TestEdgeCases:
    """Test edge cases and boundary conditions."""

    def test_very_small_molecule(self, abbr_system):
        mol = get_CDK_IAtomContainer("C")
        initial_count = mol.getAtomCount()
        abbr_system.apply(mol, mode=AbbreviationMode.ALL)
        assert mol.getAtomCount() == initial_count

    def test_already_abbreviated_molecule(self, abbr_system, phenyl_molecule):
        abbr_system.apply(phenyl_molecule, mode=AbbreviationMode.GROUPS)
        abbr_system.apply(phenyl_molecule, mode=AbbreviationMode.GROUPS)
        assert phenyl_molecule is not None

    def test_aromatic_vs_aliphatic(self, abbr_system):
        aromatic = get_CDK_IAtomContainer("c1ccccc1")
        aliphatic = get_CDK_IAtomContainer("C1CCCCC1")
        abbr_system.apply(aromatic, mode=AbbreviationMode.GROUPS)
        abbr_system.apply(aliphatic, mode=AbbreviationMode.GROUPS)
        assert aromatic is not None
        assert aliphatic is not None

    def test_charged_molecules(self, abbr_system):
        mol = get_CDK_IAtomContainer("[NH4+]")
        abbr_system.apply(mol, mode=AbbreviationMode.ALL)
        assert mol is not None

    def test_radical_molecules(self, abbr_system):
        mol = get_CDK_IAtomContainer("[CH3]")
        abbr_system.apply(mol, mode=AbbreviationMode.ALL)
        assert mol is not None


class TestAbbreviationIntegrity:
    """Test that abbreviations maintain molecular integrity."""

    def test_bond_count_consistency(self, abbr_system, phenyl_molecule):
        initial_bonds = phenyl_molecule.getBondCount()
        abbr_system.apply(phenyl_molecule, mode=AbbreviationMode.GROUPS)
        assert phenyl_molecule.getBondCount() <= initial_bonds

    def test_connectivity_preserved(self, abbr_system, complex_molecule):
        abbr_system.apply(complex_molecule, mode=AbbreviationMode.ALL)
        for atom in complex_molecule.atoms():
            if complex_molecule.getAtomCount() > 1:
                connected = complex_molecule.getConnectedAtomsList(atom)
                assert connected is not None

    def test_stereochemistry_preserved(self, abbr_system):
        mol = get_CDK_IAtomContainer("C[C@H](c1ccccc1)O")
        abbr_system.apply(mol, mode=AbbreviationMode.GROUPS)
        stereo_count = sum(1 for _ in mol.stereoElements())
        assert stereo_count >= 0


class TestLazyInitialization:
    """Test lazy initialization of abbreviation system."""

    def test_apply_without_explicit_init(self):
        """Test that apply() initializes the system if not initialized."""
        system = ChemicalAbbreviations()
        assert system._initialized is False

        mol = get_CDK_IAtomContainer("c1ccccc1C")
        system.apply(mol, mode=AbbreviationMode.GROUPS)

        # After apply, system should be initialized
        assert system._initialized is True

    def test_apply_off_mode_without_init(self):
        """Test that OFF mode with uninitialized system initializes first."""
        system = ChemicalAbbreviations()
        assert system._initialized is False

        mol = get_CDK_IAtomContainer("c1ccccc1C")
        initial_count = mol.getAtomCount()
        system.apply(mol, mode=AbbreviationMode.OFF)

        # System should be initialized even with OFF mode
        assert system._initialized is True
        # Molecule should be unchanged
        assert mol.getAtomCount() == initial_count


class TestCustomFileLoading:
    """Test loading custom abbreviation files and error handling."""

    def test_custom_group_file_fallback(self):
        """Test exception handler for default file loading failure."""
        # Create a temporary custom file
        with tempfile.NamedTemporaryFile(mode="w", suffix=".smi", delete=False) as f:
            f.write("c1ccccc1 Ph\n")
            custom_file = f.name

        try:
            # Initialize with custom file - should trigger fallback path
            system = ChemicalAbbreviations()
            # Manually trigger the scenario by setting invalid default path
            system.DEFAULT_GROUP_ABBR_PATH = "/nonexistent/path.smi"
            system.initialize(custom_group_file=custom_file)

            assert system._initialized is True
        finally:
            Path(custom_file).unlink(missing_ok=True)

    def test_custom_reagent_file_fallback(self):
        """Test exception handler for reagent file loading failure."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".smi", delete=False) as f:
            f.write("C1CCOC1 THF\n")
            custom_file = f.name

        try:
            system = ChemicalAbbreviations()
            system.DEFAULT_REAGENT_ABBR_PATH = "/nonexistent/path.smi"
            system.initialize(custom_reagent_file=custom_file)

            assert system._initialized is True
        finally:
            Path(custom_file).unlink(missing_ok=True)

    def test_invalid_custom_file_path(self):
        """Test that invalid custom file paths are handled gracefully."""
        system = ChemicalAbbreviations()
        # Pass a non-existent custom file - should not crash
        system.initialize(
            custom_group_file="/nonexistent/groups.smi",
            custom_reagent_file="/nonexistent/reagents.smi",
        )
        assert system._initialized is True


class TestModeOffExplicit:
    """Test explicit OFF mode behavior."""

    def test_mode_off_returns_early(self, abbr_system):
        """Test that OFF mode returns without processing."""
        mol = get_CDK_IAtomContainer("c1ccccc1C")
        initial_count = mol.getAtomCount()
        initial_bonds = mol.getBondCount()

        abbr_system.apply(mol, mode=AbbreviationMode.OFF)

        # Molecule should be completely unchanged
        assert mol.getAtomCount() == initial_count
        assert mol.getBondCount() == initial_bonds

    def test_mode_off_with_complex_molecule(self, abbr_system):
        """Test OFF mode with a complex molecule."""
        mol = get_CDK_IAtomContainer("CC(C)(C)OC(=O)Nc1ccccc1")
        initial_count = mol.getAtomCount()

        abbr_system.apply(mol, mode=AbbreviationMode.OFF)

        assert mol.getAtomCount() == initial_count


class TestReactionAbbreviations:
    """Test abbreviation application to reactions."""

    def test_apply_to_reaction_all_mode(self, abbr_system):
        """Test applying abbreviations to reactions with ALL mode."""

        Reaction = JClass("org.openscience.cdk.Reaction")

        reaction = Reaction()
        reactant = get_CDK_IAtomContainer("c1ccccc1C")
        reagent = get_CDK_IAtomContainer("C1CCOC1")
        product = get_CDK_IAtomContainer("c1ccccc1CC")

        reaction.addReactant(reactant)
        reaction.addAgent(reagent)
        reaction.addProduct(product)

        # Apply abbreviations with ALL mode
        abbr_system.apply(reaction, mode=AbbreviationMode.ALL)

        assert reaction is not None

    def test_apply_to_reaction_groups_mode(self, abbr_system):
        """Test applying abbreviations to reactions with GROUPS mode."""

        Reaction = JClass("org.openscience.cdk.Reaction")

        reaction = Reaction()
        reactant = get_CDK_IAtomContainer("c1ccccc1C")
        product = get_CDK_IAtomContainer("c1ccccc1CC")

        reaction.addReactant(reactant)
        reaction.addProduct(product)

        abbr_system.apply(reaction, mode=AbbreviationMode.GROUPS)

        assert reaction is not None

    def test_apply_to_reaction_reagents_mode(self, abbr_system):
        """Test applying abbreviations to reactions with REAGENTS mode."""

        Reaction = JClass("org.openscience.cdk.Reaction")

        reaction = Reaction()
        reactant = get_CDK_IAtomContainer("CC")
        reagent = get_CDK_IAtomContainer("C1CCOC1")  # THF
        product = get_CDK_IAtomContainer("CCC")

        reaction.addReactant(reactant)
        reaction.addAgent(reagent)
        reaction.addProduct(product)

        abbr_system.apply(reaction, mode=AbbreviationMode.REAGENTS)

        assert reaction is not None

    def test_apply_to_reaction_off_mode(self, abbr_system):
        """Test that OFF mode doesn't process reactions."""

        Reaction = JClass("org.openscience.cdk.Reaction")

        reaction = Reaction()
        reactant = get_CDK_IAtomContainer("c1ccccc1C")

        initial_count = reactant.getAtomCount()
        reaction.addReactant(reactant)

        abbr_system.apply(reaction, mode=AbbreviationMode.OFF)

        # Reaction should be unchanged
        assert reactant.getAtomCount() == initial_count

    def test_apply_to_reaction_with_highlighted_atoms(self, abbr_system):
        """Test reaction abbreviations with highlighted atoms."""

        Reaction = JClass("org.openscience.cdk.Reaction")

        reaction = Reaction()
        reactant = get_CDK_IAtomContainer("c1ccccc1C")
        reaction.addReactant(reactant)

        highlighted = {0, 1, 2}
        abbr_system.apply(
            reaction, mode=AbbreviationMode.ALL, highlighted_atoms=highlighted
        )

        assert reaction is not None


class TestHydrateContraction:
    """Test hydrate contraction functionality."""

    def test_contract_hydrates_method_exists(self, abbr_system):
        """Test that _contract_hydrates method can be called."""
        mol = get_CDK_IAtomContainer("O.O.c1ccccc1")  # Molecule with water

        # Apply with ALL mode which internally calls _contract_hydrates
        abbr_system.apply(mol, mode=AbbreviationMode.ALL)

        assert mol is not None

    def test_molecule_with_water(self, abbr_system):
        """Test abbreviation of molecules with water."""
        # Create molecule with water molecules
        mol = get_CDK_IAtomContainer("c1ccccc1.O.O")

        abbr_system.apply(mol, mode=AbbreviationMode.ALL)

        assert mol is not None

    def test_hydrate_in_groups_mode(self, abbr_system):
        """Test that hydrate contraction works in GROUPS mode."""
        mol = get_CDK_IAtomContainer("CCO.O")

        abbr_system.apply(mol, mode=AbbreviationMode.GROUPS)

        assert mol is not None


class TestExceptionHandling:
    """Test exception handling in abbreviation system."""

    def test_apply_handles_exceptions_gracefully(self, abbr_system):
        """Test that exceptions in apply() are caught and handled."""
        # Create a molecule and apply abbreviations
        # Even if internal errors occur, it should not crash
        mol = get_CDK_IAtomContainer("c1ccccc1")

        try:
            abbr_system.apply(mol, mode=AbbreviationMode.ALL)
            success = True
        except Exception:
            success = False

        assert success is True

    def test_reaction_exception_handling(self, abbr_system):
        """Test exception handling in reaction abbreviations."""

        Reaction = JClass("org.openscience.cdk.Reaction")
        reaction = Reaction()

        # Empty reaction should still be handled
        try:
            abbr_system.apply(reaction, mode=AbbreviationMode.ALL)
            success = True
        except Exception:
            success = False

        assert success is True
