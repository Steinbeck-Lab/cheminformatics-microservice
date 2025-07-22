# Functional Group Highlighting Feature

## Overview

The Functional Group Detection (Ertl) tool now includes interactive highlighting capabilities that allow users to visualize detected functional groups directly on the molecular structure.

## Features

### Interactive Functional Group List
- Click on any functional group in the results list to highlight it in the molecular structure
- Visual feedback shows which group is currently selected
- Hover effects provide clear interaction cues
- Atom count information displayed for each functional group

### Visual Highlighting
- Functional groups are highlighted using SMARTS pattern matching
- Common functional groups have predefined highlighting patterns
- Blue highlighting overlay shows the selected functional group atoms
- Clear highlighting button to remove selection

### Enhanced Data Structure
The backend now returns structured functional group data including:
- `atomIds`: Array of atom indices involved in the functional group
- `atoms`: Atom symbols in the functional group
- `type`: Functional group type identifier
- `description`: Human-readable description

## How to Use

1. **Enter a SMILES string** in the input field
2. **Click "Detect Functional Groups"** to analyze the molecule
3. **View results** in the split-panel layout:
   - Left panel: Interactive molecular structure
   - Right panel: List of detected functional groups
4. **Click on any functional group** in the list to highlight it in the structure
5. **Click "Clear Highlighting"** or click the same group again to remove highlighting

## Technical Implementation

### Backend Changes
- Modified `get_ertl_functional_groups()` in `rdkit_wrapper.py` to return structured data
- Each functional group now includes atom IDs for precise highlighting
- Maintains backward compatibility with existing functionality

### Frontend Changes
- New `HighlightedMoleculeCard` component for enhanced visualization
- Integration with 2D depiction API's highlight parameter
- Interactive functional group selection interface
- SMARTS pattern generation from functional group data

### Highlighting Patterns
The system uses common SMARTS patterns for highlighting:
- Hydroxyl groups: `[OH]`
- Carbonyl groups: `[C]=[O]` 
- Carboxyl groups: `[C](=[O])[OH]`
- Amino groups: `[NH2]`
- Aromatic rings: `c1ccccc1`
- And many more...

## Benefits

1. **Educational**: Helps users understand molecular structure and functional groups
2. **Analytical**: Visual confirmation of detected functional groups
3. **Interactive**: Engaging user experience with immediate visual feedback
4. **Accurate**: Uses established SMARTS pattern matching for precise highlighting

## Example Usage

For caffeine (SMILES: `CN1C=NC2=C1C(=O)N(C(=O)N2C)C`):
- Detects 6 functional groups including carbonyls and nitrogen-containing groups
- Each group can be individually highlighted
- Shows atom count and type information
- Provides clear visual distinction of different functional groups

## Future Enhancements

- Custom SMARTS pattern input for highlighting
- Multiple simultaneous group highlighting
- Color-coded highlighting for different group types
- Export highlighted structures
- Integration with other analysis tools
