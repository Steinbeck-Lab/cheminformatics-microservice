---
phase: 05-loading-states-ux
plan: 02
subsystem: frontend-ux
tags: [loading-states, error-handling, skeletons, toast, refactoring]
dependency_graph:
  requires: [05-01]
  provides: [glass-skeleton-integration, error-card-integration, loadingscreen-removal]
  affects: [all-tool-views, molecule-cards]
tech_stack:
  added: []
  patterns: [ToolSkeleton-per-component, GlassErrorCard-with-retry, getErrorMessage-domain-mapping, sonner-toast-error]
key_files:
  created: []
  modified:
    - frontend/src/components/chem/AllFiltersView.tsx
    - frontend/src/components/chem/ClassyfireView.tsx
    - frontend/src/components/chem/CoconutPreProcessingView.tsx
    - frontend/src/components/chem/DescriptorsView.tsx
    - frontend/src/components/chem/ErtlFunctionalGroupView.tsx
    - frontend/src/components/chem/HOSECodeView.tsx
    - frontend/src/components/chem/NPlikenessView.tsx
    - frontend/src/components/chem/PubChemLookupView.tsx
    - frontend/src/components/chem/StandardizeView.tsx
    - frontend/src/components/chem/StereoisomersView.tsx
    - frontend/src/components/chem/StructureErrorView.tsx
    - frontend/src/components/chem/TanimotoView.tsx
    - frontend/src/components/convert/FormatConversionView.tsx
    - frontend/src/components/convert/Mol2DView.tsx
    - frontend/src/components/convert/Mol3DView.tsx
    - frontend/src/components/depict/Depict2DMultiView.tsx
    - frontend/src/components/depict/StructureVisualizerView.tsx
    - frontend/src/components/tools/StructureGenView.tsx
    - frontend/src/components/tools/SugarRemovalView.tsx
    - frontend/src/components/common/MoleculeCard.tsx
    - frontend/src/components/common/HighlightedMoleculeCard.tsx
  deleted:
    - frontend/src/components/common/LoadingScreen.tsx
decisions:
  - "StandardizeView keeps domain-specific molblock error messages with getErrorMessage as fallback"
  - "StructureVisualizerView retry focuses identifier-input (not smiles-input) since it uses a text identifier"
  - "SugarRemovalView ToolSkeleton always shows during loading (no array-state guard) since results are object-based"
metrics:
  duration: 16min
  completed: 2026-03-13
---

# Phase 05 Plan 02: Tool View LoadingScreen Refactoring Summary

Replaced full-page LoadingScreen modal overlay with per-component ToolSkeleton glass skeletons, GlassErrorCard error display with retry, and domain-mapped getErrorMessage in all 21 tool views across 4 directories (chem/, convert/, depict/, tools/). Replaced alert() with sonner toast.error in MoleculeCard and HighlightedMoleculeCard. Deleted LoadingScreen.tsx.

## What Was Done

### Task 1a: Refactor 14 chem/ tool views (commit 08ba5b3)

Refactored all 14 files in `frontend/src/components/chem/`:
- Removed `import LoadingScreen` from each file
- Added imports for `ToolSkeleton`, `GlassErrorCard`, `EmptyState`, and `getErrorMessage`
- Replaced `{loading && <LoadingScreen text="..." />}` with `{loading && !result && <ToolSkeleton variant="..." />}`
- Replaced inline error `<div>` blocks with `<GlassErrorCard message={error} onRetry={...} />`
- Updated all catch blocks from template literals to `setError(getErrorMessage("chem", err))`
- Fixed ToolSkeleton conditions for array-based state vars (`.length === 0` instead of `!var`)
- StandardizeView kept domain-specific molblock error messages with getErrorMessage as fallback

ToolSkeleton variants matched per tool:
- `descriptors`: DescriptorsView, NPlikenessView, AllFiltersView, ClassyfireView
- `molecule`: StereoisomersView, StructureErrorView, ErtlFunctionalGroupView, HOSECodeView, PubChemLookupView, TanimotoView, StandardizeView, CoconutPreProcessingView

### Task 1b: Refactor 7 convert/depict/tools views (commit a846494)

Refactored all 7 files across convert/, depict/, and tools/:
- Same pattern: removed LoadingScreen, added feedback component imports
- Domain-mapped error messages: `getErrorMessage("convert", err)`, `getErrorMessage("depict", err)`, `getErrorMessage("tools", err)`
- FormatConversionView: variant="conversion"
- Mol2DView, Mol3DView: variant="conversion"
- Depict2DMultiView: variant="molecule", also updated secondary catch blocks (download errors)
- StructureVisualizerView: variant="molecule", retry focuses identifier-input
- StructureGenView: variant="molecule", retry focuses formula-input
- SugarRemovalView: variant="molecule"

### Task 2: Replace alert() with toast.error, delete LoadingScreen (commit 0024f39)

- MoleculeCard.tsx: added `import { toast } from "sonner"`, replaced `alert(...)` with `toast.error(...)`
- HighlightedMoleculeCard.tsx: same toast.error replacement
- Deleted `frontend/src/components/common/LoadingScreen.tsx` after verifying zero imports remain

## Verification

- TypeScript compilation: PASSED (zero errors after each task)
- Grep for `import LoadingScreen` across frontend/src: zero matches
- Grep for `alert(` in MoleculeCard and HighlightedMoleculeCard: zero matches
- LoadingScreen.tsx file: confirmed deleted
- GlassErrorCard present in all 21 tool views (import + usage)
- getErrorMessage present in all 21 tool views
- ToolSkeleton present in all 21 tool views

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] Fixed ToolSkeleton conditions for array-based state vars**
- **Found during:** Task 1a
- **Issue:** `{loading && !stereoisomers && <ToolSkeleton />}` never shows because `![]` is `false`
- **Fix:** Changed to `{loading && stereoisomers.length === 0 && <ToolSkeleton />}` for all array-based state vars
- **Files modified:** AllFiltersView, ErtlFunctionalGroupView, HOSECodeView, StereoisomersView
- **Commit:** 08ba5b3

**2. [Rule 2 - Missing functionality] StructureVisualizerView retry focuses identifier-input**
- **Found during:** Task 1b
- **Issue:** StructureVisualizerView uses a text identifier input (not SMILES), so retry should focus `#identifier-input`
- **Fix:** Used `document.getElementById("identifier-input")?.focus()` instead of `smiles-input`
- **Files modified:** StructureVisualizerView.tsx
- **Commit:** a846494

## Decisions Made

1. **StandardizeView error handling**: Kept domain-specific molblock error messages (cannot access local variable, invalid molblock format, parsing failure, server error) with `getErrorMessage("chem", err)` as the generic fallback
2. **StructureVisualizerView retry target**: Focuses `identifier-input` element since this view uses a generic text identifier, not SMILES input
3. **SugarRemovalView ToolSkeleton**: Always shows during loading without array guard since results are stored as an object, not array

## Commits

| Task | Commit | Description |
|------|--------|-------------|
| 1a | 08ba5b3 | Refactor 14 chem/ tool views to glass skeletons + error cards |
| 1b | a846494 | Refactor 7 convert/depict/tools views to glass skeletons + error cards |
| 2 | 0024f39 | Replace alert() with toast.error and delete LoadingScreen |
