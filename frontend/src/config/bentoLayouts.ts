/**
 * Bento grid layout sizing configuration per page.
 *
 * Each page's tools are assigned a size category (hero/medium/compact)
 * that determines their visual prominence in the bento grid.
 * Hero cards span 2x2, medium cards span 1x2 or 2x1, compact cards are 1x1.
 */

/** Size categories for bento grid cards. */
export type BentoSize = "hero" | "medium" | "compact";

/** Layout entry for a single tool in the bento grid. */
export interface BentoLayoutEntry {
  id: string;
  size: BentoSize;
}

/** Bento layout configuration keyed by page. */
export const BENTO_LAYOUTS: Record<string, BentoLayoutEntry[]> = {
  chem: [
    // Hero -- most-used tools
    { id: "descriptors", size: "hero" },
    { id: "structure-finder", size: "hero" },
    // Medium -- secondary tools
    { id: "stereoisomers", size: "medium" },
    { id: "standardize", size: "medium" },
    { id: "similarity", size: "medium" },
    { id: "classyfire", size: "medium" },
    // Compact -- minor tools
    { id: "hose-codes", size: "compact" },
    { id: "functional-groups", size: "compact" },
    { id: "tautomers", size: "compact" },
    { id: "fix-radicals", size: "compact" },
    { id: "np-likeness", size: "compact" },
    { id: "check-structure", size: "compact" },
    { id: "all-filters", size: "compact" },
    { id: "coconut-preprocessing", size: "compact" },
  ],
  tools: [
    // Hero -- primary tools
    { id: "sugardetection", size: "hero" },
    { id: "structuregeneration", size: "hero" },
    // Medium -- secondary tools
    { id: "inchiconverter", size: "medium" },
    { id: "rinchiconverter", size: "medium" },
  ],
  convert: [
    // Hero -- primary tool
    { id: "format-conversion", size: "hero" },
    // Medium -- secondary tools
    { id: "2d-coordinates", size: "medium" },
    { id: "3d-coordinates", size: "medium" },
  ],
  depict: [
    // Hero -- primary tool
    { id: "structure-explorer", size: "hero" },
    // Medium -- secondary tools
    { id: "2d-depiction", size: "medium" },
    { id: "3d-depiction", size: "medium" },
    // Compact -- minor tools
    { id: "draw-structure", size: "compact" },
  ],
  home: [
    // HomePage feature cards -- all medium (editorial grid, no hero distinction)
    { id: "chemical-analysis", size: "medium" },
    { id: "format-conversion", size: "medium" },
    { id: "structure-depiction", size: "medium" },
    { id: "ocsr", size: "medium" },
    { id: "tools", size: "medium" },
  ],
} as const;

/**
 * Get the bento size for a specific tool on a specific page.
 *
 * @param pageKey - The page identifier (e.g., 'chem', 'tools')
 * @param toolId - The tool identifier
 * @returns The bento size, or 'compact' as the default fallback
 */
export function getBentoSize(pageKey: string, toolId: string): BentoSize {
  const layout = BENTO_LAYOUTS[pageKey];
  if (!layout) return "compact";
  const entry = layout.find((item) => item.id === toolId);
  return entry?.size ?? "compact";
}
