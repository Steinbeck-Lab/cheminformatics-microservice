/**
 * Navigation Registry -- single source of truth for all pages, tools, and example molecules.
 *
 * Used by CommandPalette (Cmd+K search) and Breadcrumbs (section > tool wayfinding).
 */
import type { ComponentType } from "react";
import {
  Camera,
  FlaskConical,
  House,
  LineChart,
  RefreshCw,
  SlidersHorizontal,
  Users,
  // Chem tool icons
  Copy,
  Code,
  LayoutList,
  Zap,
  BarChart3,
  Fingerprint,
  Globe,
  CheckCircle,
  Filter,
  Database,
  Tag,
  // Convert tool icons
  LayoutGrid,
  Box,
  // Depict tool icons
  Search,
  Pencil,
  // Tools tool icons
  Puzzle,
  FileText,
  ArrowLeftRight,
} from "lucide-react";

// --- Types ---

export interface NavEntry {
  id: string;
  name: string;
  path: string;
  section?: string;
  sectionPath?: string;
  icon?: ComponentType<{ className?: string }>;
  keywords?: string[];
}

export interface MoleculeEntry {
  name: string;
  smiles: string;
  description?: string;
}

// --- Pages ---

export const pages: NavEntry[] = [
  { id: "home", name: "Home", path: "/home", icon: House },
  { id: "chem", name: "Chemical Analysis", path: "/chem", icon: FlaskConical },
  { id: "convert", name: "Format Conversion", path: "/convert", icon: RefreshCw },
  { id: "depict", name: "Depiction", path: "/depict", icon: LineChart },
  { id: "tools", name: "Tools", path: "/tools", icon: SlidersHorizontal },
  { id: "ocsr", name: "OCSR", path: "/ocsr", icon: Camera },
  { id: "about", name: "About", path: "/about", icon: Users },
];

// --- Tools ---

export const tools: NavEntry[] = [
  // Chemical Analysis (14 tools)
  {
    id: "stereoisomers",
    name: "Stereoisomers",
    path: "/chem/stereoisomers",
    section: "Chemical Analysis",
    sectionPath: "/chem",
    icon: Copy,
    keywords: ["stereo", "isomers", "chirality", "enantiomers"],
  },
  {
    id: "hosecodes",
    name: "HOSE Codes",
    path: "/chem/hosecodes",
    section: "Chemical Analysis",
    sectionPath: "/chem",
    icon: Code,
    keywords: ["hose", "codes", "NMR", "environment"],
  },
  {
    id: "standardize",
    name: "Standardize",
    path: "/chem/standardize",
    section: "Chemical Analysis",
    sectionPath: "/chem",
    icon: RefreshCw,
    keywords: ["normalize", "canonical", "clean"],
  },
  {
    id: "functional-groups",
    name: "Functional Groups",
    path: "/chem/functional-groups",
    section: "Chemical Analysis",
    sectionPath: "/chem",
    icon: LayoutList,
    keywords: ["ertl", "functional", "groups", "moieties"],
  },
  {
    id: "tautomer",
    name: "Tautomers",
    path: "/chem/tautomer",
    section: "Chemical Analysis",
    sectionPath: "/chem",
    icon: RefreshCw,
    keywords: ["tautomer", "tautomers", "proton shift"],
  },
  {
    id: "fixradicals",
    name: "Fix Radicals",
    path: "/chem/fixradicals",
    section: "Chemical Analysis",
    sectionPath: "/chem",
    icon: Zap,
    keywords: ["radical", "electrons", "fix", "repair"],
  },
  {
    id: "descriptors",
    name: "Descriptors",
    path: "/chem/descriptors",
    section: "Chemical Analysis",
    sectionPath: "/chem",
    icon: BarChart3,
    keywords: ["molecular weight", "logP", "properties", "descriptors", "TPSA"],
  },
  {
    id: "nplikeness",
    name: "NP-likeness",
    path: "/chem/nplikeness",
    section: "Chemical Analysis",
    sectionPath: "/chem",
    icon: FlaskConical,
    keywords: ["natural product", "NP score", "likeness"],
  },
  {
    id: "similarity",
    name: "Similarity",
    path: "/chem/similarity",
    section: "Chemical Analysis",
    sectionPath: "/chem",
    icon: Fingerprint,
    keywords: ["tanimoto", "fingerprint", "similarity", "compare"],
  },
  {
    id: "structure-finder",
    name: "Structure Finder",
    path: "/chem/structure-finder",
    section: "Chemical Analysis",
    sectionPath: "/chem",
    icon: Globe,
    keywords: ["pubchem", "lookup", "search", "identifier", "name"],
  },
  {
    id: "structureerror",
    name: "Check Structure",
    path: "/chem/structureerror",
    section: "Chemical Analysis",
    sectionPath: "/chem",
    icon: CheckCircle,
    keywords: ["validate", "error", "check", "verify"],
  },
  {
    id: "filters",
    name: "All Filters",
    path: "/chem/filters",
    section: "Chemical Analysis",
    sectionPath: "/chem",
    icon: Filter,
    keywords: ["lipinski", "veber", "PAINS", "drug-like", "filter"],
  },
  {
    id: "coconut",
    name: "COCONUT Preprocessing",
    path: "/chem/coconut",
    section: "Chemical Analysis",
    sectionPath: "/chem",
    icon: Database,
    keywords: ["coconut", "preprocessing", "natural products", "database"],
  },
  {
    id: "classyfire",
    name: "ClassyFire",
    path: "/chem/classyfire",
    section: "Chemical Analysis",
    sectionPath: "/chem",
    icon: Tag,
    keywords: ["classification", "taxonomy", "classyfire", "chemical class"],
  },

  // Format Conversion (3 tools)
  {
    id: "format-conversion",
    name: "Format Conversion",
    path: "/convert/format-conversion",
    section: "Format Conversion",
    sectionPath: "/convert",
    icon: Copy,
    keywords: ["SMILES", "InChI", "mol", "SDF", "convert", "format"],
  },
  {
    id: "2d-coordinates",
    name: "2D Coordinates",
    path: "/convert/2d-coordinates",
    section: "Format Conversion",
    sectionPath: "/convert",
    icon: LayoutGrid,
    keywords: ["2D", "coordinates", "layout", "generate"],
  },
  {
    id: "3d-coordinates",
    name: "3D Coordinates",
    path: "/convert/3d-coordinates",
    section: "Format Conversion",
    sectionPath: "/convert",
    icon: Box,
    keywords: ["3D", "coordinates", "conformer", "generate"],
  },

  // Depiction (4 tools)
  {
    id: "structureexplorer",
    name: "Structure Explorer",
    path: "/depict/structureexplorer",
    section: "Depiction",
    sectionPath: "/depict",
    icon: Search,
    keywords: ["visualize", "explore", "identifier", "name", "search"],
  },
  {
    id: "2ddepiction",
    name: "2D Depiction",
    path: "/depict/2ddepiction",
    section: "Depiction",
    sectionPath: "/depict",
    icon: LayoutGrid,
    keywords: ["2D", "depiction", "image", "SVG", "render"],
  },
  {
    id: "3ddepiction",
    name: "3D Depiction",
    path: "/depict/3ddepiction",
    section: "Depiction",
    sectionPath: "/depict",
    icon: Box,
    keywords: ["3D", "visualization", "interactive", "3Dmol"],
  },
  {
    id: "structuredraw",
    name: "Draw a Structure",
    path: "/depict/structuredraw",
    section: "Depiction",
    sectionPath: "/depict",
    icon: Pencil,
    keywords: ["draw", "editor", "sketch", "JSME", "ketcher"],
  },

  // Tools (4 tools)
  {
    id: "sugardetection",
    name: "Sugar Detection",
    path: "/tools/sugardetection",
    section: "Tools",
    sectionPath: "/tools",
    icon: Box,
    keywords: ["sugar", "removal", "deglycosylation", "moiety"],
  },
  {
    id: "structuregeneration",
    name: "Structure Generation",
    path: "/tools/structuregeneration",
    section: "Tools",
    sectionPath: "/tools",
    icon: Puzzle,
    keywords: ["generate", "structure", "virtual screening", "enumeration"],
  },
  {
    id: "inchiconverter",
    name: "InChI Converter",
    path: "/tools/inchiconverter",
    section: "Tools",
    sectionPath: "/tools",
    icon: FileText,
    keywords: ["InChI", "InChIKey", "converter", "notation"],
  },
  {
    id: "rinchiconverter",
    name: "RInChI Converter",
    path: "/tools/rinchiconverter",
    section: "Tools",
    sectionPath: "/tools",
    icon: ArrowLeftRight,
    keywords: ["RInChI", "reaction", "converter", "transformation"],
  },
];

// --- Example Molecules ---

export const exampleMolecules: MoleculeEntry[] = [
  {
    name: "Caffeine",
    smiles: "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
    description: "Stimulant found in coffee",
  },
  {
    name: "Aspirin",
    smiles: "CC(=O)OC1=CC=CC=C1C(=O)O",
    description: "Analgesic & anti-inflammatory",
  },
  {
    name: "Sucrose",
    smiles: "C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O[C@]2([C@H]([C@@H]([C@H](O2)CO)O)O)CO)O)O)O)O",
    description: "Common sugar",
  },
  {
    name: "Cholesterol",
    smiles: "C[C@H](CCCC(C)C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC=C4[C@@]3(CC[C@@H](C4)O)C)C",
    description: "Steroid in cell membranes",
  },
  {
    name: "Paracetamol",
    smiles: "CC(=O)NC1=CC=C(C=C1)O",
    description: "Pain reliever & fever reducer",
  },
  {
    name: "Ibuprofen",
    smiles: "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O",
    description: "Anti-inflammatory drug",
  },
];

// --- Helpers ---

/**
 * Look up a tool by its route pathname and optional toolId param.
 * For example, getToolByPath("/chem", "descriptors") returns the Descriptors entry.
 */
export function getToolByPath(pathname: string, toolId?: string): NavEntry | undefined {
  if (!toolId) return undefined;
  // Find a tool whose sectionPath matches the base route and whose id matches the toolId
  return tools.find((t) => t.sectionPath === pathname.replace(/\/[^/]+$/, "") || t.id === toolId);
}

/**
 * Look up the section (page) for a given pathname.
 * For example, getSectionByPath("/chem/descriptors") returns { name: "Chemical Analysis", path: "/chem" }.
 */
export function getSectionByPath(pathname: string): { name: string; path: string } | undefined {
  // Match on the first path segment
  const basePath = "/" + pathname.split("/").filter(Boolean)[0];
  const page = pages.find((p) => p.path === basePath);
  if (!page) return undefined;
  return { name: page.name, path: page.path };
}
