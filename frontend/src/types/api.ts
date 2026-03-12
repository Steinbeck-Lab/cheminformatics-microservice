/**
 * Shared API response type definitions for the Cheminformatics Microservice.
 */

// --- Chem Service Types ---

export interface StructureErrorResult {
  messages: string;
  standardized?: {
    smi: string;
  };
}

export interface DescriptorResult {
  [key: string]: string | number | boolean;
}

export interface MultipleDescriptorResult {
  [smiles: string]: DescriptorResult;
}

export interface FixRadicalsResult {
  fixed_smiles: string;
  radicals_detected: boolean;
  radicals_fixed: boolean;
}

export interface StandardizeResult {
  [key: string]: unknown;
}

export interface PubChemLookupResult {
  canonical_smiles: string;
  input: string;
  input_type: string;
  success: boolean;
}

export interface CoconutPreprocessingResult {
  [key: string]: unknown;
}

// --- Convert Service Types ---

export interface MultipleFormatsResult {
  [format: string]: string;
}

export interface CdxConversionResult {
  molblock: string;
}

// --- Depict Service Types ---

export interface DepictionOptions {
  toolkit?: string;
  width?: number;
  height?: number;
  rotate?: number;
  CIP?: boolean;
  unicolor?: boolean;
  highlight?: string;
  atomIds?: number[] | number[][] | null;
  showAtomNumbers?: boolean;
  hydrogen_display?: string;
}

export interface EnhancedDepictionOptions extends DepictionOptions {
  abbreviate?: string;
  dative?: string;
  multicenter?: string;
  annotate?: string;
  style?: string;
  donuts?: boolean;
  arrow?: string;
  alignrxnmap?: boolean;
  showtitle?: boolean;
  title?: string | null;
  bgcolor?: string | null;
  fgcolor?: string | null;
  zoom?: number;
  ratio?: number;
  flip?: boolean;
  anon?: boolean;
  smalim?: number;
  svgunits?: string;
  perceive_radicals?: boolean;
  apply_mdl_highlighting?: boolean;
  format?: string;
}

export interface BatchDepictionResult {
  smiles: string;
  svg: string;
}

export interface ParsedSmilesEntry {
  smiles: string;
  title: string;
}

// --- OCSR Service Types ---

export interface OCSRResult {
  [key: string]: unknown;
}

// --- Tools Service Types ---

export interface SugarRemovalOptions {
  gly_bond?: boolean;
  oxygen_atoms?: boolean;
  oxygen_atoms_threshold?: number;
  linear_sugars_in_rings?: boolean;
  linear_sugars_min_size?: number;
  linear_sugars_max_size?: number;
  linear_acidic_sugars?: boolean;
  spiro_sugars?: boolean;
  keto_sugars?: boolean;
  only_terminal?: boolean;
  preservation_mode?: number;
  preservation_threshold?: number;
  mark_attach_points?: boolean;
}

export interface ExtractionOptions extends SugarRemovalOptions {
  extract_circular_sugars?: boolean;
  extract_linear_sugars?: boolean;
  post_process_sugars?: boolean;
  limit_post_process_by_size?: boolean;
}

export interface ChemicalFilterOptions {
  pains?: boolean;
  lipinski?: boolean;
  veber?: boolean;
  reos?: boolean;
  ghose?: boolean;
  ruleofthree?: boolean;
  qedscore?: string;
  sascore?: string;
  nplikeness?: string;
  filterOperator?: string;
}

export interface ClassyFireResult {
  [key: string]: unknown;
}

// --- InChI Utils Types ---

export interface InChIVersionConfig {
  label: string;
  scriptSrc: string;
  moduleName: string;
  default?: boolean;
}

export interface InChIResult {
  [key: string]: unknown;
}

export interface InChIKeyResult {
  [key: string]: unknown;
}

// --- RInChI Utils Types ---

export interface RInChIVersionConfig {
  label: string;
  scriptSrc: string;
  moduleName: string;
  default?: boolean;
}

export interface RInChIResult {
  return_code: number;
  rinchi: string;
  rauxinfo: string;
  error: string;
}

export interface RInChIKeyResult {
  return_code: number;
  rinchikey: string;
  error: string;
}

export interface RInChIFileResult {
  return_code: number;
  fileText: string;
  error: string;
}
