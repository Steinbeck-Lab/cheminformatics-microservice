import {
  BarChart3,
  CheckCircle,
  Code,
  Copy,
  Database,
  Filter,
  Fingerprint,
  FlaskConical,
  Globe,
  LayoutList,
  RefreshCw,
  Tag,
  Zap,
} from "lucide-react";
import StereoisomersView from "../components/chem/StereoisomersView";
import HOSECodeView from "../components/chem/HOSECodeView";
import StandardizeView from "../components/chem/StandardizeView";
import ErtlFunctionalGroupView from "../components/chem/ErtlFunctionalGroupView";
import StandardizedTautomerView from "../components/chem/StandardizedTautomerView";
import FixRadicalsView from "../components/chem/FixRadicalsView";
import PubChemLookupView from "../components/chem/PubChemLookupView";
import DescriptorsView from "../components/chem/DescriptorsView";
import NPlikenessView from "../components/chem/NPlikenessView";
import TanimotoView from "../components/chem/TanimotoView";
import StructureErrorView from "../components/chem/StructureErrorView";
import AllFiltersView from "../components/chem/AllFiltersView";
import CoconutPreProcessingView from "../components/chem/CoconutPreProcessingView";
import ClassyfireView from "../components/chem/ClassyfireView";
import { SectionPage, type SectionTool } from "../components/common/SectionPage";

// Grouped in sidebar order: structure, search, analysis, comparison, validation, advanced.
const tools: SectionTool[] = [
  {
    id: "stereoisomers",
    name: "Stereoisomers",
    description: "Generate all possible stereoisomers.",
    category: "Structure tools",
    icon: Copy,
    component: StereoisomersView,
  },
  {
    id: "hosecodes",
    name: "HOSE Codes",
    description: "Generate HOSE codes.",
    category: "Structure tools",
    icon: Code,
    component: HOSECodeView,
  },
  {
    id: "standardize",
    name: "Standardize",
    description: "Standardize structures with the ChEMBL structure pipeline.",
    category: "Structure tools",
    icon: RefreshCw,
    component: StandardizeView,
  },
  {
    id: "functional-groups",
    name: "Functional Groups",
    description: "Identify functional groups (Ertl).",
    category: "Structure tools",
    icon: LayoutList,
    component: ErtlFunctionalGroupView,
  },
  {
    id: "tautomer",
    name: "Tautomers",
    description: "Generate standardized tautomers.",
    category: "Structure tools",
    icon: RefreshCw,
    component: StandardizedTautomerView,
  },
  {
    id: "fixradicals",
    name: "Fix Radicals",
    description: "Fix radical electrons in molecules.",
    category: "Structure tools",
    icon: Zap,
    component: FixRadicalsView,
  },
  {
    id: "structure-finder",
    name: "Structure Finder",
    description: "Find chemical structures by name or identifier.",
    category: "Search and retrieval",
    icon: Globe,
    component: PubChemLookupView,
  },
  {
    id: "descriptors",
    name: "Descriptors",
    description: "Calculate molecular descriptors.",
    category: "Property analysis",
    icon: BarChart3,
    component: DescriptorsView,
  },
  {
    id: "nplikeness",
    name: "NP-likeness",
    description: "Calculate the natural-product likeness score.",
    category: "Property analysis",
    icon: FlaskConical,
    component: NPlikenessView,
  },
  {
    id: "similarity",
    name: "Similarity",
    description: "Compare structures with Tanimoto similarity.",
    category: "Comparison",
    icon: Fingerprint,
    component: TanimotoView,
  },
  {
    id: "structureerror",
    name: "Check Structure",
    description: "Validate chemical structures.",
    category: "Validation and filtering",
    icon: CheckCircle,
    component: StructureErrorView,
  },
  {
    id: "filters",
    name: "All Filters",
    description: "Apply multiple chemical filters.",
    category: "Validation and filtering",
    icon: Filter,
    component: AllFiltersView,
  },
  {
    id: "coconut",
    name: "COCONUT Preprocessing",
    description: "Prepare data for the COCONUT database.",
    category: "Advanced",
    icon: Database,
    component: CoconutPreProcessingView,
  },
  {
    id: "classyfire",
    name: "ClassyFire",
    description: "Chemical classification through the ClassyFire service.",
    category: "Advanced",
    icon: Tag,
    component: ClassyfireView,
  },
];

const ChemPage = () => (
  <SectionPage label="Analysis tools" basePath="/chem" tools={tools} layout="sidebar" />
);

export default ChemPage;
