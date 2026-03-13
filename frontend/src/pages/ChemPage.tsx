/**
 * ChemPage -- bento grid layout for 14 chemical analysis tools.
 *
 * Replaced sidebar-tab navigation with BentoGrid expand-in-place UX.
 * Each tool is a card that expands to show its full UI when clicked.
 */
import { useState, useEffect } from "react";
import { motion } from "motion/react";
import { useParams, useNavigate } from "react-router-dom";

// Import Views
import StereoisomersView from "../components/chem/StereoisomersView";
import DescriptorsView from "../components/chem/DescriptorsView";
import NPlikenessView from "../components/chem/NPlikenessView";
import TanimotoView from "../components/chem/TanimotoView";
import HOSECodeView from "../components/chem/HOSECodeView";
import StructureErrorView from "../components/chem/StructureErrorView";
import CoconutPreProcessingView from "../components/chem/CoconutPreProcessingView";
import StandardizeView from "../components/chem/StandardizeView";
import ClassyfireView from "../components/chem/ClassyfireView";
import ErtlFunctionalGroupView from "../components/chem/ErtlFunctionalGroupView";
import StandardizedTautomerView from "../components/chem/StandardizedTautomerView";
import AllFiltersView from "../components/chem/AllFiltersView";
import PubChemLookupView from "../components/chem/PubChemLookupView";
import FixRadicalsView from "../components/chem/FixRadicalsView";
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

import { BentoGrid, BentoCard } from "@/components/common/BentoGrid";
import { getBentoSize } from "@/config/bentoLayouts";
import { GradientMesh } from "@/components/common/GradientMesh";

// Define tab data with icons and components
const tabs = [
  // Structure
  {
    id: "stereoisomers",
    name: "Stereoisomers",
    description: "Generate all possible stereoisomers",
    icon: Copy,
    component: StereoisomersView,
  },
  {
    id: "hosecodes",
    name: "HOSE Codes",
    description: "Generate HOSE codes",
    icon: Code,
    component: HOSECodeView,
  },
  {
    id: "standardize",
    name: "Standardize",
    description: "Standardize structures",
    icon: RefreshCw,
    component: StandardizeView,
  },
  {
    id: "functional-groups",
    name: "Functional Groups",
    description: "Identify functional groups",
    icon: LayoutList,
    component: ErtlFunctionalGroupView,
  },
  {
    id: "tautomer",
    name: "Tautomers",
    description: "Generate standardized tautomers",
    icon: RefreshCw,
    component: StandardizedTautomerView,
  },
  {
    id: "fixradicals",
    name: "Fix Radicals",
    description: "Fix radical electrons in molecules",
    icon: Zap,
    component: FixRadicalsView,
  },
  // Analysis
  {
    id: "descriptors",
    name: "Descriptors",
    description: "Calculate molecular descriptors",
    icon: BarChart3,
    component: DescriptorsView,
  },
  {
    id: "nplikeness",
    name: "NP-likeness",
    description: "Calculate NP-likeness score",
    icon: FlaskConical,
    component: NPlikenessView,
  },
  // Comparison
  {
    id: "similarity",
    name: "Similarity",
    description: "Compare structures (Tanimoto)",
    icon: Fingerprint,
    component: TanimotoView,
  },
  // Search
  {
    id: "structure-finder",
    name: "Structure Finder",
    description: "Find chemical structures by name or identifier",
    icon: Globe,
    component: PubChemLookupView,
  },
  // Validation
  {
    id: "structureerror",
    name: "Check Structure",
    description: "Validate chemical structures",
    icon: CheckCircle,
    component: StructureErrorView,
  },
  {
    id: "filters",
    name: "All Filters",
    description: "Apply multiple chemical filters",
    icon: Filter,
    component: AllFiltersView,
  },
  // Advanced
  {
    id: "coconut",
    name: "COCONUT Preprocessing",
    description: "Prepare data for COCONUT DB",
    icon: Database,
    component: CoconutPreProcessingView,
  },
  {
    id: "classyfire",
    name: "ClassyFire",
    description: "Chemical classification",
    icon: Tag,
    component: ClassyfireView,
  },
];

// --- Animation Variants ---
const pageVariants = {
  hidden: { opacity: 0 },
  visible: { opacity: 1, transition: { duration: 0.5 } },
};

const ChemPage = () => {
  const { toolId } = useParams<{ toolId?: string }>();
  const navigate = useNavigate();
  const [expandedToolId, setExpandedToolId] = useState<string | null>(null);

  // Initialize expanded tool from URL parameter on mount
  useEffect(() => {
    if (toolId) {
      const matchingTab = tabs.find((tab) => tab.id === toolId);
      if (matchingTab) {
        setExpandedToolId(toolId);
      }
    }
  }, [toolId]);

  const handleToggle = (tabId: string) => {
    const newId = expandedToolId === tabId ? null : tabId;
    setExpandedToolId(newId);
    // Update URL without adding to history stack
    navigate(newId ? `/chem/${newId}` : "/chem", { replace: true });
  };

  return (
    <motion.div
      className="relative min-h-screen"
      variants={pageVariants}
      initial="hidden"
      animate="visible"
    >
      <GradientMesh page="chem" />

      {/* Page header */}
      <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8 py-8">
        <h1 className="text-3xl font-bold text-foreground">Chemical Analysis</h1>
        <p className="text-muted-foreground mt-2">
          Analyze molecules with descriptors, stereoisomers, similarity, and more.
        </p>
      </div>

      {/* Bento grid */}
      <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8 pb-12">
        <BentoGrid pageKey="chem" expandedId={expandedToolId}>
          {tabs.map((tab) => (
            <BentoCard
              key={tab.id}
              id={tab.id}
              pageKey="chem"
              name={tab.name}
              description={tab.description}
              icon={tab.icon}
              size={getBentoSize("chem", tab.id)}
              isExpanded={expandedToolId === tab.id}
              isOtherExpanded={expandedToolId !== null && expandedToolId !== tab.id}
              onToggle={() => handleToggle(tab.id)}
            >
              <tab.component />
            </BentoCard>
          ))}
        </BentoGrid>
      </div>
    </motion.div>
  );
};

export default ChemPage;
