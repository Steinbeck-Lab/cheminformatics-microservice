/**
 * ToolsPage -- bento grid layout for chemistry tools.
 *
 * Replaced tab navigation with BentoGrid expand-in-place UX.
 * Each tool is a card that expands to show its full UI when clicked.
 */
import { useState, useEffect } from "react";
import { motion } from "motion/react";
import { useParams, useNavigate } from "react-router-dom";

import SugarRemovalView from "../components/tools/SugarRemovalView";
import StructureGenView from "../components/tools/StructureGenView";
import InChIView from "../components/tools/InChIView";
import RInChIView from "../components/tools/RInChIView";
import { ArrowLeftRight, Box, FileText, Puzzle } from "lucide-react";

import { BentoGrid, BentoCard } from "@/components/common/BentoGrid";
import { getBentoSize } from "@/config/bentoLayouts";
import { GradientMesh } from "@/components/common/GradientMesh";

// Tab data with icons and descriptions
const tabs = [
  {
    id: "sugardetection",
    name: "Sugar Detection",
    component: SugarRemovalView,
    icon: Box,
    description: "Detect and remove sugar moieties from complex molecular structures.",
  },
  {
    id: "structuregeneration",
    name: "Structure Generation",
    component: StructureGenView,
    icon: Puzzle,
    description:
      "Generate chemical structures based on specified parameters and constraints for virtual screening.",
  },
  {
    id: "inchiconverter",
    name: "InChI Converter",
    component: InChIView,
    icon: FileText,
    description:
      "Draw, edit, and convert chemical structures to InChI notation with full support for various InChI versions and options.",
  },
  {
    id: "rinchiconverter",
    name: "RInChI Converter",
    component: RInChIView,
    icon: ArrowLeftRight,
    description:
      "Convert chemical reactions to RInChI notation, enabling standardized representation of chemical transformations.",
  },
];

// --- Animation Variants ---
const pageVariants = {
  hidden: { opacity: 0 },
  visible: { opacity: 1, transition: { duration: 0.5 } },
};

const ToolsPage = () => {
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
    navigate(newId ? `/tools/${newId}` : "/tools", { replace: true });
  };

  return (
    <motion.div
      className="relative min-h-screen"
      variants={pageVariants}
      initial="hidden"
      animate="visible"
    >
      <GradientMesh page="tools" />

      {/* Page header */}
      <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8 py-8">
        <h1 className="text-3xl font-bold text-foreground">Chemistry Tools</h1>
        <p className="text-muted-foreground mt-2">
          Specialized tools for structure generation, sugar moiety removal, and format conversion.
        </p>
      </div>

      {/* Bento grid */}
      <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8 pb-12">
        <BentoGrid pageKey="tools" expandedId={expandedToolId}>
          {tabs.map((tab) => (
            <BentoCard
              key={tab.id}
              id={tab.id}
              pageKey="tools"
              name={tab.name}
              description={tab.description}
              icon={tab.icon}
              size={getBentoSize("tools", tab.id)}
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

export default ToolsPage;
