/**
 * DepictPage -- bento grid layout for chemical structure depiction tools.
 *
 * Replaced tab navigation with BentoGrid expand-in-place UX.
 * Each tool is a card that expands to show its full UI when clicked.
 */
import { useState, useEffect } from "react";
import { motion } from "motion/react";
import { useParams, useNavigate } from "react-router-dom";

import Depict3DView from "../components/depict/Depict3DView";
import Depict2DMultiView from "../components/depict/Depict2DMultiView";
import StructureVisualizerView from "../components/depict/StructureVisualizerView";
import StructureDrawView from "../components/depict/StructureDrawView";
import { Box, LayoutGrid, Pencil, Search } from "lucide-react";

import { BentoGrid, BentoCard } from "@/components/common/BentoGrid";
import { getBentoSize } from "@/config/bentoLayouts";
import { GradientMesh } from "@/components/common/GradientMesh";

// Tab data
const tabs = [
  {
    id: "structureexplorer",
    name: "Structure Explorer",
    component: StructureVisualizerView,
    icon: Search,
    description: "Find structures by name or identifier and visualize them in 2D and 3D",
  },
  {
    id: "2ddepiction",
    name: "2D Depiction",
    component: Depict2DMultiView,
    icon: LayoutGrid,
    description: "Generate 2D depictions for multiple molecules at once.",
  },
  {
    id: "3ddepiction",
    name: "3D Depiction",
    component: Depict3DView,
    icon: Box,
    description: "Create interactive 3D visualizations",
  },
  {
    id: "structuredraw",
    name: "Draw a Structure",
    component: StructureDrawView,
    icon: Pencil,
    description: "Draw and edit chemical structures using a user-friendly interface",
  },
];

// --- Animation Variants ---
const pageVariants = {
  hidden: { opacity: 0 },
  visible: { opacity: 1, transition: { duration: 0.5 } },
};

const DepictPage = () => {
  const { depictId } = useParams<{ depictId?: string }>();
  const navigate = useNavigate();
  const [expandedToolId, setExpandedToolId] = useState<string | null>(null);

  // Initialize expanded tool from URL parameter on mount
  useEffect(() => {
    if (depictId) {
      const matchingTab = tabs.find((tab) => tab.id === depictId);
      if (matchingTab) {
        setExpandedToolId(depictId);
      }
    }
  }, [depictId]);

  const handleToggle = (tabId: string) => {
    const newId = expandedToolId === tabId ? null : tabId;
    setExpandedToolId(newId);
    // Update URL without adding to history stack
    navigate(newId ? `/depict/${newId}` : "/depict", { replace: true });
  };

  return (
    <motion.div
      className="relative min-h-screen"
      variants={pageVariants}
      initial="hidden"
      animate="visible"
    >
      <GradientMesh page="depict" />

      {/* Page header */}
      <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8 py-8">
        <h1 className="text-3xl font-bold text-foreground">Chemical Structure Depiction</h1>
        <p className="text-muted-foreground mt-2">
          Generate customizable 2D and interactive 3D visualizations of chemical structures.
        </p>
      </div>

      {/* Bento grid */}
      <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8 pb-12">
        <BentoGrid pageKey="depict" expandedId={expandedToolId}>
          {tabs.map((tab) => (
            <BentoCard
              key={tab.id}
              id={tab.id}
              pageKey="depict"
              name={tab.name}
              description={tab.description}
              icon={tab.icon}
              size={getBentoSize("depict", tab.id)}
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

export default DepictPage;
