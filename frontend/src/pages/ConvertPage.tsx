/**
 * ConvertPage -- bento grid layout for format conversion tools.
 *
 * Replaced tab navigation with BentoGrid expand-in-place UX.
 * Each tool is a card that expands to show its full UI when clicked.
 */
import { useState, useEffect } from "react";
import { motion } from "motion/react";
import { useParams, useNavigate } from "react-router-dom";

import FormatConversionView from "../components/convert/FormatConversionView";
import Mol2DView from "../components/convert/Mol2DView";
import Mol3DView from "../components/convert/Mol3DView";
import { Box, Copy, LayoutGrid } from "lucide-react";

import { BentoGrid, BentoCard } from "@/components/common/BentoGrid";
import { getBentoSize } from "@/config/bentoLayouts";
import { GradientMesh } from "@/components/common/GradientMesh";

// Tab data with icons and descriptions
const tabs = [
  {
    id: "format-conversion",
    name: "Format Conversion",
    component: FormatConversionView,
    icon: Copy,
    description: "Convert between different chemical file formats",
  },
  {
    id: "2d-coordinates",
    name: "2D Coordinates",
    component: Mol2DView,
    icon: LayoutGrid,
    description: "Generate 2D coordinates for molecular structures",
  },
  {
    id: "3d-coordinates",
    name: "3D Coordinates",
    component: Mol3DView,
    icon: Box,
    description: "Generate 3D coordinates for molecular structures",
  },
];

// --- Animation Variants ---
const pageVariants = {
  hidden: { opacity: 0 },
  visible: { opacity: 1, transition: { duration: 0.5 } },
};

const ConvertPage = () => {
  const { convertId } = useParams<{ convertId?: string }>();
  const navigate = useNavigate();
  const [expandedToolId, setExpandedToolId] = useState<string | null>(null);

  // Initialize expanded tool from URL parameter on mount
  useEffect(() => {
    if (convertId) {
      const matchingTab = tabs.find((tab) => tab.id === convertId);
      if (matchingTab) {
        setExpandedToolId(convertId);
      }
    }
  }, [convertId]);

  const handleToggle = (tabId: string) => {
    const newId = expandedToolId === tabId ? null : tabId;
    setExpandedToolId(newId);
    // Update URL without adding to history stack
    navigate(newId ? `/convert/${newId}` : "/convert", { replace: true });
  };

  return (
    <motion.div
      className="relative min-h-screen"
      variants={pageVariants}
      initial="hidden"
      animate="visible"
    >
      <GradientMesh page="convert" />

      {/* Page header */}
      <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8 py-8">
        <h1 className="text-3xl font-bold text-foreground">
          Format Conversion & Coordinate Generation
        </h1>
        <p className="text-muted-foreground mt-2">
          Convert between chemical file formats, generate 2D and 3D coordinates.
        </p>
      </div>

      {/* Bento grid */}
      <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8 pb-12">
        <BentoGrid pageKey="convert" expandedId={expandedToolId}>
          {tabs.map((tab) => (
            <BentoCard
              key={tab.id}
              id={tab.id}
              pageKey="convert"
              name={tab.name}
              description={tab.description}
              icon={tab.icon}
              size={getBentoSize("convert", tab.id)}
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

export default ConvertPage;
