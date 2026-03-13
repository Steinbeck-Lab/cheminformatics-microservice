/**
 * ToolsPage -- tabbed interface for chemistry tools.
 *
 * Reverted from bento grid back to tab navigation per user feedback.
 * Uses GradientMesh for background and glass-bold for glass surfaces.
 */
import React, { useState, useEffect } from "react";
import { useParams, useNavigate } from "react-router-dom";
import { motion, AnimatePresence, LayoutGroup } from "motion/react";
import SugarRemovalView from "../components/tools/SugarRemovalView";
import StructureGenView from "../components/tools/StructureGenView";
import InChIView from "../components/tools/InChIView";
import RInChIView from "../components/tools/RInChIView";
import { ArrowLeftRight, Box, ChevronDown, FileText, Puzzle } from "lucide-react";
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
  visible: { opacity: 1, transition: { duration: 0.6, ease: "easeInOut" } },
};
const headerContainerVariants = {
  hidden: { opacity: 0 },
  visible: {
    opacity: 1,
    transition: { staggerChildren: 0.1, delayChildren: 0.15 },
  },
};
const headerItemVariants = {
  hidden: { opacity: 0, y: -15 },
  visible: {
    opacity: 1,
    y: 0,
    transition: { duration: 0.6, ease: [0.2, 0.65, 0.3, 0.9] },
  },
};
const tabContainerVariant = {
  hidden: { opacity: 0, y: 25, scale: 0.97 },
  visible: {
    opacity: 1,
    y: 0,
    scale: 1,
    transition: { duration: 0.6, delay: 0.25, ease: [0.2, 0.65, 0.3, 0.9] },
  },
};
const contentVariants = {
  hidden: { opacity: 0, x: -20 },
  visible: {
    opacity: 1,
    x: 0,
    transition: { type: "spring", stiffness: 500, damping: 30 },
  },
  exit: { opacity: 0, x: 20, transition: { duration: 0.3, ease: "easeIn" } },
};
const mobileMenuVariants = {
  closed: { opacity: 0, y: -10, height: 0 },
  open: {
    opacity: 1,
    y: 0,
    height: "auto",
    transition: { duration: 0.3, ease: "easeOut" },
  },
};

const ToolsPage = () => {
  const { toolId } = useParams<{ toolId?: string }>();
  const navigate = useNavigate();

  const [activeTabId, setActiveTabId] = useState(toolId || tabs[0].id);
  const [isMobileMenuOpen, setIsMobileMenuOpen] = useState(false);
  const [isMobile, setIsMobile] = useState(false);

  // Redirect to first tab if no toolId, or update active tab when URL changes
  useEffect(() => {
    if (!toolId) {
      navigate(`/tools/${tabs[0].id}`, { replace: true });
    } else {
      const tab = tabs.find((t) => t.id === toolId);
      if (tab) {
        setActiveTabId(toolId);
      } else {
        navigate(`/tools/${tabs[0].id}`, { replace: true });
      }
    }
  }, [toolId, navigate]);

  // Check if the window is mobile size
  useEffect(() => {
    const checkIfMobile = () => {
      setIsMobile(window.innerWidth < 768);
    };
    checkIfMobile();
    window.addEventListener("resize", checkIfMobile);
    return () => window.removeEventListener("resize", checkIfMobile);
  }, []);

  const activeTab = tabs.find((tab) => tab.id === activeTabId);
  const ActiveComponent = activeTab ? activeTab.component : null;

  const handleTabSelection = (tabId: string) => {
    setActiveTabId(tabId);
    setIsMobileMenuOpen(false);
    navigate(`/tools/${tabId}`);
  };

  return (
    <motion.div
      className="relative min-h-screen"
      variants={pageVariants}
      initial="hidden"
      animate="visible"
    >
      <GradientMesh page="tools" />

      {/* Content Area */}
      <div className="relative w-full mx-auto px-4 sm:px-6 lg:px-8 py-8 md:py-12 z-10">
        <div className="w-full lg:w-3/4 mx-auto">
          {/* Page Header */}
          <motion.div
            className="mb-6 md:mb-8 max-w-5xl mx-auto text-center"
            variants={headerContainerVariants}
            initial="hidden"
            animate="visible"
          >
            <motion.h1
              variants={headerItemVariants}
              className="text-2xl md:text-4xl lg:text-5xl font-bold text-foreground mb-2 md:mb-3"
            >
              Chemistry Tools
            </motion.h1>
            <motion.p
              variants={headerItemVariants}
              className="text-muted-foreground text-sm md:text-lg max-w-3xl mx-auto"
            >
              Specialized tools for structure generation, sugar moiety removal, and format
              conversion.
            </motion.p>
          </motion.div>

          {/* Tab Container */}
          <motion.div
            className="glass-bold rounded-xl shadow-2xl overflow-hidden"
            variants={tabContainerVariant}
            initial="hidden"
            animate="visible"
          >
            {/* Tab Navigation */}
            <div className="relative border-b border-white/10 bg-white/5 dark:bg-slate-800/20">
              {/* Mobile Tab Navigation (Dropdown) */}
              <div className="block md:hidden">
                <button
                  onClick={() => setIsMobileMenuOpen(!isMobileMenuOpen)}
                  className="flex items-center justify-between w-full py-3 px-4 text-left focus:outline-hidden"
                  aria-expanded={isMobileMenuOpen}
                >
                  <div className="flex items-center">
                    {activeTab && (
                      <>
                        <activeTab.icon className="h-5 w-5 mr-2 text-sky-600 dark:text-sky-400" />
                        <span className="font-medium text-sky-700 dark:text-white">
                          {activeTab.name}
                        </span>
                      </>
                    )}
                  </div>
                  <motion.span
                    animate={{ rotate: isMobileMenuOpen ? 180 : 0 }}
                    transition={{ duration: 0.3 }}
                  >
                    <ChevronDown className="h-5 w-5 text-slate-500 dark:text-slate-400" />
                  </motion.span>
                </button>

                <AnimatePresence>
                  {isMobileMenuOpen && (
                    <motion.div
                      variants={mobileMenuVariants}
                      initial="closed"
                      animate="open"
                      exit="closed"
                      className="border-t border-white/10 overflow-hidden"
                    >
                      {tabs.map((tab) => (
                        <button
                          key={tab.id}
                          onClick={() => handleTabSelection(tab.id)}
                          className={`w-full flex items-center py-3 px-4 ${
                            activeTabId === tab.id
                              ? "bg-sky-100/20 dark:bg-slate-700/50 text-sky-700 dark:text-white"
                              : "text-slate-600 dark:text-slate-400 hover:bg-white/10"
                          }`}
                        >
                          <tab.icon
                            className={`h-5 w-5 mr-3 ${
                              activeTabId === tab.id
                                ? "text-sky-600 dark:text-sky-400"
                                : "text-slate-500 dark:text-slate-400"
                            }`}
                          />
                          <span>{tab.name}</span>
                        </button>
                      ))}
                    </motion.div>
                  )}
                </AnimatePresence>
              </div>

              {/* Desktop Tab Navigation (Horizontal Tabs) */}
              <div className="hidden md:block">
                <div className="flex justify-center overflow-x-auto py-3 px-4 space-x-3">
                  <LayoutGroup id="tools-tabs">
                    {tabs.map((tab) => {
                      const isActive = activeTabId === tab.id;
                      return (
                        <motion.button
                          key={tab.id}
                          onClick={() => handleTabSelection(tab.id)}
                          className={`relative flex items-center px-5 py-2.5 text-sm font-medium rounded-lg transition-colors duration-200 whitespace-nowrap focus:outline-hidden focus-visible:ring-2 focus-visible:ring-offset-2 focus-visible:ring-offset-background focus-visible:ring-primary ${
                            isActive
                              ? "text-sky-700 dark:text-white"
                              : "text-slate-600 dark:text-slate-400 hover:text-slate-900 dark:hover:text-white"
                          }`}
                          aria-selected={isActive}
                          role="tab"
                          whileHover={isMobile ? {} : { scale: 1.04 }}
                          whileTap={{ scale: 0.96 }}
                          transition={{
                            type: "spring",
                            stiffness: 400,
                            damping: 20,
                            mass: 0.5,
                          }}
                        >
                          {isActive && (
                            <motion.div
                              className="absolute inset-0 bg-white/80 dark:bg-slate-700/90 rounded-lg shadow-xs"
                              layoutId="activeTabPill"
                              transition={{
                                type: "spring",
                                stiffness: 380,
                                damping: 30,
                                mass: 0.6,
                              }}
                              style={{ zIndex: 0 }}
                            />
                          )}
                          <span className="relative z-10 flex items-center">
                            <tab.icon
                              className={`h-5 w-5 mr-2 shrink-0 transition-colors duration-200 ${
                                isActive
                                  ? "text-sky-600 dark:text-sky-400"
                                  : "text-slate-500 dark:text-slate-400"
                              }`}
                            />
                            <span>{tab.name}</span>
                          </span>
                        </motion.button>
                      );
                    })}
                  </LayoutGroup>
                </div>
              </div>
            </div>

            {/* Tab Content Header */}
            {activeTab && (
              <motion.div
                key={`${activeTabId}-header`}
                initial={{ opacity: 0, y: -10 }}
                animate={{ opacity: 1, y: 0 }}
                transition={{ duration: 0.4, delay: 0.15 }}
                className="mx-3 sm:mx-4 mt-3 sm:mt-4 mb-2"
              >
                <div className="bg-linear-to-r from-green-600 to-emerald-700 dark:from-green-700 dark:to-emerald-900 rounded-xl shadow-xl overflow-hidden relative">
                  <div className="absolute inset-0 overflow-hidden opacity-10">
                    <div className="absolute left-0 top-0 w-32 h-32 rounded-full bg-white dark:bg-slate-300 transform -translate-x-1/2 -translate-y-1/2"></div>
                    <div className="absolute right-16 bottom-0 w-48 h-48 rounded-full bg-white dark:bg-slate-300 transform translate-x-1/2 translate-y-1/2"></div>
                    <div className="absolute left-1/3 top-1/2 w-24 h-24 rounded-full bg-white dark:bg-slate-300 transform -translate-y-1/2"></div>
                  </div>

                  <div className="relative p-4 md:p-6 flex flex-col md:flex-row items-center justify-between gap-3 md:gap-4">
                    <div className="space-y-1 text-center md:text-left max-w-3xl">
                      <h1 className="text-xl md:text-3xl font-bold leading-tight text-white dark:text-slate-50">
                        {activeTab.name}
                      </h1>
                      <p className="text-green-50 dark:text-emerald-100 text-xs md:text-base opacity-90">
                        {activeTab.description}
                      </p>
                    </div>

                    <div className="hidden md:flex items-center justify-center bg-white/15 dark:bg-slate-200/10 backdrop-blur-sm p-4 rounded-full w-16 h-16 border border-white/20 dark:border-slate-300/20 shadow-lg">
                      <activeTab.icon className="w-8 h-8 text-white dark:text-slate-50" />
                    </div>
                  </div>
                </div>
              </motion.div>
            )}

            {/* Tab Content Area */}
            <div className="relative overflow-hidden min-h-[50vh]">
              <AnimatePresence mode="wait">
                <motion.div
                  key={activeTabId}
                  initial="hidden"
                  animate="visible"
                  exit="exit"
                  variants={contentVariants}
                  className="p-3 sm:p-5 md:p-6"
                >
                  {ActiveComponent && <ActiveComponent isActive={true} isMobile={isMobile} />}
                </motion.div>
              </AnimatePresence>
            </div>
          </motion.div>
        </div>
      </div>
    </motion.div>
  );
};
export default ToolsPage;
