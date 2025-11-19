// Description: This file contains the ChemPage component, which provides a user interface for various chemical analysis tools. It includes a sidebar for navigation and displays different views based on the selected tool.
import React, { useState, useEffect, useMemo } from "react";
import { motion, AnimatePresence, LayoutGroup } from "framer-motion";

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

// Import icons
import {
  HiOutlineBeaker,
  HiOutlineDatabase,
  HiOutlineChartSquareBar,
  HiOutlineFingerPrint,
  HiOutlineTag,
  HiOutlineRefresh,
  HiOutlineFilter,
  HiOutlineCollection,
  HiOutlineCode,
  HiOutlineCheckCircle,
  HiOutlineViewList,
  HiOutlineTemplate,
  HiOutlineMenuAlt1,
  HiOutlineX,
  HiOutlineDocumentDuplicate,
  HiOutlineGlobeAlt,
  HiOutlineLightningBolt,
} from "react-icons/hi";

// Define tab data with icons, categories, and components
const tabs = [
  // Structure
  {
    id: "stereoisomers",
    name: "Stereoisomers",
    description: "Generate all possible stereoisomers",
    icon: HiOutlineDocumentDuplicate,
    category: "structure",
    component: StereoisomersView,
  },
  {
    id: "hosecodes",
    name: "HOSE Codes",
    description: "Generate HOSE codes",
    icon: HiOutlineCode,
    category: "structure",
    component: HOSECodeView,
  },
  {
    id: "standardize",
    name: "Standardize",
    description: "Standardize structures",
    icon: HiOutlineRefresh,
    category: "structure",
    component: StandardizeView,
  },
  {
    id: "functional-groups",
    name: "Functional Groups",
    description: "Identify functional groups",
    icon: HiOutlineCollection,
    category: "structure",
    component: ErtlFunctionalGroupView,
  },
  {
    id: "tautomer",
    name: "Tautomers",
    description: "Generate standardized tautomers",
    icon: HiOutlineRefresh,
    category: "structure",
    component: StandardizedTautomerView,
  },
  {
    id: "fixradicals",
    name: "Fix Radicals",
    description: "Fix radical electrons in molecules",
    icon: HiOutlineLightningBolt,
    category: "structure",
    component: FixRadicalsView,
  },
  // Analysis
  {
    id: "descriptors",
    name: "Descriptors",
    description: "Calculate molecular descriptors",
    icon: HiOutlineChartSquareBar,
    category: "analysis",
    component: DescriptorsView,
  },
  {
    id: "nplikeness",
    name: "NP-likeness",
    description: "Calculate NP-likeness score",
    icon: HiOutlineBeaker,
    category: "analysis",
    component: NPlikenessView,
  },
  // Comparison
  {
    id: "similarity",
    name: "Similarity",
    description: "Compare structures (Tanimoto)",
    icon: HiOutlineFingerPrint,
    category: "comparison",
    component: TanimotoView,
  },
  // Search
  {
    id: "structure-finder",
    name: "Structure Finder",
    description: "Find chemical structures by name or identifier",
    icon: HiOutlineGlobeAlt,
    category: "search",
    component: PubChemLookupView,
  },
  // Validation
  {
    id: "structureerror",
    name: "Check Structure",
    description: "Validate chemical structures",
    icon: HiOutlineCheckCircle,
    category: "validation",
    component: StructureErrorView,
  },
  {
    id: "filters",
    name: "All Filters",
    description: "Apply multiple chemical filters",
    icon: HiOutlineFilter,
    category: "validation",
    component: AllFiltersView,
  },
  // Advanced
  {
    id: "coconut",
    name: "COCONUT Preprocessing",
    description: "Prepare data for COCONUT DB",
    icon: HiOutlineDatabase,
    category: "advanced",
    component: CoconutPreProcessingView,
  },
  {
    id: "classyfire",
    name: "ClassyFire",
    description: "Chemical classification",
    icon: HiOutlineTag,
    category: "advanced",
    component: ClassyfireView,
  },
];

// Define categories with colors
const categories = {
  structure: { name: "Structure Tools", color: "blue" },
  analysis: { name: "Property Analysis", color: "purple" },
  comparison: { name: "Comparison", color: "indigo" },
  validation: { name: "Validation & Filtering", color: "green" },
  advanced: { name: "Advanced / Specific", color: "orange" },
  search: { name: "Search & Retrieval", color: "yellow" },
};
const categoryOrder = ["structure", "search", "analysis", "comparison", "validation", "advanced"];

// --- Animation Variants ---
const pageVariants = {
  hidden: { opacity: 0 },
  visible: { opacity: 1, transition: { duration: 0.5 } },
};
const contentVariants = {
  hidden: { opacity: 0, y: 15, scale: 0.98 },
  visible: {
    opacity: 1,
    y: 0,
    scale: 1,
    transition: { duration: 0.5, ease: [0.25, 1, 0.5, 1] },
  },
};
// FIX: Removed unused categoryContentVariants
// const categoryContentVariants = { ... };
const sidebarItemVariant = {
  hidden: { opacity: 0, x: -15 },
  visible: { opacity: 1, x: 0, transition: { duration: 0.3, ease: "easeOut" } },
};
const sidebarStaggerContainer = {
  hidden: { opacity: 0 },
  visible: {
    opacity: 1,
    transition: { staggerChildren: 0.05, delayChildren: 0.1 },
  },
};

const ChemPage = () => {
  const [activeTabId, setActiveTabId] = useState(tabs[0].id);
  const [viewMode, setViewMode] = useState("categories");
  const [isMobile, setIsMobile] = useState(false);
  const [sidebarOpen, setSidebarOpen] = useState(true);
  // FIX: Removed unused expandedCategories state
  // const [expandedCategories, setExpandedCategories] = useState(...);

  const activeTab = tabs.find((tab) => tab.id === activeTabId);
  const ActiveComponent = activeTab ? activeTab.component : () => <div>Select a tool</div>;

  const tabsByCategory = useMemo(() => {
    return tabs.reduce((acc, tab) => {
      acc[tab.category] = acc[tab.category] || [];
      acc[tab.category].push(tab);
      return acc;
    }, {});
  }, []);

  useEffect(() => {
    const checkMobile = () => {
      const mobile = window.innerWidth < 768;
      setIsMobile(mobile);
      setSidebarOpen(!mobile);
    };
    checkMobile();
    window.addEventListener("resize", checkMobile);
    return () => window.removeEventListener("resize", checkMobile);
  }, []);

  const handleTabClick = (tabId) => {
    setActiveTabId(tabId);
    if (isMobile) setSidebarOpen(false);
  };

  const handleSidebarToggle = () => {
    setSidebarOpen((prev) => !prev);
  };

  // FIX: Removed unused toggleCategory function
  // const toggleCategory = (category) => { ... };

  const getCategoryColorClass = (category, type = "text") => {
    const colorName = categories[category]?.color || "blue";
    if (type === "text") return `text-${colorName}-600 dark:text-${colorName}-400`;
    if (type === "border") return `border-${colorName}-500 dark:border-${colorName}-400`;
    return "";
  };

  return (
    <motion.div
      className="relative min-h-screen w-full bg-slate-100 dark:bg-gray-950 text-slate-900 dark:text-slate-100 font-sans overflow-hidden isolate"
      variants={pageVariants}
      initial="hidden"
      animate="visible"
    >
      {/* Background Effects */}
      <div className="absolute inset-0 -z-20 overflow-hidden dark:opacity-100 opacity-0 transition-opacity duration-500">
        <div className="absolute inset-0 bg-gradient-to-br from-gray-950 via-slate-900 to-indigo-950"></div>
        <div className="animated-mesh-gradient"></div>
      </div>
      <div className="absolute inset-0 -z-20 overflow-hidden dark:opacity-0 opacity-100 transition-opacity duration-500">
        <div className="absolute inset-0 bg-gradient-to-br from-sky-50 via-white to-indigo-50"></div>
        <div
          className="absolute inset-0 opacity-[0.025]"
          style={{
            backgroundImage: `url("data:image/svg+xml,%3Csvg width='70' height='70' viewBox='0 0 70 70' xmlns='http://www.w3.org/2000/svg'%3E%3Cg fill='none' fill-rule='evenodd'%3E%3Cg fill='%23BBB' fill-opacity='.2'%3E%3Cpath d='M35 0v70M0 35h70' stroke-width='1'/%3E%3C/g%3E%3C/g%3E%3C/svg%3E")`,
          }}
        />
      </div>
      <div
        className="absolute inset-0 -z-10 opacity-[0.02] dark:opacity-[0.03]"
        style={{ backgroundImage: "url(/noise.svg)" }}
      ></div>
      {/* Mobile toggle button */}
      <motion.button
        className={`md:hidden fixed top-20 left-4 z-[60] p-2 rounded-full shadow-lg transition-colors duration-200 ${sidebarOpen ? "bg-slate-300/80 dark:bg-slate-700/80 text-slate-800 dark:text-slate-100 hover:bg-slate-400/80 dark:hover:bg-slate-600/80" : "bg-sky-600/80 dark:bg-blue-600/80 text-white hover:bg-sky-700/90 dark:hover:bg-blue-700/90"}`}
        onClick={handleSidebarToggle}
        aria-label={sidebarOpen ? "Close sidebar" : "Open sidebar"}
        whileTap={{ scale: 0.95 }}
      >
        <AnimatePresence mode="wait" initial={false}>
          {sidebarOpen ? (
            <motion.div
              key="close-sb"
              initial={{ rotate: -90, opacity: 0 }}
              animate={{ rotate: 0, opacity: 1 }}
              exit={{ rotate: 90, opacity: 0 }}
              transition={{ duration: 0.2 }}
            >
              {" "}
              <HiOutlineX className="h-6 w-6" />{" "}
            </motion.div>
          ) : (
            <motion.div
              key="open-sb"
              initial={{ rotate: 90, opacity: 0 }}
              animate={{ rotate: 0, opacity: 1 }}
              exit={{ rotate: -90, opacity: 0 }}
              transition={{ duration: 0.2 }}
            >
              {" "}
              <HiOutlineMenuAlt1 className="h-6 w-6" />{" "}
            </motion.div>
          )}
        </AnimatePresence>
      </motion.button>
      {/* Main Layout Container */}
      <div className="relative w-full mx-auto px-4 sm:px-6 lg:px-8 py-12 md:py-16 z-10">
        <div className="w-full lg:w-3/4 mx-auto">
          <motion.div
            className="mb-8 md:mb-10 text-center mx-auto"
            initial={{ opacity: 0 }}
            animate={{ opacity: 1 }}
            transition={{ delay: 0.1 }}
          >
            <h1 className="text-3xl md:text-4xl lg:text-5xl font-bold text-[var(--text-primary)] mb-3">
              Chemical Analysis Tools
            </h1>
            <p className="text-[var(--text-secondary)] text-base md:text-lg max-w-3xl mx-auto">
              Comprehensive toolkit for molecular analysis, classification, and validation.
            </p>
          </motion.div>
          {/* Flex container for Sidebar + Main Content */}
          <div className="flex flex-col md:flex-row md:gap-x-4 lg:gap-x-6">
            {/* Sidebar */}
            <div
              className={`glass md:w-64 lg:w-72 flex-shrink-0 rounded-xl shadow-lg border border-slate-200 dark:border-slate-700/50 overflow-hidden
                        fixed inset-y-0 left-0 z-40 w-64 mt-16 md:mt-0 md:sticky md:top-20 md:h-[calc(100vh-10rem)] transform transition-transform duration-300 ease-in-out md:transform-none ${
                          sidebarOpen ? "translate-x-0" : "-translate-x-full"
                        }`}
            >
              <div className="h-full overflow-y-auto">
                {/* View Mode Toggle */}
                <div className="sticky top-0 z-10 p-2 bg-white/50 dark:bg-slate-800/50 backdrop-blur-sm border-b border-slate-200 dark:border-slate-700/50">
                  <LayoutGroup id="viewModeToggle-chem">
                    <div className="flex p-1 bg-slate-200 dark:bg-slate-900 rounded-lg">
                      <motion.button
                        className={`relative flex-1 py-1.5 px-3 rounded-md font-medium text-xs sm:text-sm transition-colors ${viewMode === "categories" ? "text-sky-800 dark:text-white" : "text-slate-600 dark:text-slate-400 hover:text-slate-900 dark:hover:text-white"}`}
                        onClick={() => setViewMode("categories")}
                      >
                        {viewMode === "categories" && (
                          <motion.div
                            layoutId="viewModePill-chem"
                            className="absolute inset-0 bg-white dark:bg-slate-700 rounded-md shadow-sm"
                            transition={{
                              type: "spring",
                              stiffness: 400,
                              damping: 35,
                            }}
                          />
                        )}
                        <span className="relative z-10 flex items-center justify-center gap-1">
                          <HiOutlineTemplate className="h-4 w-4" /> Categories
                        </span>
                      </motion.button>
                      <motion.button
                        className={`relative flex-1 py-1.5 px-3 rounded-md font-medium text-xs sm:text-sm transition-colors ${viewMode === "all" ? "text-sky-800 dark:text-white" : "text-slate-600 dark:text-slate-400 hover:text-slate-900 dark:hover:text-white"}`}
                        onClick={() => setViewMode("all")}
                      >
                        {viewMode === "all" && (
                          <motion.div
                            layoutId="viewModePill-chem"
                            className="absolute inset-0 bg-white dark:bg-slate-700 rounded-md shadow-sm"
                            transition={{
                              type: "spring",
                              stiffness: 400,
                              damping: 35,
                            }}
                          />
                        )}
                        <span className="relative z-10 flex items-center justify-center gap-1">
                          <HiOutlineViewList className="h-4 w-4" /> All Tools
                        </span>
                      </motion.button>
                    </div>
                  </LayoutGroup>
                </div>

                {/* Navigation List */}
                <motion.div
                  className="p-3"
                  variants={sidebarStaggerContainer}
                  initial="hidden"
                  animate="visible"
                >
                  <LayoutGroup id="toolNav-chem">
                    {viewMode === "categories" ? (
                      <div className="space-y-3">
                        {categoryOrder.map(
                          (category) =>
                            tabsByCategory[category] && (
                              <motion.div
                                key={category}
                                variants={sidebarItemVariant}
                                className="border-b border-slate-300/70 dark:border-slate-700/50 pb-3 last:border-b-0"
                              >
                                {/* Removed category toggle button */}
                                <h3
                                  className={`text-xs font-semibold px-2 py-1 ${getCategoryColorClass(category, "text")} uppercase tracking-wider mb-1`}
                                >
                                  {" "}
                                  {categories[category]?.name || category}{" "}
                                </h3>
                                {/* Always show items, removed AnimatePresence & variants */}
                                <div className="mt-2 space-y-1">
                                  {tabsByCategory[category].map((tab) => {
                                    const isActive = activeTabId === tab.id;
                                    return (
                                      <button
                                        key={tab.id}
                                        onClick={() => handleTabClick(tab.id)}
                                        className={`relative w-full flex items-center px-3 py-2 text-sm rounded-md transition-colors group ${isActive ? "text-sky-700 dark:text-white" : "text-slate-600 dark:text-slate-300 hover:text-slate-900 dark:hover:text-white hover:bg-slate-200/40 dark:hover:bg-slate-700/40"}`}
                                      >
                                        {" "}
                                        {isActive && (
                                          <motion.div
                                            className="absolute inset-0 bg-sky-100 dark:bg-slate-700 rounded-md -z-10"
                                            layoutId="activeToolIndicator-chem"
                                            transition={{
                                              type: "spring",
                                              stiffness: 350,
                                              damping: 30,
                                            }}
                                          />
                                        )}{" "}
                                        <tab.icon
                                          className={`h-5 w-5 mr-2 flex-shrink-0 ${isActive ? getCategoryColorClass(category, "text") : "text-slate-500 dark:text-slate-400 group-hover:text-slate-600 dark:group-hover:text-slate-300"}`}
                                        />{" "}
                                        <span className="truncate">{tab.name}</span>{" "}
                                      </button>
                                    );
                                  })}
                                </div>
                              </motion.div>
                            )
                        )}
                      </div>
                    ) : (
                      // All Tools View
                      <div className="space-y-1">
                        {tabs.map((tab) => {
                          const isActive = activeTabId === tab.id;
                          return (
                            <motion.button
                              key={tab.id}
                              variants={sidebarItemVariant}
                              onClick={() => handleTabClick(tab.id)}
                              className={`relative w-full flex items-center px-3 py-2 text-sm rounded-md transition-colors group ${isActive ? "text-sky-700 dark:text-white" : "text-slate-600 dark:text-slate-300 hover:text-slate-900 dark:hover:text-white hover:bg-slate-200/40 dark:hover:bg-slate-700/40"}`}
                            >
                              {" "}
                              {isActive && (
                                <motion.div
                                  className="absolute inset-0 bg-sky-100 dark:bg-slate-700 rounded-md -z-10"
                                  layoutId="activeToolIndicator-chem"
                                  transition={{
                                    type: "spring",
                                    stiffness: 350,
                                    damping: 30,
                                  }}
                                />
                              )}{" "}
                              <tab.icon
                                className={`h-5 w-5 mr-2 flex-shrink-0 ${isActive ? getCategoryColorClass(tab.category, "text") : "text-slate-500 dark:text-slate-400 group-hover:text-slate-600 dark:group-hover:text-slate-300"}`}
                              />{" "}
                              <span className="truncate">{tab.name}</span>{" "}
                            </motion.button>
                          );
                        })}
                      </div>
                    )}
                  </LayoutGroup>
                </motion.div>
              </div>
            </div>{" "}
            {/* End Sidebar Div */}
            {/* Main content Area */}
            <div className="flex-1 min-w-0 mt-6 md:mt-0">
              <AnimatePresence mode="wait">
                <motion.div
                  key={activeTabId}
                  className="glass rounded-xl shadow-lg border border-slate-200 dark:border-slate-700/50 overflow-hidden min-h-[calc(100vh-12rem)]"
                  initial="hidden"
                  animate="visible"
                  exit="hidden"
                  variants={contentVariants}
                >
                  {/* Tab Content Header */}
                  {activeTab && (
                    <div className="p-5 sm:p-6 border-b border-slate-200/80 dark:border-slate-700/50 bg-white/30 dark:bg-slate-800/20">
                      <div className="flex items-start sm:items-center">
                        <activeTab.icon
                          className={`h-6 w-6 sm:h-7 sm:w-7 mr-3 sm:mr-4 flex-shrink-0 mt-0.5 sm:mt-0 ${getCategoryColorClass(activeTab.category, "text")}`}
                        />
                        <div>
                          <h1 className="text-xl sm:text-2xl font-bold text-[var(--text-primary)]">
                            {activeTab.name}
                          </h1>
                          <p className="text-sm text-[var(--text-secondary)] mt-1">
                            {activeTab.description}
                          </p>
                        </div>
                      </div>
                    </div>
                  )}

                  {/* Tab Content */}
                  <div className="p-5 sm:p-6">
                    <ActiveComponent />
                  </div>
                </motion.div>
              </AnimatePresence>
            </div>{" "}
            {/* End Main Content Flex Item */}
          </div>{" "}
          {/* End Main Layout Flex Container */}
        </div>{" "}
        {/* End Width Constraint Wrapper */}
      </div>{" "}
      {/* End Padding Container */}
      {/* Global Styles */}
      <style jsx="true" global="true">{`
        /* Ensure required CSS variables and utilities like .glass are defined */
        /* --- Animated Mesh Gradient (Dark Mode Only) --- */
        @keyframes mesh-gradient-move {
          /* ... */
        }
        .animated-mesh-gradient {
          /* ... */
        }
        .dark .animated-mesh-gradient {
          /* ... */
        }
        /* Noise overlay needs noise.svg in public folder */
      `}</style>
    </motion.div> // End Page Container
  );
};

export default ChemPage;
