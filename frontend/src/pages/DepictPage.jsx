// Description: DepictPage component for generating 2D and 3D depictions of chemical structures
import React, { useState } from "react";
import {
  motion,
  AnimatePresence,
  LayoutGroup,
  useScroll,
  useTransform,
} from "framer-motion";
import Depict3DView from "../components/depict/Depict3DView";
import Depict2DMultiView from "../components/depict/Depict2DMultiView";
import StructureVisualizerView from "../components/depict/StructureVisualizerView";
//import StructureSketcherView from '../components/depict/StructureDrawView';
// Import icons
import {
  HiOutlineViewGrid,
  HiOutlineCube,
  HiOutlineSearch,
} from "react-icons/hi";

// Tab data
const tabs = [
  {
    id: "batch-depiction",
    name: "2D Depiction",
    component: Depict2DMultiView,
    icon: HiOutlineViewGrid,
    description: `Generate 2D depictions for multiple molecules at once. 
                 Enter one SMILES string per line. 
                 Optional: add a title after each SMILES separated by a space or tab.`,
  },
  {
    id: "3d-depiction",
    name: "3D Depiction",
    component: Depict3DView,
    icon: HiOutlineCube,

    description: "Create interactive 3D visualizations",
  },
  {
    id: "structure-explorer",
    name: "Structure Explorer",
    component: StructureVisualizerView,
    icon: HiOutlineSearch,
    description:
      "Find structures by name or identifier and visualize them in 2D and 3D",
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
    transition: { duration: 0.45, ease: [0.25, 1, 0.5, 1] },
  },
  exit: { opacity: 0, x: 20, transition: { duration: 0.3, ease: "easeIn" } },
};

const DepictPage = () => {
  const [activeTabId, setActiveTabId] = useState(tabs[0].id);

  const activeTab = tabs.find((tab) => tab.id === activeTabId);
  const ActiveComponent = activeTab ? activeTab.component : null;

  const { scrollYProgress } = useScroll();
  const meshY = useTransform(scrollYProgress, [0, 1], ["0%", "15%"]);
  const noiseY = useTransform(scrollYProgress, [0, 1], ["0%", "30%"]);

  return (
    <motion.div
      className="relative min-h-screen w-full bg-slate-100 dark:bg-gray-950 text-slate-900 dark:text-slate-100 font-sans overflow-x-hidden isolate"
      variants={pageVariants}
      initial="hidden"
      animate="visible"
    >
      {/* --- Adaptive Background Effects --- */}
      <div className="absolute inset-0 -z-20 overflow-hidden dark:opacity-100 opacity-0 transition-opacity duration-500">
        <div className="absolute inset-0 bg-gradient-to-br from-gray-950 via-slate-900 to-indigo-950"></div>
        <motion.div className="animated-mesh-gradient" style={{ y: meshY }} />
      </div>
      <div className="absolute inset-0 -z-20 overflow-hidden dark:opacity-0 opacity-100 transition-opacity duration-500">
        <div className="absolute inset-0 bg-gradient-to-br from-sky-50 via-white to-indigo-100"></div>
        <motion.div
          className="absolute inset-0 opacity-[0.03]"
          style={{
            y: noiseY,
            backgroundImage: `url("data:image/svg+xml,%3Csvg width='70' height='70' viewBox='0 0 70 70' xmlns='http://www.w3.org/2000/svg'%3E%3Cg fill='none' fill-rule='evenodd'%3E%3Cg fill='%23AAB' fill-opacity='.25'%3E%3Cpath d='M35 0v70M0 35h70' stroke-width='1.5'/%3E%3C/g%3E%3C/g%3E%3C/svg%3E")`,
          }}
        />
      </div>
      <motion.div
        className="absolute inset-0 -z-10 opacity-[0.02] dark:opacity-[0.03]"
        style={{ y: noiseY, backgroundImage: "url(/noise.svg)" }}
      ></motion.div>
      {/* Content Area - Outer padding container */}
      <div className="relative w-full mx-auto px-4 sm:px-6 lg:px-8 py-12 md:py-16 z-10">
        {/* FIX: Wrapper div for width constraint (matches ChemPage: lg:w-3/4) */}
        <div className="w-full lg:w-3/4 mx-auto">
          {/* Page Header - Animated */}
          <motion.div
            className="mb-8 md:mb-10 max-w-5xl mx-auto text-center" // Centered within the 3/4 width
            variants={headerContainerVariants}
            initial="hidden"
            animate="visible"
          >
            <motion.h1
              variants={headerItemVariants}
              className="text-3xl md:text-4xl lg:text-5xl font-bold text-[var(--text-primary)] mb-3"
            >
              Chemical Structure Depiction
            </motion.h1>
            <motion.p
              variants={headerItemVariants}
              className="text-[var(--text-secondary)] text-base md:text-lg max-w-3xl mx-auto"
            >
              Generate customizable 2D and interactive 3D visualizations of
              chemical structures.
            </motion.p>
          </motion.div>
          {/* Tab Container - Animated */}
          <motion.div
            className="glass bg-white/75 dark:bg-slate-800/80 backdrop-blur-xl dark:backdrop-blur-2xl rounded-xl shadow-xl dark:shadow-2xl border border-slate-200/80 dark:border-slate-700/60 overflow-hidden"
            variants={tabContainerVariant}
            initial="hidden"
            animate="visible"
          >
            {/* Tab Navigation */}
            <div className="relative border-b border-slate-200/80 dark:border-slate-700/50 bg-slate-100/60 dark:bg-slate-800/40">
              <div className="flex justify-center overflow-x-auto p-2 space-x-2">
                <LayoutGroup id="depict-tabs-enhanced">
                  {tabs.map((tab) => {
                    const isActive = activeTabId === tab.id;
                    return (
                      <motion.button
                        key={tab.id}
                        onClick={() => setActiveTabId(tab.id)}
                        className={`tab-button relative flex items-center px-5 py-3 text-sm font-medium rounded-lg transition-colors duration-200 whitespace-nowrap focus:outline-none focus-visible:ring-2 focus-visible:ring-offset-2 focus-visible:ring-offset-[var(--bg-secondary)] focus-visible:ring-[var(--text-accent)] ${
                          isActive
                            ? "text-sky-700 dark:text-white"
                            : "text-slate-600 dark:text-slate-400 hover:text-slate-900 dark:hover:text-white"
                        }`}
                        aria-selected={isActive}
                        role="tab"
                        whileHover={{ scale: 1.04 }}
                        whileTap={{ scale: 0.96 }}
                        transition={{
                          type: "spring",
                          stiffness: 400,
                          damping: 15,
                        }}
                      >
                        {isActive && (
                          <motion.div
                            className="absolute inset-0 bg-white dark:bg-slate-700 rounded-lg shadow-sm"
                            layoutId="activeTabPill"
                            transition={{
                              type: "spring",
                              stiffness: 380,
                              damping: 35,
                              mass: 0.8,
                            }}
                            style={{ zIndex: 0 }}
                          />
                        )}
                        <motion.span
                          className="relative z-10 flex items-center"
                          animate={{ opacity: 1, scale: 1 }}
                          whileHover={{ scale: isActive ? 1 : 1.02 }}
                        >
                          <tab.icon
                            className={`h-5 w-5 mr-2 flex-shrink-0 transition-colors duration-200 ${
                              isActive
                                ? "text-sky-600 dark:text-sky-400"
                                : "text-slate-500 dark:text-slate-400 group-hover:text-slate-600 dark:group-hover:text-slate-300"
                            }`}
                          />
                          <span>{tab.name}</span>
                        </motion.span>
                      </motion.button>
                    );
                  })}
                </LayoutGroup>
              </div>
            </div>

            {/* Tab Content Header */}
            {activeTab && (
              <motion.div
                className="p-5 sm:p-6 border-b border-slate-200/80 dark:border-slate-700/50 bg-white/20 dark:bg-slate-800/10"
                key={`${activeTabId}-header`}
                initial={{ opacity: 0, y: -10 }}
                animate={{ opacity: 1, y: 0 }}
                transition={{ duration: 0.4, delay: 0.15 }}
              >
                {/* This content is now inside the 3/4 width container */}
                <div className="flex items-start sm:items-center max-w-5xl mx-auto">
                  <motion.div
                    initial={{ scale: 0.8 }}
                    animate={{ scale: 1 }}
                    transition={{ delay: 0.2 }}
                  >
                    <activeTab.icon className="h-7 w-7 sm:h-8 sm:w-8 text-[var(--text-accent)] mr-3 sm:mr-4 flex-shrink-0 mt-0.5 sm:mt-0" />
                  </motion.div>
                  <div>
                    <h2 className="text-xl sm:text-2xl font-bold text-[var(--text-primary)]">
                      {activeTab.name}
                    </h2>
                    <p className="text-sm text-[var(--text-secondary)] mt-1">
                      {activeTab.description}
                    </p>
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
                  className="p-5 sm:p-6"
                >
                  {/* ActiveComponent will render within the 3/4 width */}
                  {ActiveComponent && (
                    <ActiveComponent
                      isActive={activeTabId === "3d-depiction"}
                    />
                  )}
                </motion.div>
              </AnimatePresence>
            </div>
          </motion.div>{" "}
          {/* End Tab Container */}
        </div>{" "}
        {/* End Width Constraint Wrapper */}
      </div>{" "}
      {/* End Content Area */}
      {/* Global Styles */}
      <style jsx global>{`
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
        /* --- Inactive Tab Border Reveal --- */
        .tab-button {
          position: relative;
        }
        .tab-button::before {
          /* ... same as before ... */
        }
        .tab-button:not([aria-selected="true"]):hover::before {
          /* ... same as before ... */
        }
      `}</style>
    </motion.div> // End Main Container
  );
};

export default DepictPage;
