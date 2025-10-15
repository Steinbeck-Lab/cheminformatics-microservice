// Description: This page provides a tabbed interface for format conversion and coordinate generation.
import React, { useState } from "react";
import { motion, AnimatePresence, LayoutGroup } from "framer-motion";
import FormatConversionView from "../components/convert/FormatConversionView"; // Assuming path is correct
import Mol2DView from "../components/convert/Mol2DView"; // Assuming path is correct
import Mol3DView from "../components/convert/Mol3DView"; // Assuming path is correct

// Import icons for the tabs
import { HiOutlineDocumentDuplicate, HiOutlineViewGrid, HiOutlineCube } from "react-icons/hi";

// Tab data with added icons and descriptions
const tabs = [
  {
    id: "format-conversion",
    name: "Format Conversion",
    component: FormatConversionView,
    icon: HiOutlineDocumentDuplicate,
    description: "Convert between different chemical file formats",
  },
  {
    id: "2d-coordinates",
    name: "2D Coordinates",
    component: Mol2DView,
    icon: HiOutlineViewGrid,
    description: "Generate 2D coordinates for molecular structures",
  },
  {
    id: "3d-coordinates",
    name: "3D Coordinates",
    component: Mol3DView,
    icon: HiOutlineCube,
    description: "Generate 3D coordinates for molecular structures",
  },
];

// --- Animation Variants (Consistent with other pages) ---
const pageVariants = {
  hidden: { opacity: 0 },
  visible: { opacity: 1, transition: { duration: 0.5, ease: "easeOut" } },
};
const headerContainerVariants = {
  hidden: { opacity: 0 },
  visible: {
    opacity: 1,
    transition: { staggerChildren: 0.1, delayChildren: 0.1 },
  },
};
const headerItemVariants = {
  hidden: { opacity: 0, y: -10 },
  visible: {
    opacity: 1,
    y: 0,
    transition: { duration: 0.5, ease: [0.2, 0.65, 0.3, 0.9] },
  },
};
const tabContainerVariant = {
  hidden: { opacity: 0, y: 20, scale: 0.98 },
  visible: {
    opacity: 1,
    y: 0,
    scale: 1,
    transition: { duration: 0.5, delay: 0.2, ease: [0.2, 0.65, 0.3, 0.9] },
  },
};
const contentVariants = {
  hidden: { opacity: 0, x: -15 }, // Adjusted slide distance
  visible: {
    opacity: 1,
    x: 0,
    transition: { duration: 0.4, ease: [0.25, 1, 0.5, 1] },
  },
  exit: { opacity: 0, x: 15, transition: { duration: 0.25, ease: "easeIn" } },
};

const ConvertPage = () => {
  const [activeTabId, setActiveTabId] = useState(tabs[0].id);

  const activeTab = tabs.find((tab) => tab.id === activeTabId);
  const ActiveComponent = activeTab ? activeTab.component : null;

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
      {/* Content Area - Outer padding container */}
      <div className="relative w-full mx-auto px-4 sm:px-6 lg:px-8 py-12 md:py-16 z-10">
        {/* Inner wrapper for width constraint (lg:w-3/4) */}
        <div className="w-full lg:w-3/4 mx-auto">
          {/* Page Header - Animated */}
          <motion.div
            className="mb-8 md:mb-10 max-w-5xl mx-auto text-center"
            variants={headerContainerVariants}
            initial="hidden"
            animate="visible"
          >
            <motion.h1
              variants={headerItemVariants}
              className="text-3xl md:text-4xl lg:text-5xl font-bold text-[var(--text-primary)] mb-3"
            >
              Format Conversion & Coordinate Generation
            </motion.h1>
            <motion.p
              variants={headerItemVariants}
              className="text-[var(--text-secondary)] text-base md:text-lg max-w-3xl mx-auto"
            >
              Convert between chemical file formats, generate 2D and 3D coordinates.
            </motion.p>
          </motion.div>
          {/* Tab Container - Animated */}
          <motion.div
            className="glass rounded-xl shadow-lg border border-slate-200 dark:border-slate-700/50 overflow-hidden"
            variants={tabContainerVariant}
            initial="hidden"
            animate="visible"
          >
            {/* Tab Navigation */}
            <div className="relative border-b border-slate-200/80 dark:border-slate-700/50 bg-slate-100/50 dark:bg-slate-800/30">
              <div className="flex justify-center overflow-x-auto p-2 space-x-1">
                <LayoutGroup id="convert-tabs">
                  {tabs.map((tab) => {
                    const isActive = activeTabId === tab.id;
                    return (
                      <motion.button
                        key={tab.id}
                        onClick={() => setActiveTabId(tab.id)}
                        // Apply tab-button class for potential border reveal via global CSS
                        className={`tab-button relative flex items-center px-4 py-2.5 text-sm rounded-lg transition-colors whitespace-nowrap focus:outline-none focus-visible:ring-2 focus-visible:ring-offset-2 focus-visible:ring-offset-[var(--bg-secondary)] focus-visible:ring-[var(--text-accent)] ${
                          isActive
                            ? "text-sky-700 dark:text-white"
                            : "text-slate-600 dark:text-slate-400 hover:text-slate-900 dark:hover:text-white" // Adjusted inactive colors slightly
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
                        {/* Animated Pill Background */}
                        {isActive && (
                          <motion.div
                            className="absolute inset-0 bg-white dark:bg-slate-700 rounded-lg shadow-sm"
                            layoutId="activeTabPill-convert" // Unique layoutId for this tab group
                            transition={{
                              type: "spring",
                              stiffness: 380,
                              damping: 35,
                              mass: 0.8,
                            }}
                            style={{ zIndex: 0 }}
                          />
                        )}
                        {/* Icon and Text */}
                        <motion.span
                          className="relative z-10 flex items-center"
                          animate={{ opacity: 1, scale: 1 }}
                          whileHover={{ scale: isActive ? 1 : 1.02 }}
                        >
                          <tab.icon
                            className={`h-5 w-5 mr-2 flex-shrink-0 transition-colors duration-200 ${isActive ? "text-sky-600 dark:text-sky-400" : "text-slate-500 dark:text-slate-400 group-hover:text-slate-600 dark:group-hover:text-slate-300"}`}
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
                  {/* Render active component */}
                  {/* Ensure child components handle layout and theme */}
                  {ActiveComponent && <ActiveComponent />}
                </motion.div>
              </AnimatePresence>
            </div>
          </motion.div>{" "}
          {/* End Tab Container */}
        </div>{" "}
        {/* End Width Constraint Wrapper */}
      </div>{" "}
      {/* End Content Area */}
    </motion.div> // End Main Container
  );
};

export default ConvertPage;
