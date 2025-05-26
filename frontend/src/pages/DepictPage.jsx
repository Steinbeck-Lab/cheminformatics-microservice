// Description: DepictPage component for generating 2D and 3D depictions of chemical structures
import React, { useState, useEffect } from "react";
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
import StructureDrawView from "../components/depict/StructureDrawView";

// Import icons
import {
  HiOutlineViewGrid,
  HiOutlineCube,
  HiOutlineSearch,
  HiOutlinePencil,
  HiChevronDown,
} from "react-icons/hi";

// Tab data
const tabs = [
  {
    id: "batch-depiction",
    name: "2D Depiction",
    component: Depict2DMultiView,
    icon: HiOutlineViewGrid,
    description: "Generate 2D depictions for multiple molecules at once.",
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
  {
    id: "structure-draw",
    name: "Draw a Structure",
    component: StructureDrawView,
    icon: HiOutlinePencil,
    description:
      "Draw and edit chemical structures using a user-friendly interface",
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
const mobileMenuVariants = {
  closed: { opacity: 0, y: -10, height: 0 },
  open: {
    opacity: 1,
    y: 0,
    height: "auto",
    transition: { duration: 0.3, ease: "easeOut" },
  },
};

const DepictPage = () => {
  const [activeTabId, setActiveTabId] = useState(tabs[0].id);
  const [isMobileMenuOpen, setIsMobileMenuOpen] = useState(false);
  const [isMobile, setIsMobile] = useState(false);

  // Check if the window is mobile size
  useEffect(() => {
    const checkIfMobile = () => {
      setIsMobile(window.innerWidth < 768);
    };

    // Initial check
    checkIfMobile();

    // Add event listener
    window.addEventListener("resize", checkIfMobile);

    // Clean up
    return () => window.removeEventListener("resize", checkIfMobile);
  }, []);

  const activeTab = tabs.find((tab) => tab.id === activeTabId);
  const ActiveComponent = activeTab ? activeTab.component : null;

  const { scrollYProgress } = useScroll();
  const meshY = useTransform(scrollYProgress, [0, 1], ["0%", "15%"]);
  const noiseY = useTransform(scrollYProgress, [0, 1], ["0%", "30%"]);

  // Function to handle tab selection and close the mobile menu
  const handleTabSelection = (tabId) => {
    setActiveTabId(tabId);
    setIsMobileMenuOpen(false);
  };

  return (
    <motion.div
      className="relative min-h-screen w-full bg-gradient-to-b from-slate-100 to-slate-200 dark:from-gray-950 dark:to-indigo-950/60 text-slate-900 dark:text-slate-100 font-sans overflow-x-hidden isolate"
      variants={pageVariants}
      initial="hidden"
      animate="visible"
    >
      {/* --- Enhanced Background Effects --- */}
      <div className="absolute inset-0 -z-20 overflow-hidden dark:opacity-100 opacity-0 transition-opacity duration-500">
        <div className="absolute inset-0 bg-gradient-to-br from-gray-950 via-slate-900 to-indigo-950"></div>
        <motion.div
          className="animated-mesh-gradient"
          style={{ y: meshY }}
          animate={{
            filter: ["blur(40px)", "blur(60px)", "blur(40px)"],
            opacity: [0.4, 0.5, 0.4],
          }}
          transition={{
            duration: 10,
            repeat: Infinity,
            repeatType: "reverse",
          }}
        />
      </div>
      <div className="absolute inset-0 -z-20 overflow-hidden dark:opacity-0 opacity-100 transition-opacity duration-500">
        <div className="absolute inset-0 bg-gradient-to-br from-sky-50 via-white to-indigo-100"></div>
        <motion.div
          className="absolute inset-0 opacity-[0.05]"
          style={{
            y: noiseY,
            backgroundImage: `url("data:image/svg+xml,%3Csvg width='70' height='70' viewBox='0 0 70 70' xmlns='http://www.w3.org/2000/svg'%3E%3Cg fill='none' fill-rule='evenodd'%3E%3Cg fill='%23AAB' fill-opacity='.25'%3E%3Cpath d='M35 0v70M0 35h70' stroke-width='1.5'/%3E%3C/g%3E%3C/g%3E%3C/svg%3E")`,
          }}
        />
      </div>
      <motion.div
        className="absolute inset-0 -z-10 opacity-[0.03] dark:opacity-[0.04]"
        style={{ y: noiseY, backgroundImage: "url(/noise.svg)" }}
        animate={{
          backgroundPosition: ["0% 0%", "100% 100%"],
        }}
        transition={{
          duration: 120,
          repeat: Infinity,
          repeatType: "reverse",
          ease: "linear",
        }}
      ></motion.div>
      {/* Floating particles in dark mode */}
      <div className="absolute inset-0 -z-10 dark:opacity-40 opacity-0 transition-opacity duration-500 overflow-hidden">
        <div className="absolute w-2 h-2 bg-blue-400 rounded-full top-1/4 left-1/4"></div>
        <motion.div
          className="absolute w-3 h-3 bg-purple-500 rounded-full top-1/3 right-1/4"
          animate={{ y: [0, -30, 0], opacity: [0.4, 0.8, 0.4] }}
          transition={{ duration: 6, repeat: Infinity, repeatType: "reverse" }}
        ></motion.div>
        <motion.div
          className="absolute w-2 h-2 bg-indigo-400 rounded-full bottom-1/4 left-1/3"
          animate={{ y: [0, 20, 0], opacity: [0.3, 0.7, 0.3] }}
          transition={{ duration: 8, repeat: Infinity, repeatType: "reverse" }}
        ></motion.div>
        <motion.div
          className="absolute w-4 h-4 bg-sky-400 rounded-full bottom-1/3 right-1/3 opacity-20"
          animate={{ y: [0, -40, 0], opacity: [0.2, 0.5, 0.2] }}
          transition={{ duration: 10, repeat: Infinity, repeatType: "reverse" }}
        ></motion.div>
      </div>
      {/* Content Area - Outer padding container */}
      <div className="relative w-full mx-auto px-4 sm:px-6 lg:px-8 py-8 md:py-12 z-10">
        {/* Wrapper div for width constraint */}
        <div className="w-full lg:w-3/4 mx-auto">
          {/* Page Header - Animated */}
          <motion.div
            className="mb-6 md:mb-8 max-w-5xl mx-auto text-center"
            variants={headerContainerVariants}
            initial="hidden"
            animate="visible"
          >
            <motion.h1
              variants={headerItemVariants}
              className="text-2xl md:text-4xl lg:text-5xl font-bold text-[var(--text-primary)] mb-2 md:mb-3"
            >
              Chemical Structure Depiction
            </motion.h1>
            <motion.p
              variants={headerItemVariants}
              className="text-[var(--text-secondary)] text-sm md:text-lg max-w-3xl mx-auto"
            >
              Generate customizable 2D and interactive 3D visualizations of
              chemical structures.
            </motion.p>
          </motion.div>
          {/* Tab Container - Enhanced visual design */}
          <motion.div
            className="glass bg-white/90 dark:bg-slate-800/90 backdrop-blur-xl dark:backdrop-blur-2xl rounded-xl shadow-2xl dark:shadow-2xl border border-slate-200/80 dark:border-slate-700/60 overflow-hidden"
            variants={tabContainerVariant}
            initial="hidden"
            animate="visible"
            whileHover={{ boxShadow: "0 25px 50px -12px rgba(0, 0, 0, 0.25)" }}
            transition={{ duration: 0.5 }}
          >
            {/* Tab Navigation - Mobile Optimized */}
            <div className="relative border-b border-slate-200/80 dark:border-slate-700/50 bg-gradient-to-r from-slate-100/80 to-slate-200/80 dark:from-slate-800/60 dark:to-slate-900/60">
              {/* Mobile Tab Navigation (Dropdown Style) */}
              <div className="block md:hidden">
                <button
                  onClick={() => setIsMobileMenuOpen(!isMobileMenuOpen)}
                  className="flex items-center justify-between w-full py-3 px-4 text-left bg-white/50 dark:bg-slate-800/50 focus:outline-none"
                  aria-expanded={isMobileMenuOpen}
                >
                  <div className="flex items-center">
                    <activeTab.icon className="h-5 w-5 mr-2 text-sky-600 dark:text-sky-400" />
                    <span className="font-medium text-sky-700 dark:text-white">
                      {activeTab.name}
                    </span>
                  </div>
                  <motion.span
                    animate={{ rotate: isMobileMenuOpen ? 180 : 0 }}
                    transition={{ duration: 0.3 }}
                  >
                    <HiChevronDown className="h-5 w-5 text-slate-500 dark:text-slate-400" />
                  </motion.span>
                </button>

                <AnimatePresence>
                  {isMobileMenuOpen && (
                    <motion.div
                      variants={mobileMenuVariants}
                      initial="closed"
                      animate="open"
                      exit="closed"
                      className="border-t border-slate-200 dark:border-slate-700/50 bg-white/90 dark:bg-slate-800/90 overflow-hidden"
                    >
                      {tabs.map((tab) => (
                        <button
                          key={tab.id}
                          onClick={() => handleTabSelection(tab.id)}
                          className={`w-full flex items-center py-3 px-4 ${
                            activeTabId === tab.id
                              ? "bg-blue-50 dark:bg-slate-700/50 text-sky-700 dark:text-white"
                              : "text-slate-600 dark:text-slate-400 hover:bg-slate-100 dark:hover:bg-slate-700/30"
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
                  <LayoutGroup id="depict-tabs-enhanced">
                    {tabs.map((tab) => {
                      const isActive = activeTabId === tab.id;
                      return (
                        <motion.button
                          key={tab.id}
                          onClick={() => setActiveTabId(tab.id)}
                          className={`tab-button relative flex items-center px-5 py-2.5 text-sm font-medium rounded-lg transition-colors duration-200 whitespace-nowrap focus:outline-none focus-visible:ring-2 focus-visible:ring-offset-2 focus-visible:ring-offset-[var(--bg-secondary)] focus-visible:ring-[var(--text-accent)] ${
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
                              className="absolute inset-0 bg-white dark:bg-slate-700/90 rounded-lg shadow-sm"
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
            </div>

            {/* Tab Content Header - Adaptive for mobile */}
            {activeTab && (
              <motion.div
                key={`${activeTabId}-header`}
                initial={{ opacity: 0, y: -10 }}
                animate={{ opacity: 1, y: 0 }}
                transition={{ duration: 0.4, delay: 0.15 }}
                className="mx-3 sm:mx-4 mt-3 sm:mt-4 mb-2"
              >
                <div className="bg-gradient-to-r from-blue-500 to-purple-600 dark:from-blue-800 dark:to-purple-950 rounded-lg shadow-lg overflow-hidden relative">
                  {/* Animated background elements with theme-adaptive colors */}
                  <div className="absolute inset-0 overflow-hidden opacity-10">
                    <div className="absolute left-0 top-0 w-32 h-32 rounded-full bg-white dark:bg-slate-300 transform -translate-x-1/2 -translate-y-1/2"></div>
                    <div className="absolute right-16 bottom-0 w-48 h-48 rounded-full bg-white dark:bg-slate-300 transform translate-x-1/2 translate-y-1/2"></div>
                    <div className="absolute left-1/3 top-1/2 w-24 h-24 rounded-full bg-white dark:bg-slate-300 transform -translate-y-1/2"></div>
                  </div>

                  <div className="relative p-4 md:p-6 flex flex-col md:flex-row items-center justify-between gap-3 md:gap-4">
                    {/* Title and description with improved text contrast */}
                    <div className="space-y-1 text-center md:text-left max-w-3xl">
                      <h1 className="text-xl md:text-3xl font-bold leading-tight text-white dark:text-slate-50">
                        {activeTab.name}
                      </h1>
                      <p className="text-blue-50 dark:text-blue-100 text-xs md:text-base opacity-90">
                        {activeTab.description}
                      </p>
                    </div>

                    {/* Decorative icon with adaptive background */}
                    <div className="hidden md:flex items-center justify-center bg-white/15 dark:bg-slate-200/10 backdrop-blur-sm p-4 rounded-full w-16 h-16 border border-white/20 dark:border-slate-300/20 shadow-lg">
                      <activeTab.icon className="w-8 h-8 text-white dark:text-slate-50" />
                    </div>
                  </div>
                </div>
              </motion.div>
            )}

            {/* Tab Content Area - Optimized for mobile */}
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
                  {/* Decorative elements with theme-adaptive colors */}
                  <div className="absolute -top-20 -right-20 w-40 h-40 bg-blue-100/30 dark:bg-blue-700/20 rounded-full filter blur-3xl"></div>
                  <div className="absolute -bottom-20 -left-20 w-60 h-60 bg-purple-100/30 dark:bg-purple-700/20 rounded-full filter blur-3xl"></div>

                  {/* ActiveComponent rendered here */}
                  {ActiveComponent && (
                    <ActiveComponent
                      isActive={activeTabId === "3d-depiction"}
                      isMobile={isMobile}
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
        :root {
          --text-primary: #1e293b;
          --text-secondary: #64748b;
          --text-accent: #0284c7;
          --bg-primary: #ffffff;
          --bg-secondary: #f1f5f9;
        }

        .dark {
          --text-primary: #f1f5f9;
          --text-secondary: #94a3b8;
          --text-accent: #38bdf8;
          --bg-primary: #0f172a;
          --bg-secondary: #1e293b;
        }

        /* Glass effect */
        .glass {
          backdrop-filter: blur(16px);
          -webkit-backdrop-filter: blur(16px);
        }

        /* Animated Mesh Gradient */
        @keyframes mesh-gradient-move {
          0% {
            background-position: 0% 50%;
          }
          50% {
            background-position: 100% 50%;
          }
          100% {
            background-position: 0% 50%;
          }
        }

        .animated-mesh-gradient {
          position: absolute;
          inset: -100%;
          background-image: radial-gradient(
              circle at 25% 25%,
              rgba(60, 90, 180, 0.4) 0%,
              transparent 50%
            ),
            radial-gradient(
              circle at 75% 75%,
              rgba(120, 50, 200, 0.4) 0%,
              transparent 50%
            ),
            radial-gradient(
              circle at 85% 15%,
              rgba(90, 40, 180, 0.4) 0%,
              transparent 50%
            );
          filter: blur(60px);
          opacity: 0.5;
          transform-origin: center;
          animation: mesh-gradient-move 30s ease infinite;
          overflow: hidden;
        }

        .dark .animated-mesh-gradient {
          background-image: radial-gradient(
              circle at 25% 25%,
              rgba(75, 100, 255, 0.5) 0%,
              transparent 50%
            ),
            radial-gradient(
              circle at 75% 75%,
              rgba(160, 70, 240, 0.5) 0%,
              transparent 50%
            ),
            radial-gradient(
              circle at 85% 15%,
              rgba(100, 60, 200, 0.5) 0%,
              transparent 50%
            );
          filter: blur(50px);
          opacity: 0.6;
        }

        /* Tab button styling */
        .tab-button {
          position: relative;
        }

        .tab-button::before {
          content: "";
          position: absolute;
          bottom: -2px;
          left: 0;
          right: 0;
          height: 2px;
          background: linear-gradient(to right, #3b82f6, #8b5cf6);
          opacity: 0;
          transition: opacity 0.2s ease;
        }

        .tab-button:not([aria-selected="true"]):hover::before {
          opacity: 0.5;
        }

        /* Active tab indicator for mobile */
        .mobile-tab-active {
          position: relative;
        }

        .mobile-tab-active::after {
          content: "";
          position: absolute;
          left: 0;
          top: 0;
          bottom: 0;
          width: 3px;
          background: linear-gradient(to bottom, #3b82f6, #8b5cf6);
        }

        /* Fade in animation for content */
        @keyframes fadeIn {
          from {
            opacity: 0;
          }
          to {
            opacity: 1;
          }
        }

        .fadeIn {
          animation: fadeIn 0.5s ease forwards;
        }

        /* Mobile optimizations */
        @media (max-width: 640px) {
          .tab-button {
            padding-left: 0.75rem;
            padding-right: 0.75rem;
          }
        }

        /* Improve touch targets */
        @media (max-width: 767px) {
          button,
          [role="button"] {
            min-height: 44px;
            min-width: 44px;
          }
        }
      `}</style>
    </motion.div> // End Main Container
  );
};

export default DepictPage;
