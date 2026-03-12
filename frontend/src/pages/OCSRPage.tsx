// Description: This page implements the Optical Chemical Structure Recognition (OCSR) tool using DECIMER.
import React from "react";
import { motion } from "framer-motion";
import OCRView from "../components/ocsr/OCRView"; // Assuming path is correct
import { HiOutlineDocumentSearch, HiOutlineExternalLink } from "react-icons/hi";

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
const contentContainerVariant = {
  // Similar to tabContainerVariant
  hidden: { opacity: 0, y: 20, scale: 0.98 },
  visible: {
    opacity: 1,
    y: 0,
    scale: 1,
    transition: { duration: 0.5, delay: 0.2, ease: [0.2, 0.65, 0.3, 0.9] },
  },
};
const buttonVariant = {
  hidden: { opacity: 0, scale: 0.95 },
  visible: {
    opacity: 1,
    scale: 1,
    transition: { duration: 0.5, ease: [0.2, 0.65, 0.3, 0.9] },
  },
  hover: { scale: 1.03, transition: { duration: 0.2 } },
};

const OCSRPage = () => {
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
              Optical Chemical Structure Recognition
            </motion.h1>
            <motion.p
              variants={headerItemVariants}
              className="text-[var(--text-secondary)] text-base md:text-lg max-w-3xl mx-auto"
            >
              Extract chemical structures from images using deep learning (DECIMER & MARCUS).
            </motion.p>
          </motion.div>
          {/* DECIMER & MARCUS Buttons */}
          <motion.div
            className="flex flex-col md:flex-row gap-4 md:gap-6 justify-center items-center mb-8 md:mb-10"
            variants={headerContainerVariants}
            initial="hidden"
            animate="visible"
          >
            <motion.a
              href="https://decimer.ai"
              target="_blank"
              rel="noopener noreferrer"
              className="flex items-center justify-center px-6 py-4 bg-gradient-to-r from-blue-600 to-blue-700 hover:from-blue-500 hover:to-blue-600 text-white font-semibold rounded-lg shadow-lg w-full md:w-auto"
              variants={buttonVariant}
              whileHover="hover"
            >
              <span className="mr-2 text-lg">Visit DECIMER</span>
              <HiOutlineExternalLink className="h-5 w-5" />
            </motion.a>
            <motion.a
              href="https://marcus.decimer.ai"
              target="_blank"
              rel="noopener noreferrer"
              className="flex items-center justify-center px-6 py-4 bg-gradient-to-r from-purple-600 to-purple-700 hover:from-purple-500 hover:to-purple-600 text-white font-semibold rounded-lg shadow-lg w-full md:w-auto"
              variants={buttonVariant}
              whileHover="hover"
            >
              <span className="mr-2 text-lg">Visit MARCUS</span>
              <HiOutlineExternalLink className="h-5 w-5" />
            </motion.a>
          </motion.div>
          {/* DECIMER & MARCUS Info Cards */}
          <motion.div
            className="grid grid-cols-1 md:grid-cols-2 gap-4 md:gap-6 mb-8 md:mb-10"
            variants={headerContainerVariants}
            initial="hidden"
            animate="visible"
          >
            <motion.div
              className="bg-white dark:bg-slate-800/80 rounded-lg p-5 border border-slate-200 dark:border-slate-700/50 shadow-md"
              variants={headerItemVariants}
            >
              <h3 className="text-lg font-semibold text-blue-600 dark:text-blue-400 mb-2">
                What is DECIMER?
              </h3>
              <p className="text-[var(--text-secondary)] mb-3">
                DECIMER (Deep lEarning for Chemical IMagE Recognition) is an advanced deep learning
                model developed by the Steinbeck Lab for extracting chemical structures from images.
              </p>
              <p className="text-[var(--text-secondary)]">
                It can process images from scientific publications, patents, and other documents,
                converting visual chemical structures into machine-readable formats like SMILES.
              </p>
            </motion.div>
            <motion.div
              className="bg-white dark:bg-slate-800/80 rounded-lg p-5 border border-slate-200 dark:border-slate-700/50 shadow-md"
              variants={headerItemVariants}
            >
              <h3 className="text-lg font-semibold text-purple-600 dark:text-purple-400 mb-2">
                What is MARCUS?
              </h3>
              <p className="text-[var(--text-secondary)] mb-3">
                MARCUS (Molecular Annotation and Recognition for Curating Unravelled Structures) is
                an integrated web-based platform designed for natural product literature curation.
                It combines automated text annotation, multi-engine OCSR, and direct database
                submission capabilities.
              </p>
              <p className="text-[var(--text-secondary)]">
                MARCUS employs a Human-in-the-loop ensemble approach integrating DECIMER, MolNexTR,
                and MolScribe for structure recognition, significantly streamlining the workflow
                from PDF upload to database submission.
              </p>
            </motion.div>
          </motion.div>
          {/* Main Content Container - Animated */}
          <motion.div
            className="glass rounded-xl shadow-lg border border-slate-200 dark:border-slate-700/50 overflow-hidden"
            variants={contentContainerVariant}
            initial="hidden"
            animate="visible"
          >
            {/* Content Header */}
            <div className="p-5 sm:p-6 border-b border-slate-200/80 dark:border-slate-700/50 bg-white/30 dark:bg-slate-800/20">
              <div className="flex items-start sm:items-center">
                <HiOutlineDocumentSearch className="h-7 w-7 sm:h-8 sm:w-8 text-[var(--text-accent)] mr-3 sm:mr-4 flex-shrink-0 mt-0.5 sm:mt-0" />
                <div>
                  <h2 className="text-xl sm:text-2xl font-bold text-[var(--text-primary)]">
                    OCSR Tool
                  </h2>
                  <p className="text-sm text-[var(--text-secondary)] mt-1">
                    Analyze images to extract chemical structures
                  </p>
                </div>
              </div>
            </div>
            {/* Main Content Body with Grid */}
            <motion.div
              className="p-5 sm:p-6"
              variants={headerContainerVariants} // Re-use for simple stagger of grid children
              initial="hidden"
              animate="visible"
            >
              <div className="grid grid-cols-1 lg:grid-cols-3 gap-6">
                {/* Main OCSR Tool */}
                <motion.div
                  className="lg:col-span-2"
                  variants={headerItemVariants} // Use headerItemVariants for grid item entrance
                >
                  <OCRView />
                </motion.div>

                {/* Info Sidebar */}
                <motion.div
                  className="lg:col-span-1 space-y-6"
                  variants={headerItemVariants} // Use headerItemVariants for grid item entrance
                >
                  {/* Info Box 1 */}
                  <div className="bg-slate-100 dark:bg-slate-800/60 rounded-lg p-5 border border-slate-200 dark:border-slate-700">
                    <h3 className="text-lg font-semibold text-[var(--text-primary)] mb-3">
                      About OCSR Technology
                    </h3>
                    <p className="text-[var(--text-secondary)] text-sm leading-relaxed space-y-3">
                      <span className="block mb-3">
                        The OCSR functionality is powered by DECIMER (Deep lEarning for Chemical
                        IMagE Recognition), a deep learning model developed by the Steinbeck Lab. It
                        enables the extraction of chemical structures from images such as scientific
                        publications, patents, and other documents.
                      </span>
                      <span className="block mb-4">
                        For optimal results, provide clear images of chemical structures with good
                        contrast. The system works best with standard 2D chemical depictions.
                      </span>
                    </p>

                    {/* References Section */}
                    <div className="mt-4 pt-4 border-t border-slate-300 dark:border-slate-600">
                      <h4 className="text-sm font-semibold text-[var(--text-primary)] mb-3">
                        References
                      </h4>
                      <div className="space-y-3 text-xs">
                        <div className="border-l-2 border-blue-500 pl-3 py-1">
                          <p className="text-[var(--text-secondary)] leading-relaxed">
                            <span className="font-medium">[1]</span> Rajan, K., Brinkhaus, H.O.,
                            Agea, I.A., Zielesny, A., Steinbeck, C. (2023). DECIMER.ai: an open
                            platform for automated optical chemical structure identification,
                            segmentation and recognition in scientific publications.
                            <em className="text-blue-600 dark:text-blue-400">
                              {" "}
                              Nat Commun, 14, 5045
                            </em>
                            .
                            <a
                              href="https://doi.org/10.1038/s41467-023-40782-0"
                              target="_blank"
                              rel="noopener noreferrer"
                              className="text-blue-600 dark:text-blue-400 hover:underline ml-1 text-xs"
                            >
                              doi.org/10.1038/s41467-023-40782-0
                            </a>
                          </p>
                        </div>
                        <div className="border-l-2 border-purple-500 pl-3 py-1">
                          <p className="text-[var(--text-secondary)] leading-relaxed">
                            <span className="font-medium">[2]</span> Rajan K, Weissenborn VK,
                            Lederer L, Zielesny A, Steinbeck C (2025). MARCUS: Molecular annotation
                            and recognition for curating unravelled structures.
                            <em className="text-purple-600 dark:text-purple-400">
                              {" "}
                              Digit Discovery
                            </em>
                            .
                            <a
                              href="https://doi.org/10.1039/d5dd00313j"
                              target="_blank"
                              rel="noopener noreferrer"
                              className="text-purple-600 dark:text-purple-400 hover:underline ml-1 text-xs"
                            >
                              doi.org/10.1039/d5dd00313j
                            </a>
                          </p>
                        </div>
                      </div>
                    </div>
                  </div>
                </motion.div>
              </div>
            </motion.div>{" "}
            {/* End Main Content Body */}
          </motion.div>{" "}
          {/* End Main Content Container */}
        </div>{" "}
        {/* End Width Constraint Wrapper */}
      </div>{" "}
      {/* End Outer Padding Container */}
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
      `}</style>
    </motion.div> // End Page Container
  );
};

export default OCSRPage;
