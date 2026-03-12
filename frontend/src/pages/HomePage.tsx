// Description: HomePage component with a modern design, featuring a hero section, features grid, and recent molecules section. It uses Framer Motion for animations and Tailwind CSS for styling.
import React, { useState } from "react";
import { Link } from "react-router-dom";
import { motion } from "framer-motion";
import {
  HiOutlineArrowRight,
  HiOutlineBeaker,
  HiOutlineArrowsRightLeft,
  HiOutlineEye,
  HiOutlineWrenchScrewdriver,
  HiOutlineCamera,
  HiOutlineCubeTransparent,
} from "react-icons/hi2";
import { FaBook, FaCode } from "react-icons/fa";
import { useAppContext } from "../context/AppContext";
import MoleculeCard from "../components/common/MoleculeCard";

// --- Feature Configuration ---
const features = [
  {
    title: "Chemical Analysis",
    description:
      "Analyze molecules with descriptors, stereoisomer generation, and similarity calculations.",
    icon: HiOutlineBeaker,
    link: "/chem",
  },
  {
    title: "Format Conversion",
    description:
      "Convert between different chemical file formats (SMILES, InChI, 2D/3D coordinates).",
    icon: HiOutlineArrowsRightLeft,
    link: "/convert",
  },
  {
    title: "Visualization",
    description: "Generate customizable 2D and 3D visualizations of chemical structures.",
    icon: HiOutlineEye,
    link: "/depict",
  },
  {
    title: "Structure Tools",
    description: "Specialized tools for structure generation, sugar removal, and more.",
    icon: HiOutlineWrenchScrewdriver,
    link: "/tools",
  },
  {
    title: "OCSR",
    description: "Optical Chemical Structure Recognition to extract structures from images.",
    icon: HiOutlineCamera,
    link: "/ocsr",
  },
];

// --- Animation Variants ---
const itemVariants = {
  hidden: { y: 25, opacity: 0, scale: 0.97 },
  visible: {
    y: 0,
    opacity: 1,
    scale: 1,
    transition: { duration: 0.6, ease: [0.25, 1, 0.5, 1] },
  },
};

const containerVariants = {
  hidden: { opacity: 0 },
  visible: {
    opacity: 1,
    transition: { staggerChildren: 0.08, delayChildren: 0.15 },
  },
};

const heroTextVariants = {
  hidden: { opacity: 0, y: -20 },
  visible: (i = 1) => ({
    opacity: 1,
    y: 0,
    transition: { duration: 0.8, delay: i * 0.07, ease: [0.2, 0.8, 0.2, 1] },
  }),
};

const sectionFadeIn = {
  hidden: { opacity: 0, y: 20 },
  visible: { opacity: 1, y: 0, transition: { duration: 0.9, ease: "easeOut" } },
};

const ribbonContentVariants = {
  hidden: { opacity: 0, y: 20 },
  visible: (i = 1) => ({
    opacity: 1,
    y: 0,
    transition: { duration: 0.7, delay: i * 0.1, ease: [0.2, 0.65, 0.3, 0.9] },
  }),
};

// --- Helper Component for 3D Tilt Effect (Subtle) ---
const TiltCard = ({ children, className = "", tiltIntensity = 7 }) => {
  const [rotate, setRotate] = useState({ x: 0, y: 0 });

  const onMouseMove = (e) => {
    if ("ontouchstart" in window && e.nativeEvent instanceof TouchEvent) return;
    const card = e.currentTarget;
    const { width, height, left, top } = card.getBoundingClientRect();
    const x = e.clientX - left;
    const y = e.clientY - top;
    const rotateX = (y / height - 0.5) * -tiltIntensity;
    const rotateY = (x / width - 0.5) * tiltIntensity;
    setRotate({ x: rotateX, y: rotateY });
  };

  const onMouseLeave = () => {
    setRotate({ x: 0, y: 0 });
  };

  return (
    <motion.div
      className={`transform-style-3d transition-transform duration-700 ease-out ${className}`}
      style={{
        transform: `perspective(1800px) rotateX(${rotate.x}deg) rotateY(${rotate.y}deg)`,
        willChange: "transform",
      }}
      onMouseMove={onMouseMove}
      onMouseLeave={onMouseLeave}
      whileHover={{
        scale: 1.035,
        transition: { type: "spring", stiffness: 200, damping: 20 },
      }}
    >
      {children}
    </motion.div>
  );
};

// --- HomePage Component ---
const HomePage = () => {
  const { recentMolecules } = useAppContext();

  return (
    // Main container
    <div className="relative min-h-screen w-full bg-slate-100 dark:bg-gray-950 text-slate-900 dark:text-slate-100 font-sans overflow-x-hidden isolate">
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
      {/* --- Content Container --- */}
      {/* Outer container for padding ONLY */}
      <div className="relative w-full mx-auto px-4 sm:px-6 lg:px-8 py-12 md:py-16 z-10">
        {/* Inner wrapper - using full width */}
        <div className="w-full mx-auto">
          {/* --- Hero Section (Constrained) --- */}
          <motion.section
            className="relative text-center mb-24 md:mb-40 pt-8 sm:pt-12 md:pt-16 pb-10"
            initial="hidden"
            animate="visible"
          >
            {/* ... Glows ... */}
            <motion.div
              className="absolute -top-20 -left-20 w-80 h-80 sm:w-96 sm:h-96 bg-gradient-radial from-sky-400/20 via-transparent dark:from-blue-600/20 dark:via-transparent to-transparent rounded-full blur-3xl opacity-80 dark:opacity-60 pointer-events-none"
              style={{ x: "-10%", y: "-10%" }}
              whileInView={{
                x: "0%",
                y: "0%",
                scale: [1, 1.05, 1],
                transition: { duration: 1.8, ease: "easeInOut" },
              }}
              viewport={{ amount: 0.3 }}
            />
            <motion.div
              className="absolute -bottom-20 -right-20 w-72 h-72 sm:w-80 sm:h-80 bg-gradient-radial from-indigo-400/20 via-transparent dark:from-cyan-600/20 dark:via-transparent to-transparent rounded-full blur-3xl opacity-70 dark:opacity-50 pointer-events-none"
              style={{ x: "10%", y: "10%" }}
              whileInView={{
                x: "0%",
                y: "0%",
                scale: [1, 1.03, 1],
                transition: { duration: 1.8, ease: "easeInOut", delay: 0.2 },
              }}
              viewport={{ amount: 0.3 }}
            />

            {/* Headline */}
            <motion.h1
              variants={heroTextVariants}
              custom={0}
              className="text-4xl sm:text-6xl md:text-7xl lg:text-8xl font-extrabold text-slate-900 dark:text-white mb-5 sm:mb-6 leading-tight tracking-tighter"
            >
              Seamless Access to
              <motion.span
                variants={heroTextVariants}
                custom={1}
                className="block text-transparent bg-clip-text bg-gradient-to-r from-sky-600 via-cyan-500 to-indigo-600 dark:from-sky-400 dark:via-cyan-400 dark:to-emerald-400 mt-1 md:mt-2"
              >
                Open Cheminformatics Tools
              </motion.span>
            </motion.h1>

            {/* Sub-headline */}
            <motion.p
              variants={heroTextVariants}
              custom={2}
              className="text-base sm:text-lg md:text-xl text-slate-600 dark:text-slate-300 max-w-4xl mx-auto mb-10 sm:mb-12"
            >
              Visualize, analyze, and manipulate chemical structures seamlessly using the powerful
              Cheminformatics Microservice API.
            </motion.p>

            {/* Button Container */}
            <motion.div
              variants={heroTextVariants}
              custom={3}
              className="flex flex-col sm:flex-row flex-wrap justify-center items-center gap-4 md:gap-6"
            >
              {/* Primary Button */}
              <Link
                to="/depict"
                className="group primary-button w-full sm:w-auto relative inline-flex items-center justify-center px-7 py-3 sm:px-9 sm:py-4 text-base sm:text-lg font-semibold text-white bg-gradient-to-r from-sky-600 to-cyan-500 dark:from-blue-600 dark:to-cyan-500 rounded-xl shadow-lg hover:shadow-xl hover:brightness-110 dark:hover:shadow-cyan-500/30 transition-all duration-300 ease-out transform hover:scale-[1.06] focus:outline-none focus:ring-4 ring-offset-2 ring-offset-slate-100 dark:ring-offset-gray-950 focus:ring-cyan-500/50 dark:focus:ring-cyan-400/50"
              >
                <span className="absolute inset-0 rounded-xl bg-gradient-to-r from-white/10 to-transparent opacity-0 group-hover:opacity-100 dark:from-white/5 transition-opacity duration-300 blur-sm"></span>
                <span className="relative z-10 flex items-center">
                  Get Started{" "}
                  <HiOutlineArrowRight className="ml-2.5 h-5 w-5 transition-transform duration-300 group-hover:translate-x-1.5" />
                </span>
              </Link>
              {/* Secondary Buttons */}
              <a
                href="https://docs.api.naturalproducts.net/introduction.html"
                target="_blank"
                rel="noopener noreferrer"
                className="group w-full sm:w-auto relative inline-flex items-center justify-center px-7 py-3 sm:px-9 sm:py-4 text-base sm:text-lg font-semibold text-slate-700 dark:text-slate-200 bg-white/90 dark:bg-slate-800/80 backdrop-blur-md border border-slate-300 dark:border-slate-700 hover:bg-white dark:hover:bg-slate-700/90 hover:border-slate-400 dark:hover:border-slate-500 rounded-xl shadow-md hover:shadow-lg dark:hover:shadow-slate-700/40 transition-all duration-300 ease-out transform hover:scale-[1.06] focus:outline-none focus:ring-4 focus:ring-slate-400/50 dark:focus:ring-slate-600/50 ring-offset-2 ring-offset-slate-100 dark:ring-offset-gray-950"
              >
                <span className="absolute inset-0 rounded-xl bg-gradient-to-t from-black/5 to-transparent opacity-0 group-hover:opacity-50 dark:from-white/5 dark:group-hover:opacity-100 transition-opacity duration-300"></span>
                <span className="relative z-10 flex items-center">
                  <FaCode className="mr-2 h-5 w-5 text-green-600 dark:text-green-400" /> Guides
                </span>
              </a>
              <a
                href="https://api.naturalproducts.net/latest/docs"
                target="_blank"
                rel="noopener noreferrer"
                className="group w-full sm:w-auto relative inline-flex items-center justify-center px-7 py-3 sm:px-9 sm:py-4 text-base sm:text-lg font-semibold text-slate-700 dark:text-slate-200 bg-white/90 dark:bg-slate-800/80 backdrop-blur-md border border-slate-300 dark:border-slate-700 hover:bg-white dark:hover:bg-slate-700/90 hover:border-slate-400 dark:hover:border-slate-500 rounded-xl shadow-md hover:shadow-lg dark:hover:shadow-slate-700/40 transition-all duration-300 ease-out transform hover:scale-[1.06] focus:outline-none focus:ring-4 focus:ring-slate-400/50 dark:focus:ring-slate-600/50 ring-offset-2 ring-offset-slate-100 dark:ring-offset-gray-950"
              >
                <span className="absolute inset-0 rounded-xl bg-gradient-to-t from-black/5 to-transparent opacity-0 group-hover:opacity-50 dark:from-white/5 dark:group-hover:opacity-100 transition-opacity duration-300"></span>
                <span className="relative z-10 flex items-center">
                  <FaBook className="mr-2 h-5 w-5 text-blue-600 dark:text-blue-400" /> API Docs
                </span>
              </a>
            </motion.div>
          </motion.section>

          {/* --- Features Grid (Constrained) --- */}
          <motion.section
            className="mb-24 md:mb-40"
            variants={containerVariants}
            initial="hidden"
            whileInView="visible"
            viewport={{ once: true, amount: 0.1 }}
          >
            <h2 className="text-3xl sm:text-4xl md:text-5xl font-bold text-center mb-12 md:mb-16 tracking-tight text-slate-900 dark:text-white">
              Core
              <span className="text-transparent bg-clip-text bg-gradient-to-r from-sky-600 to-indigo-600 dark:from-blue-400 dark:to-cyan-400">
                {" "}
                Capabilities
              </span>
            </h2>
            <div className="grid grid-cols-1 sm:grid-cols-2 md:grid-cols-3 lg:grid-cols-5 gap-1.0">
              {/* ... features.map ... */}
              {features.map((feature, index) => {
                const Icon = feature.icon;
                return (
                  <motion.div key={index} variants={itemVariants} className="mx-auto w-[85%]">
                    <TiltCard tiltIntensity={7}>
                      <Link
                        to={feature.link}
                        className="feature-card-enhanced group relative flex flex-col h-full min-h-[270px] sm:min-h-[290px] bg-white/80 dark:bg-slate-800/75 backdrop-blur-xl border border-transparent rounded-3xl p-4 sm:p-6 transition-all duration-400 ease-out shadow-lg dark:shadow-xl hover:shadow-xl dark:hover:shadow-blue-900/30 overflow-hidden transform-style-3d"
                        style={{
                          transform: "translateZ(0)",
                          willChange: "border-color, backdrop-filter",
                        }}
                      >
                        <div className="absolute inset-0 rounded-3xl bg-gradient-radial from-[var(--primary-accent-faint)] via-transparent to-transparent opacity-0 group-hover:opacity-60 dark:group-hover:opacity-100 transition-opacity duration-500 pointer-events-none"></div>
                        <motion.div
                          className="mb-5 sm:mb-6 flex justify-center items-center h-14 w-14 sm:h-16 sm:w-16 rounded-2xl bg-gradient-to-br from-slate-100 to-sky-100 dark:from-blue-700/60 dark:to-cyan-700/60 group-hover:from-sky-100 group-hover:to-cyan-100 dark:group-hover:from-blue-600/70 dark:group-hover:to-cyan-600/70 transition-all duration-300 shadow-md dark:shadow-inner transform group-hover:scale-105"
                          style={{ transform: "translateZ(40px)" }}
                        >
                          <Icon className="h-7 w-7 sm:h-8 sm:w-8 text-sky-700 dark:text-sky-200 group-hover:text-cyan-800 dark:group-hover:text-cyan-100 transition-colors duration-300" />
                        </motion.div>
                        <div style={{ transform: "translateZ(25px)" }}>
                          <h3 className="text-lg sm:text-xl font-semibold text-slate-800 dark:text-white mb-2 sm:mb-3">
                            {feature.title}
                          </h3>
                          <p className="text-sm sm:text-[15px] text-slate-600 dark:text-slate-300 flex-grow mb-5 sm:mb-6 leading-relaxed">
                            {feature.description}
                          </p>
                        </div>
                        <div
                          className="mt-auto flex items-center text-sky-700 dark:text-sky-400 group-hover:text-cyan-800 dark:group-hover:text-cyan-200 font-medium transition-colors duration-300 text-sm"
                          style={{ transform: "translateZ(15px)" }}
                        >
                          <span>Explore Feature</span>
                          <HiOutlineArrowRight className="ml-2 h-4 w-4 transition-transform duration-300 group-hover:translate-x-1.5" />
                        </div>
                      </Link>
                    </TiltCard>
                  </motion.div>
                );
              })}
            </div>
          </motion.section>

          {/* --- Recent Molecules Section (Constrained) --- */}
          {recentMolecules && recentMolecules.length > 0 && (
            <motion.section
              className="mb-24 md:mb-40"
              variants={sectionFadeIn}
              initial="hidden"
              whileInView="visible"
              viewport={{ once: true, amount: 0.15 }}
            >
              <div className="flex flex-col sm:flex-row justify-between items-start sm:items-center mb-10 sm:mb-12 gap-3 sm:gap-6">
                <h2 className="text-2xl sm:text-3xl md:text-4xl font-bold text-slate-900 dark:text-white tracking-tight">
                  Recently Viewed
                </h2>
                <Link
                  to="/chem"
                  className="text-sky-700 dark:text-sky-400 hover:text-indigo-700 dark:hover:text-cyan-200 flex-shrink-0 flex items-center font-medium transition-colors duration-300 group text-sm hover-underline-link"
                >
                  <span>View all</span>
                  <HiOutlineArrowRight className="ml-1.5 h-4 w-4 transition-transform duration-300 group-hover:translate-x-1" />
                </Link>
              </div>
              <motion.div
                className="grid grid-cols-1 sm:grid-cols-2 md:grid-cols-3 lg:grid-cols-5 gap-0.5"
                variants={containerVariants}
                initial="hidden"
                whileInView="visible"
                viewport={{ once: true, amount: 0.2 }}
              >
                {/* ... recentMolecules.map ... */}
                {recentMolecules.slice(0, 3).map((molecule, index) => (
                  <motion.div key={index} variants={itemVariants} className="mx-auto w-[85%]">
                    <TiltCard tiltIntensity={5}>
                      <div className="relative bg-white/80 dark:bg-slate-800/70 backdrop-blur-lg sm:backdrop-blur-xl rounded-2xl shadow-lg dark:shadow-xl border border-slate-200 dark:border-slate-700/80 overflow-hidden transition-all duration-300 hover:border-slate-300 dark:hover:border-slate-600 hover:bg-slate-50 dark:hover:bg-slate-700/80 p-1">
                        <div className="absolute inset-0 rounded-2xl bg-gradient-radial from-[var(--primary-accent-fainter)] via-transparent to-transparent opacity-0 hover:opacity-70 transition-opacity duration-500 pointer-events-none"></div>
                        <MoleculeCard
                          smiles={molecule.smiles}
                          title={molecule.name || `Molecule ${index + 1}`}
                          size="sm"
                          className="p-4 text-slate-800 dark:text-white bg-transparent"
                        />
                      </div>
                    </TiltCard>
                  </motion.div>
                ))}
              </motion.div>
            </motion.section>
          )}
        </div>{" "}
        {/* End Width Constraint Wrapper */}
        {/*
          The API Information Section, previously here (lines approx. 363-409),
          is moved out of this container to achieve full-width background.
        */}
      </div>{" "}
      {/* End Outer Padding Container */}
      {/* --- API Information Section (Ribbon Style - Full Width Background) --- */}
      {/* This section's background spans full width; content is centered by its own container. */}
      <motion.section
        className="relative bg-gradient-to-r from-sky-100 via-indigo-100 to-fuchsia-100 dark:from-gray-900 dark:via-indigo-950/80 dark:to-fuchsia-950/70 py-12 sm:py-16 md:py-20"
        variants={sectionFadeIn}
        initial="hidden"
        whileInView="visible"
        viewport={{ once: true, amount: 0.2 }}
      >
        {/* Internal container centers content within the ribbon */}
        <div className="max-w-screen-lg mx-auto flex flex-col lg:flex-row items-center text-center lg:text-left px-4 sm:px-6 lg:px-8">
          <motion.div
            className="mb-8 sm:mb-10 lg:mb-0 lg:mr-16 flex-shrink-0 text-sky-600 dark:text-blue-500/90"
            variants={ribbonContentVariants}
            custom={0}
            initial="hidden"
            whileInView="visible"
            viewport={{ once: true, amount: 0.5 }}
            animate={{ rotate: 360 }}
            transition={{ duration: 25, repeat: Infinity, ease: "linear" }}
          >
            <HiOutlineCubeTransparent className="h-28 w-28 sm:h-32 sm:w-32 md:h-40 md:w-40 mx-auto lg:mx-0 opacity-90 drop-shadow-lg" />
          </motion.div>
          <motion.div
            className="flex-grow lg:max-w-3xl"
            variants={ribbonContentVariants}
            custom={1}
            initial="hidden"
            whileInView="visible"
            viewport={{ once: true, amount: 0.3 }}
          >
            <h3 className="text-2xl sm:text-3xl md:text-4xl font-bold text-slate-900 dark:text-white mb-4 sm:mb-5 tracking-tight whitespace-nowrap">
              Designed for Chemists. Fueled by Open-Source.
            </h3>
            <p className="text-base md:text-lg text-slate-700 dark:text-slate-300 leading-relaxed whitespace-nowrap">
              Built upon the open-source Cheminformatics Microservice API developed by the{" "}
              <a
                href="https://cheminf.uni-jena.de"
                target="_blank"
                rel="noopener noreferrer"
                className="text-sky-700 dark:text-sky-300 hover:text-indigo-700 dark:hover:text-cyan-300 font-medium transition-colors duration-300 hover-underline-link"
              >
                Steinbeck Lab
              </a>{" "}
              at the{" "}
              <a
                href="https://www.uni-jena.de/en"
                target="_blank"
                rel="noopener noreferrer"
                className="text-sky-700 dark:text-sky-300 hover:text-indigo-700 dark:hover:text-cyan-300 font-medium transition-colors duration-300 hover-underline-link"
              >
                Friedrich Schiller University Jena
              </a>
              .
            </p>
          </motion.div>
        </div>
      </motion.section>
      {/* Global Styles */}
      <style jsx="true" global="true">{`
        /* Ensure required CSS variables and utilities are defined */
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
        /* --- Feature Card Border Reveal --- */
        .feature-card-enhanced {
          /* ... */
        }
        .feature-card-enhanced::before {
          /* ... */
        }
        .feature-card-enhanced:hover::before {
          /* ... */
        }
        .feature-card-enhanced:hover {
          /* ... */
        }
        /* --- Hover Underline Link --- */
        .hover-underline-link {
          /* ... */
        }
        .hover-underline-link::after {
          /* ... */
        }
        .hover-underline-link:hover::after {
          /* ... */
        }
        /* --- General Styles --- */
        /* ... */
      `}</style>
    </div>
  );
};

export default HomePage;
