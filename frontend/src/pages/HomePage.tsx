/**
 * HomePage — Immersive hero with 3D caffeine molecule background,
 * glassmorphism card, and claymorphism buttons.
 *
 * Wide-screen optimised: hero fills the viewport, content uses max-w-7xl.
 * 3D molecule renders on desktop only; mobile gets the GradientMesh fallback.
 */
import React, { useState, lazy, Suspense } from "react";
import { Link } from "react-router-dom";
import { motion } from "motion/react";
import { useAppContext } from "../context/AppContext";
import MoleculeCard from "../components/common/MoleculeCard";
import {
  ArrowLeftRight,
  ArrowRight,
  BookOpen,
  Camera,
  Code,
  Eye,
  FlaskConical,
  Wrench,
} from "lucide-react";
import { Button } from "@/components/ui/button";
import { Card } from "@/components/ui/card";
import { GradientMesh } from "@/components/common/GradientMesh";

// Lazy-load 3D scene (heavy — Three.js bundle)
const CaffeineMolecule3D = lazy(() => import("@/components/3d/CaffeineMolecule3D"));

// --- Feature Configuration ---
const features = [
  {
    title: "Chemical Analysis",
    description:
      "Analyze molecules with descriptors, stereoisomer generation, and similarity calculations.",
    icon: FlaskConical,
    link: "/chem",
  },
  {
    title: "Format Conversion",
    description:
      "Convert between different chemical file formats (SMILES, InChI, 2D/3D coordinates).",
    icon: ArrowLeftRight,
    link: "/convert",
  },
  {
    title: "Visualization",
    description: "Generate customizable 2D and 3D visualizations of chemical structures.",
    icon: Eye,
    link: "/depict",
  },
  {
    title: "Structure Tools",
    description: "Specialized tools for structure generation, sugar removal, and more.",
    icon: Wrench,
    link: "/tools",
  },
  {
    title: "OCSR",
    description: "Optical Chemical Structure Recognition to extract structures from images.",
    icon: Camera,
    link: "/ocsr",
  },
];

// --- Animation Variants ---
const heroTextVariants = {
  hidden: { opacity: 0, y: -20 },
  visible: (i = 1) => ({
    opacity: 1,
    y: 0,
    transition: { type: "spring", stiffness: 100, damping: 20, delay: i * 0.07 },
  }),
};

const containerVariants = {
  hidden: { opacity: 0 },
  visible: {
    opacity: 1,
    transition: { staggerChildren: 0.08, delayChildren: 0.15 },
  },
};

const itemVariants = {
  hidden: { y: 25, opacity: 0, scale: 0.97 },
  visible: {
    y: 0,
    opacity: 1,
    scale: 1,
    transition: { type: "spring", stiffness: 100, damping: 20 },
  },
};

const sectionFadeIn = {
  hidden: { opacity: 0, y: 20 },
  visible: { opacity: 1, y: 0, transition: { type: "spring", stiffness: 100, damping: 20 } },
};

const ribbonContentVariants = {
  hidden: { opacity: 0, y: 20 },
  visible: (i = 1) => ({
    opacity: 1,
    y: 0,
    transition: { duration: 0.7, delay: i * 0.1, ease: [0.2, 0.65, 0.3, 0.9] },
  }),
};

// --- Tilt Card (for feature cards) ---
const TiltCard = ({
  children,
  className = "",
  tiltIntensity = 7,
}: {
  children: React.ReactNode;
  className?: string;
  tiltIntensity?: number;
}) => {
  const [rotate, setRotate] = useState({ x: 0, y: 0 });

  const onMouseMove = (e: React.MouseEvent) => {
    if ("ontouchstart" in window && e.nativeEvent instanceof TouchEvent) return;
    const card = e.currentTarget;
    const { width, height, left, top } = card.getBoundingClientRect();
    const x = e.clientX - left;
    const y = e.clientY - top;
    const rotateX = (y / height - 0.5) * -tiltIntensity;
    const rotateY = (x / width - 0.5) * tiltIntensity;
    setRotate({ x: rotateX, y: rotateY });
  };

  const onMouseLeave = () => setRotate({ x: 0, y: 0 });

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
    <div className="relative min-h-screen w-full text-slate-900 dark:text-slate-100 font-sans overflow-x-hidden isolate">
      {/* Gradient mesh — visible on mobile, sits behind 3D on desktop */}
      <GradientMesh page="home" />

      {/* 3D Caffeine molecule background — desktop only, lazy-loaded */}
      <Suspense fallback={null}>
        <CaffeineMolecule3D />
      </Suspense>

      {/* ============================================================ */}
      {/* HERO SECTION — Full viewport, glass card centered             */}
      {/* ============================================================ */}
      <div className="relative z-10 flex items-center justify-center min-h-[92vh] w-full px-4 sm:px-6 lg:px-8">
        <motion.section
          className="w-full max-w-6xl mx-auto text-center"
          initial="hidden"
          animate="visible"
        >
          {/* Glass hero card */}
          <div className="glass-bold rounded-3xl p-8 sm:p-12 md:p-16 lg:p-20 shadow-2xl">
            {/* Headline */}
            <motion.h1
              variants={heroTextVariants}
              custom={0}
              className="text-4xl sm:text-5xl md:text-6xl lg:text-7xl xl:text-8xl font-extrabold text-slate-900 dark:text-white mb-5 sm:mb-6 leading-[1.08] tracking-tighter"
            >
              Seamless Access to{" "}
              <motion.span
                variants={heroTextVariants}
                custom={1}
                className="text-transparent bg-clip-text bg-linear-to-r from-sky-600 via-cyan-500 to-indigo-600 dark:from-sky-400 dark:via-cyan-400 dark:to-emerald-400"
              >
                Open Cheminformatics Tools
              </motion.span>
            </motion.h1>

            {/* Sub-headline */}
            <motion.p
              variants={heroTextVariants}
              custom={2}
              className="text-base sm:text-lg md:text-xl lg:text-2xl text-slate-600 dark:text-slate-300 max-w-4xl mx-auto mb-10 sm:mb-12 leading-relaxed"
            >
              Visualize, analyze, and manipulate chemical structures seamlessly using the powerful
              Cheminformatics Microservice API.
            </motion.p>

            {/* CTA Buttons — claymorphism */}
            <motion.div
              variants={heroTextVariants}
              custom={3}
              className="flex flex-col sm:flex-row flex-wrap justify-center items-center gap-4 md:gap-6"
            >
              {/* Primary CTA */}
              <Button
                asChild
                size="lg"
                className="clay-interactive group w-full sm:w-auto relative h-auto inline-flex items-center justify-center px-8 py-3.5 sm:px-10 sm:py-4 text-base sm:text-lg font-semibold text-white bg-linear-to-r from-sky-600 to-cyan-500 dark:from-blue-600 dark:to-cyan-500 rounded-2xl shadow-lg hover:shadow-xl hover:brightness-110 dark:hover:shadow-cyan-500/30 transition-all duration-300 ease-out transform hover:scale-[1.04] active:scale-[0.98] focus:outline-hidden focus:ring-4 ring-offset-2 ring-offset-slate-100 dark:ring-offset-gray-950 focus:ring-cyan-500/50 cursor-pointer"
              >
                <Link to="/depict">
                  <span className="relative z-10 flex items-center">
                    Get Started{" "}
                    <ArrowRight className="ml-2.5 h-5 w-5 transition-transform duration-300 group-hover:translate-x-1.5" />
                  </span>
                </Link>
              </Button>

              {/* Secondary: Guides */}
              <Button
                asChild
                variant="outline"
                size="lg"
                className="clay-interactive group w-full sm:w-auto relative h-auto inline-flex items-center justify-center px-8 py-3.5 sm:px-10 sm:py-4 text-base sm:text-lg font-semibold text-slate-700 dark:text-slate-200 glass-bold rounded-2xl border-white/30 dark:border-slate-600/30 hover:bg-white/20 dark:hover:bg-slate-700/30 shadow-lg hover:shadow-xl transition-all duration-300 ease-out transform hover:scale-[1.04] active:scale-[0.98] focus:outline-hidden focus:ring-4 focus:ring-slate-400/50 dark:focus:ring-slate-600/50 ring-offset-2 ring-offset-slate-100 dark:ring-offset-gray-950 cursor-pointer"
              >
                <a
                  href="https://docs.api.naturalproducts.net/introduction.html"
                  target="_blank"
                  rel="noopener noreferrer"
                >
                  <span className="relative z-10 flex items-center">
                    <Code className="mr-2 h-5 w-5 text-green-600 dark:text-green-400" /> Guides
                  </span>
                </a>
              </Button>

              {/* Secondary: API Docs */}
              <Button
                asChild
                variant="outline"
                size="lg"
                className="clay-interactive group w-full sm:w-auto relative h-auto inline-flex items-center justify-center px-8 py-3.5 sm:px-10 sm:py-4 text-base sm:text-lg font-semibold text-slate-700 dark:text-slate-200 glass-bold rounded-2xl border-white/30 dark:border-slate-600/30 hover:bg-white/20 dark:hover:bg-slate-700/30 shadow-lg hover:shadow-xl transition-all duration-300 ease-out transform hover:scale-[1.04] active:scale-[0.98] focus:outline-hidden focus:ring-4 focus:ring-slate-400/50 dark:focus:ring-slate-600/50 ring-offset-2 ring-offset-slate-100 dark:ring-offset-gray-950 cursor-pointer"
              >
                <a
                  href="https://api.naturalproducts.net/latest/docs"
                  target="_blank"
                  rel="noopener noreferrer"
                >
                  <span className="relative z-10 flex items-center">
                    <BookOpen className="mr-2 h-5 w-5 text-blue-600 dark:text-blue-400" /> API Docs
                  </span>
                </a>
              </Button>
            </motion.div>
          </div>
        </motion.section>
      </div>

      {/* ============================================================ */}
      {/* BELOW THE FOLD — Feature cards, recent molecules, ribbon     */}
      {/* ============================================================ */}
      <div className="relative w-full mx-auto px-4 sm:px-6 lg:px-8 pb-12 md:pb-16 z-10">
        <div className="max-w-7xl mx-auto">
          {/* --- Feature Cards --- */}
          <motion.section
            className="mb-16 md:mb-24"
            variants={containerVariants}
            initial="hidden"
            whileInView="visible"
            viewport={{ once: true, amount: 0.1 }}
          >
            <h2 className="text-3xl sm:text-4xl md:text-5xl font-bold text-center mb-12 md:mb-16 tracking-tight text-slate-900 dark:text-white">
              Core
              <span className="text-transparent bg-clip-text bg-linear-to-r from-sky-600 to-indigo-600 dark:from-blue-400 dark:to-cyan-400">
                {" "}
                Capabilities
              </span>
            </h2>
            <div className="grid grid-cols-1 sm:grid-cols-2 lg:grid-cols-3 gap-6">
              {features.map((feature, index) => {
                const Icon = feature.icon;
                return (
                  <motion.div key={index} variants={itemVariants}>
                    <TiltCard tiltIntensity={5}>
                      <Link to={feature.link} className="block h-full">
                        <Card className="glass-bold clay-interactive rounded-2xl p-6 sm:p-8 h-full flex flex-col hover:-translate-y-1 transition-transform duration-200 border-0 gap-0 cursor-pointer">
                          <div className="mb-5 flex items-center justify-center h-14 w-14 rounded-2xl bg-linear-to-br from-slate-100 to-sky-100 dark:from-blue-700/60 dark:to-cyan-700/60 shadow-md">
                            <Icon className="h-7 w-7 text-sky-700 dark:text-sky-200" />
                          </div>
                          <h3 className="text-lg sm:text-xl font-semibold text-slate-800 dark:text-white mb-2">
                            {feature.title}
                          </h3>
                          <p className="text-sm text-slate-600 dark:text-slate-300 grow mb-4 leading-relaxed">
                            {feature.description}
                          </p>
                          <div className="mt-auto flex items-center text-sky-700 dark:text-sky-400 font-medium text-sm group">
                            <span>Learn more</span>
                            <ArrowRight className="ml-2 h-4 w-4 transition-transform duration-300 group-hover:translate-x-1.5" />
                          </div>
                        </Card>
                      </Link>
                    </TiltCard>
                  </motion.div>
                );
              })}
            </div>
          </motion.section>

          {/* --- Recent Molecules --- */}
          {recentMolecules && recentMolecules.length > 0 && (
            <motion.section
              className="mb-16 md:mb-24"
              variants={sectionFadeIn}
              initial="hidden"
              whileInView="visible"
              viewport={{ once: true, amount: 0.15 }}
            >
              <div className="glass-bold rounded-3xl p-6 sm:p-8 md:p-10">
                <div className="flex flex-col sm:flex-row justify-between items-start sm:items-center mb-8 gap-3 sm:gap-6">
                  <h2 className="text-2xl sm:text-3xl md:text-4xl font-bold text-slate-900 dark:text-white tracking-tight">
                    Recently Viewed
                  </h2>
                  <Link
                    to="/chem"
                    className="text-sky-700 dark:text-sky-400 hover:text-indigo-700 dark:hover:text-cyan-200 shrink-0 flex items-center font-medium transition-colors duration-300 group text-sm"
                  >
                    <span>View all</span>
                    <ArrowRight className="ml-1.5 h-4 w-4 transition-transform duration-300 group-hover:translate-x-1" />
                  </Link>
                </div>
                <motion.div
                  className="grid grid-cols-1 sm:grid-cols-2 md:grid-cols-3 gap-4"
                  variants={containerVariants}
                  initial="hidden"
                  whileInView="visible"
                  viewport={{ once: true, amount: 0.2 }}
                >
                  {recentMolecules.slice(0, 3).map((molecule, index) => (
                    <motion.div key={index} variants={itemVariants}>
                      <Card className="relative bg-white/80 dark:bg-slate-800/70 backdrop-blur-lg rounded-2xl shadow-lg dark:shadow-xl border-slate-200 dark:border-slate-700/80 overflow-hidden transition-all duration-300 hover:border-slate-300 dark:hover:border-slate-600 p-1 py-0 gap-0">
                        <MoleculeCard
                          smiles={molecule.smiles}
                          title={molecule.name || `Molecule ${index + 1}`}
                          size="sm"
                          className="p-4 text-slate-800 dark:text-white bg-transparent"
                        />
                      </Card>
                    </motion.div>
                  ))}
                </motion.div>
              </div>
            </motion.section>
          )}
        </div>
      </div>

      {/* --- API Ribbon Section --- */}
      <motion.section
        className="relative bg-linear-to-r from-sky-100 via-indigo-100 to-fuchsia-100 dark:from-gray-900 dark:via-indigo-950/80 dark:to-fuchsia-950/70 py-12 sm:py-16 md:py-20"
        variants={sectionFadeIn}
        initial="hidden"
        whileInView="visible"
        viewport={{ once: true, amount: 0.2 }}
      >
        <div className="max-w-5xl mx-auto flex flex-col lg:flex-row items-center text-center lg:text-left px-4 sm:px-6 lg:px-8">
          <motion.div
            className="mb-8 sm:mb-10 lg:mb-0 lg:mr-10 shrink-0 text-sky-600 dark:text-blue-500/90"
            variants={ribbonContentVariants}
            custom={0}
            initial="hidden"
            whileInView="visible"
            viewport={{ once: true, amount: 0.5 }}
            animate={{ rotate: 360 }}
            transition={{ duration: 25, repeat: Infinity, ease: "linear" }}
          >
            <svg
              className="h-28 w-28 sm:h-32 sm:w-32 md:h-36 md:w-36 mx-auto lg:mx-0 opacity-90 drop-shadow-lg"
              viewBox="0 0 24 24"
              fill="none"
              stroke="currentColor"
              strokeWidth="1.5"
              strokeLinecap="round"
              strokeLinejoin="round"
            >
              <path d="M14 10l-2 1m0 0l-2-1m2 1v2.5M20 7l-2 1m2-1l-2-1m2 1v2.5M14 4l-2-1-2 1M4 7l2-1M4 7l2 1M4 7v2.5M12 21l-2-1m2 1l2-1m-2 1v-2.5M6 18l-2-1v-2.5M18 18l2-1v-2.5" />
            </svg>
          </motion.div>
          <motion.div
            className="grow min-w-0"
            variants={ribbonContentVariants}
            custom={1}
            initial="hidden"
            whileInView="visible"
            viewport={{ once: true, amount: 0.3 }}
          >
            <h3 className="text-2xl sm:text-3xl md:text-4xl font-bold text-slate-900 dark:text-white mb-4 sm:mb-5 tracking-tight whitespace-nowrap">
              Designed for Chemists. Fueled by Open-Source.
            </h3>
            <p className="text-base md:text-lg text-slate-700 dark:text-slate-300 leading-relaxed">
              Built upon the open-source Cheminformatics Microservice API developed by the{" "}
              <a
                href="https://cheminf.uni-jena.de"
                target="_blank"
                rel="noopener noreferrer"
                className="text-sky-700 dark:text-sky-300 hover:text-indigo-700 dark:hover:text-cyan-300 font-medium transition-colors duration-300 whitespace-nowrap"
              >
                Steinbeck Lab
              </a>{" "}
              at the{" "}
              <a
                href="https://www.uni-jena.de/en"
                target="_blank"
                rel="noopener noreferrer"
                className="text-sky-700 dark:text-sky-300 hover:text-indigo-700 dark:hover:text-cyan-300 font-medium transition-colors duration-300 whitespace-nowrap"
              >
                Friedrich Schiller University Jena
              </a>
              .
            </p>
          </motion.div>
        </div>
      </motion.section>
    </div>
  );
};

export default HomePage;
