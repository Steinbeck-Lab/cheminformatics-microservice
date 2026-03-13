import React, { useEffect, useState } from "react";
import { motion } from "motion/react";
import { Atom, Code, Database, ExternalLink, FileText, FlaskConical, Network } from "lucide-react";
import { Button } from "@/components/ui/button";
import { Card } from "@/components/ui/card";
import { GradientMesh } from "@/components/common/GradientMesh";

// Sophisticated Animation Variants
const pageVariants = {
  hidden: { opacity: 0 },
  visible: { opacity: 1, transition: { duration: 0.7, ease: "easeOut" } },
};

const headerContainerVariants = {
  hidden: { opacity: 0 },
  visible: {
    opacity: 1,
    transition: { staggerChildren: 0.15, delayChildren: 0.2 },
  },
};

const headerItemVariants = {
  hidden: { opacity: 0, y: -20 },
  visible: {
    opacity: 1,
    y: 0,
    transition: { type: "spring", stiffness: 100, damping: 20 },
  },
};

const contentContainerVariant = {
  hidden: { opacity: 0, y: 30, scale: 0.97 },
  visible: {
    opacity: 1,
    y: 0,
    scale: 1,
    transition: { type: "spring", stiffness: 100, damping: 20, delay: 0.3 },
  },
};

const buttonVariant = {
  hidden: { opacity: 0, scale: 0.95 },
  visible: {
    opacity: 1,
    scale: 1,
    transition: { duration: 0.5, ease: [0.2, 0.65, 0.3, 0.9] },
  },
  hover: {
    scale: 1.05,
    boxShadow: "0 10px 25px -5px rgba(0, 0, 0, 0.1), 0 10px 10px -5px rgba(0, 0, 0, 0.04)",
    transition: { duration: 0.3 },
  },
};

// Enhanced CitationCard Component
const CitationCard = ({ title, citation }) => {
  const [copySuccess, setCopySuccess] = useState(false);

  const handleCopy = async () => {
    try {
      if (navigator.clipboard && window.isSecureContext) {
        await navigator.clipboard.writeText(citation);
        setCopySuccess(true);
        setTimeout(() => setCopySuccess(false), 2000);
        return;
      }

      const textArea = document.createElement("textarea");
      textArea.value = citation;
      textArea.style.position = "fixed";
      textArea.style.left = "-999999px";
      textArea.style.top = "-999999px";
      document.body.appendChild(textArea);
      textArea.focus();
      textArea.select();

      const successful = document.execCommand("copy");
      document.body.removeChild(textArea);

      if (successful) {
        setCopySuccess(true);
        setTimeout(() => setCopySuccess(false), 2000);
      }
    } catch (err) {
      console.error("Failed to copy text: ", err);
    }
  };

  return (
    <motion.div
      className="bg-white dark:bg-slate-800 p-6 rounded-2xl shadow-xs border border-slate-200 dark:border-slate-700 hover:shadow-xl hover:border-blue-300 dark:hover:border-blue-600 transition-all duration-300"
      whileHover={{
        y: -4,
      }}
      transition={{ duration: 0.2 }}
    >
      <div className="flex justify-between items-center mb-3">
        <h3 className="font-semibold text-slate-900 dark:text-white">{title}</h3>
        <Button
          onClick={handleCopy}
          className={`flex items-center justify-center p-2 rounded-lg text-white transition-all duration-300 ${
            copySuccess ? "bg-green-500 hover:bg-green-600" : "bg-blue-500 hover:bg-blue-600"
          }`}
          aria-label="Copy citation"
          title="Copy citation"
        >
          {copySuccess ? (
            <>
              <svg
                xmlns="http://www.w3.org/2000/svg"
                className="h-5 w-5 mr-1"
                viewBox="0 0 20 20"
                fill="currentColor"
              >
                <path
                  fillRule="evenodd"
                  d="M16.707 5.293a1 1 0 010 1.414l-8 8a1 1 0 01-1.414 0l-4-4a1 1 0 011.414-1.414L8 12.586l7.293-7.293a1 1 0 011.414 0z"
                  clipRule="evenodd"
                />
              </svg>
              <span className="text-xs font-medium">Copied!</span>
            </>
          ) : (
            <>
              <svg
                xmlns="http://www.w3.org/2000/svg"
                className="h-5 w-5 mr-1"
                viewBox="0 0 20 20"
                fill="currentColor"
              >
                <path d="M8 3a1 1 0 011-1h2a1 1 0 110 2H9a1 1 0 01-1-1z" />
                <path d="M6 3a2 2 0 00-2 2v11a2 2 0 002 2h8a2 2 0 002-2V5a2 2 0 00-2-2 3 3 0 01-3 3H9a3 3 0 01-3-3z" />
              </svg>
              <span className="text-xs font-medium">Copy</span>
            </>
          )}
        </Button>
      </div>
      <div className="text-sm text-slate-700 dark:text-slate-300 bg-slate-50 dark:bg-slate-900 p-4 rounded-lg border border-slate-300 dark:border-slate-700 font-mono whitespace-pre-wrap">
        {citation}
      </div>
    </motion.div>
  );
};

const AboutPage = () => {
  const [currentYear] = useState(new Date().getFullYear());
  const [isLoaded, setIsLoaded] = useState(false);

  useEffect(() => {
    const timer = setTimeout(() => {
      setIsLoaded(true);
    }, 100);
    return () => clearTimeout(timer);
  }, []);

  return (
    <motion.div
      className="flex flex-col min-h-screen bg-linear-to-br from-slate-50 to-slate-100 dark:from-slate-900 dark:to-gray-950 text-slate-900 dark:text-slate-100 font-sans overflow-x-hidden"
      style={{ maxHeight: "fit-content" }}
      variants={pageVariants}
      initial="hidden"
      animate="visible"
    >
      {/* Gradient Mesh Background (slate-blue) */}
      <GradientMesh page="about" />

      {/* Content Area */}
      <div className="relative w-full mx-auto px-4 sm:px-6 lg:px-8 py-12 md:py-16 pb-0 z-10 grow">
        <div className="w-full max-w-6xl mx-auto">
          {/* Page Header */}
          <motion.div
            className="mb-10 md:mb-16 max-w-5xl mx-auto text-center"
            variants={headerContainerVariants}
            initial="hidden"
            animate={isLoaded ? "visible" : "hidden"}
          >
            <motion.div
              variants={headerItemVariants}
              className="mx-auto w-64 sm:w-96 mb-8 relative group"
            >
              <img
                src="https://raw.githubusercontent.com/Steinbeck-Lab/cheminformatics-microservice/refs/heads/main/public/img/logo.svg"
                alt="Cheminformatics Microservice Logo"
                className="w-full filter dark:brightness-100 brightness-75 transition-all duration-500 transform hover:scale-105 drop-shadow-2xl"
              />
              <div className="absolute inset-0 bg-linear-to-r from-blue-500/0 via-blue-500/10 to-blue-500/0 opacity-0 group-hover:opacity-100 animate-pulse"></div>
            </motion.div>

            <motion.h1
              variants={headerItemVariants}
              className="text-4xl md:text-5xl lg:text-6xl font-bold mb-6"
            >
              <span className="text-gradient">Cheminformatics Microservice</span>
            </motion.h1>

            <motion.p
              variants={headerItemVariants}
              className="text-lg sm:text-xl text-slate-700 dark:text-slate-300 max-w-3xl mx-auto leading-relaxed"
            >
              A modern microservice platform providing seamless access to multiple open-source
              cheminformatics toolkits through a unified REST API interface.
            </motion.p>
          </motion.div>

          {/* Main Overview */}
          <motion.div
            className="grid grid-cols-1 md:grid-cols-2 gap-8 mb-16"
            variants={contentContainerVariant}
            initial="hidden"
            animate="visible"
          >
            <motion.div
              className="glass-bold rounded-2xl shadow-xl overflow-hidden h-full transition-all duration-500 hover:shadow-2xl hover:border-blue-200/50 dark:hover:border-blue-700/30"
              whileHover={{ y: -10 }}
              transition={{ duration: 0.3 }}
            >
              <div className="p-8 md:p-10">
                <div className="flex items-center mb-6">
                  <Network className="text-3xl text-blue-500 dark:text-blue-400 mr-4" />
                  <h2 className="text-2xl font-bold text-slate-900 dark:text-white">Overview</h2>
                </div>
                <p className="text-slate-700 dark:text-slate-300 mb-4 leading-relaxed">
                  The Cheminformatics Microservice offers a collection of versatile functions
                  accessible via REST endpoints that can handle chemical data and perform various
                  cheminformatics tasks. These tasks include (but are not limited to) generating
                  chemical structure depictions, 3D conformers, descriptors, IUPAC names, and
                  converting between machine-readable formats.
                </p>
                <p className="text-slate-700 dark:text-slate-300 leading-relaxed">
                  Researchers and developers can effectively access open-source cheminformatics
                  toolkits such as{" "}
                  <a
                    href="https://cdk.github.io/"
                    target="_blank"
                    rel="noopener noreferrer"
                    className="text-blue-600 dark:text-blue-400 hover:underline transition-colors"
                  >
                    CDK
                  </a>
                  ,{" "}
                  <a
                    href="https://www.rdkit.org/"
                    target="_blank"
                    rel="noopener noreferrer"
                    className="text-blue-600 dark:text-blue-400 hover:underline transition-colors"
                  >
                    RDKit
                  </a>
                  , and{" "}
                  <a
                    href="https://openbabel.org/"
                    target="_blank"
                    rel="noopener noreferrer"
                    className="text-blue-600 dark:text-blue-400 hover:underline transition-colors"
                  >
                    OpenBabel
                  </a>{" "}
                  through this microservice and extend them easily to suit their specific needs.
                </p>
              </div>
            </motion.div>

            <motion.div
              className="glass-bold rounded-2xl shadow-xl overflow-hidden h-full transition-all duration-500 hover:shadow-2xl hover:border-blue-200/50 dark:hover:border-blue-700/30"
              whileHover={{ y: -10 }}
              transition={{ duration: 0.3 }}
            >
              <div className="p-8 md:p-10">
                <div className="flex items-center mb-6">
                  <FlaskConical className="text-3xl text-indigo-500 dark:text-indigo-400 mr-4" />
                  <h2 className="text-2xl font-bold text-slate-900 dark:text-white">
                    Key Features
                  </h2>
                </div>
                <ul className="space-y-3 text-slate-700 dark:text-slate-300">
                  <li className="flex items-start">
                    <span className="text-blue-600 dark:text-blue-400 mr-2 text-xl">•</span>
                    <span>
                      Unified access to multiple cheminformatics toolkits (RDKit, CDK, OpenBabel)
                    </span>
                  </li>
                  <li className="flex items-start">
                    <span className="text-blue-600 dark:text-blue-400 mr-2 text-xl">•</span>
                    <span>High-quality 2D and 3D molecular visualizations</span>
                  </li>
                  <li className="flex items-start">
                    <span className="text-blue-600 dark:text-blue-400 mr-2 text-xl">•</span>
                    <span>Comprehensive molecular descriptor calculations</span>
                  </li>
                  <li className="flex items-start">
                    <span className="text-blue-600 dark:text-blue-400 mr-2 text-xl">•</span>
                    <span>Format conversions (SMILES, InChI, InChIKey, SELFIES, etc.)</span>
                  </li>
                  <li className="flex items-start">
                    <span className="text-blue-600 dark:text-blue-400 mr-2 text-xl">•</span>
                    <span>Structure standardization and validation</span>
                  </li>
                  <li className="flex items-start">
                    <span className="text-blue-600 dark:text-blue-400 mr-2 text-xl">•</span>
                    <span>Optical Chemical Structure Recognition (OCSR) using DECIMER</span>
                  </li>
                  <li className="flex items-start">
                    <span className="text-blue-600 dark:text-blue-400 mr-2 text-xl">•</span>
                    <span>Easy deployment via Docker containers</span>
                  </li>
                </ul>
              </div>
            </motion.div>
          </motion.div>

          {/* Integrated Tools Section */}
          <motion.div
            className="glass-bold rounded-2xl shadow-xl overflow-hidden mb-16 transition-all duration-500 hover:shadow-2xl hover:border-blue-200/50 dark:hover:border-blue-700/30"
            variants={contentContainerVariant}
            initial="hidden"
            animate="visible"
            whileHover={{ y: -5 }}
            transition={{ duration: 0.3 }}
          >
            <div className="p-8 md:p-10">
              <div className="flex items-center mb-8">
                <Code className="text-3xl text-purple-500 dark:text-purple-400 mr-4" />
                <h2 className="text-2xl font-bold text-slate-900 dark:text-white">
                  Integrated Tools
                </h2>
              </div>

              <div className="grid grid-cols-1 lg:grid-cols-2 gap-8 mt-6">
                <div className="bg-white/60 dark:bg-slate-800/60 p-6 rounded-xl shadow-md border border-slate-200/60 dark:border-slate-700/60">
                  <h3 className="font-semibold text-lg text-slate-900 dark:text-white mb-4">
                    Cheminformatics Toolkits
                  </h3>

                  <div className="space-y-6">
                    {/* RDKit */}
                    <motion.div
                      className="flex items-center space-x-4 group"
                      whileHover={{ x: 5 }}
                      transition={{ duration: 0.2 }}
                    >
                      <div className="shrink-0 w-16 h-16 bg-white dark:bg-slate-700 rounded-lg flex items-center justify-center p-2 shadow-md group-hover:shadow-lg transition-shadow">
                        <img
                          src="https://github.com/rdkit/rdkit/raw/71935ecae38c2f469ac0ed7d81dacd6cd6a19d5a/Docs/Images/logo.png"
                          alt="RDKit Logo"
                          className="w-full h-full object-contain"
                        />
                      </div>
                      <div>
                        <h4 className="font-semibold text-slate-900 dark:text-white">
                          <a
                            href="https://www.rdkit.org/"
                            target="_blank"
                            rel="noopener noreferrer"
                            className="hover:text-blue-600 dark:hover:text-blue-400 transition-colors"
                          >
                            RDKit
                          </a>
                        </h4>
                        <p className="text-slate-700 dark:text-slate-300 text-sm">
                          A collection of cheminformatics and machine learning tools
                        </p>
                      </div>
                    </motion.div>

                    {/* CDK */}
                    <motion.div
                      className="flex items-center space-x-4 group"
                      whileHover={{ x: 5 }}
                      transition={{ duration: 0.2 }}
                    >
                      <div className="shrink-0 w-16 h-16 bg-white dark:bg-slate-700 rounded-lg flex items-center justify-center p-2 shadow-md group-hover:shadow-lg transition-shadow">
                        <img
                          src="https://github.com/cdk/cdk.github.io/raw/source/src/img/logo.png"
                          alt="CDK Logo"
                          className="w-full h-full object-contain"
                        />
                      </div>
                      <div>
                        <h4 className="font-semibold text-slate-900 dark:text-white">
                          <a
                            href="https://cdk.github.io/"
                            target="_blank"
                            rel="noopener noreferrer"
                            className="hover:text-blue-600 dark:hover:text-blue-400 transition-colors"
                          >
                            CDK
                          </a>
                        </h4>
                        <p className="text-slate-700 dark:text-slate-300 text-sm">
                          Chemistry Development Kit, a Java library for structural cheminformatics
                        </p>
                      </div>
                    </motion.div>

                    {/* OpenBabel */}
                    <motion.div
                      className="flex items-center space-x-4 group"
                      whileHover={{ x: 5 }}
                      transition={{ duration: 0.2 }}
                    >
                      <div className="shrink-0 w-16 h-16 bg-white dark:bg-slate-700 rounded-lg flex items-center justify-center p-2 shadow-md group-hover:shadow-lg transition-shadow">
                        <img
                          src="https://upload.wikimedia.org/wikipedia/commons/c/ce/Open_Babel_computer_icon.png"
                          alt="OpenBabel Logo"
                          className="w-full h-full object-contain"
                        />
                      </div>
                      <div>
                        <h4 className="font-semibold text-slate-900 dark:text-white">
                          <a
                            href="https://openbabel.org/"
                            target="_blank"
                            rel="noopener noreferrer"
                            className="hover:text-blue-600 dark:hover:text-blue-400 transition-colors"
                          >
                            OpenBabel
                          </a>
                        </h4>
                        <p className="text-slate-700 dark:text-slate-300 text-sm">
                          A chemical toolbox designed to search, convert, analyze, or store data
                        </p>
                      </div>
                    </motion.div>
                  </div>
                </div>

                <div className="bg-white/60 dark:bg-slate-800/60 p-6 rounded-xl shadow-md border border-slate-200/60 dark:border-slate-700/60">
                  <h3 className="font-semibold text-lg text-slate-900 dark:text-white mb-4">
                    Additional Tools
                  </h3>
                  <ul className="space-y-5 text-slate-700 dark:text-slate-300">
                    <motion.li
                      className="flex items-start"
                      whileHover={{ x: 5 }}
                      transition={{ duration: 0.2 }}
                    >
                      <div className="shrink-0 w-12 h-12 bg-white dark:bg-slate-700 rounded-lg flex items-center justify-center p-2 shadow-md mr-4">
                        <img
                          src="https://github.com/Steinbeck-Lab/cheminf-jena-logos/blob/main/DECIMER/DECIMER-1.png?raw=true"
                          alt="DECIMER Logo"
                          className="w-full h-full object-contain"
                        />
                      </div>
                      <div>
                        <span className="font-semibold">DECIMER</span>
                        <p className="text-sm">Deep learning for Chemical Image Recognition</p>
                      </div>
                    </motion.li>
                    <motion.li
                      className="flex items-start"
                      whileHover={{ x: 5 }}
                      transition={{ duration: 0.2 }}
                    >
                      <div className="shrink-0 w-12 h-12 bg-white dark:bg-slate-700 rounded-lg flex items-center justify-center p-2 shadow-md mr-4">
                        <img
                          src="https://raw.githubusercontent.com/epam/ketcher/refs/heads/master/demo/public/favicon.ico"
                          alt="Ketcher Logo"
                          className="w-full h-full object-contain"
                        />
                      </div>
                      <div>
                        <span className="font-semibold">Ketcher</span>
                        <p className="text-sm">Web-based chemical structure editor</p>
                      </div>
                    </motion.li>
                    <motion.li
                      className="flex items-start"
                      whileHover={{ x: 5 }}
                      transition={{ duration: 0.2 }}
                    >
                      <div className="shrink-0 w-12 h-12 bg-white dark:bg-slate-700 rounded-lg flex items-center justify-center p-2 shadow-md mr-4">
                        <img
                          src="https://upload.wikimedia.org/wikipedia/commons/b/b6/PubChem_logo.svg"
                          alt="PubChem Logo"
                          className="w-full h-full object-contain"
                        />
                      </div>
                      <div>
                        <span className="font-semibold">PubChem PUG REST API</span>
                        <p className="text-sm">For compound search and retrieval</p>
                      </div>
                    </motion.li>
                    <motion.li
                      className="flex items-start"
                      whileHover={{ x: 5 }}
                      transition={{ duration: 0.2 }}
                    >
                      <div className="shrink-0 w-12 h-12 bg-white dark:bg-slate-700 rounded-lg flex items-center justify-center p-2 shadow-md mr-4">
                        <img
                          src="https://github.com/IUPAC-InChI/InChI-Web-Demo/blob/main/pages/img/InChI-logo2.png?raw=true"
                          alt="InChI Logo"
                          className="w-full h-full object-contain"
                        />
                      </div>
                      <div>
                        <span className="font-semibold">InChI Web Application</span>
                        <p className="text-sm">InChI evaluation from chemical depictions</p>
                      </div>
                    </motion.li>
                  </ul>
                </div>
              </div>
            </div>
          </motion.div>

          {/* Natural Products Project Button */}
          <motion.div
            className="mb-12"
            variants={headerContainerVariants}
            initial="hidden"
            animate="visible"
          >
            <h2 className="text-2xl font-bold text-slate-900 dark:text-white mb-6 text-center">
              Related Projects
            </h2>
            <div className="flex flex-wrap justify-center gap-4">
              <motion.div variants={buttonVariant} whileHover="hover" className="inline-block">
                <Button
                  asChild
                  size="lg"
                  className="h-auto inline-flex items-center justify-center px-8 py-4 bg-linear-to-r from-blue-600 to-indigo-700 hover:from-blue-500 hover:to-indigo-600 text-white font-semibold rounded-lg shadow-lg text-lg clay-interactive"
                >
                  <a href="https://naturalproducts.net" target="_blank" rel="noopener noreferrer">
                    <Atom className="mr-3 h-5 w-5" />
                    <span className="text-lg">Natural Products Online</span>
                    <ExternalLink className="ml-3 h-5 w-5" />
                  </a>
                </Button>
              </motion.div>
              <motion.div variants={buttonVariant} whileHover="hover" className="inline-block">
                <Button
                  asChild
                  size="lg"
                  className="h-auto inline-flex items-center justify-center px-8 py-4 bg-white hover:bg-slate-50 dark:bg-slate-800 dark:hover:bg-slate-700 font-semibold rounded-lg shadow-lg text-lg border border-slate-200 dark:border-slate-700 clay-interactive"
                >
                  <a
                    href="https://chemaudit.naturalproducts.net/"
                    target="_blank"
                    rel="noopener noreferrer"
                  >
                    <img
                      src="https://raw.githubusercontent.com/Kohulan/ChemAudit/main/docs/assets/logo.png"
                      alt=""
                      className="mr-3 h-6 w-auto"
                    />
                    <span className="text-lg">
                      <span className="text-orange-600">Chem</span>
                      <span className="text-slate-900 dark:text-white">Audit</span>
                    </span>
                    <ExternalLink className="ml-3 h-5 w-5 text-slate-500" />
                  </a>
                </Button>
              </motion.div>
            </div>
          </motion.div>

          {/* Public Instance */}
          <motion.div
            className="glass-bold rounded-2xl shadow-xl overflow-hidden mb-12 transition-all duration-500 hover:shadow-2xl hover:border-blue-200/50 dark:hover:border-blue-700/30"
            variants={contentContainerVariant}
            initial="hidden"
            animate="visible"
            whileHover={{ y: -5 }}
            transition={{ duration: 0.3 }}
          >
            <div className="p-8 md:p-10">
              <div className="flex items-center mb-6">
                <Database className="text-3xl text-green-500 dark:text-green-400 mr-4" />
                <h2 className="text-2xl font-bold text-slate-900 dark:text-white">
                  Public Instance
                </h2>
              </div>
              <div className="bg-blue-50/70 dark:bg-blue-900/20 rounded-xl p-6 border border-blue-200/70 dark:border-blue-800/50 mb-4">
                <p className="text-slate-700 dark:text-slate-300 mb-6">
                  A public instance of the Cheminformatics Microservice is hosted at{" "}
                  <a
                    href="https://www.uni-jena.de/en/"
                    target="_blank"
                    rel="noopener noreferrer"
                    className="text-blue-600 dark:text-blue-400 hover:underline font-medium"
                  >
                    Friedrich Schiller University Jena
                  </a>{" "}
                  in Germany, accessible at:
                </p>
                <div className="text-center">
                  <motion.div
                    whileHover={{ y: -5, scale: 1.03 }}
                    transition={{ duration: 0.2 }}
                    className="inline-block"
                  >
                    <Button
                      asChild
                      size="lg"
                      className="h-auto inline-flex items-center px-8 py-4 rounded-lg bg-linear-to-r from-blue-600 to-indigo-600 dark:from-blue-500 dark:to-indigo-500 text-white font-medium transition-all duration-300 hover:shadow-lg hover:from-blue-700 hover:to-indigo-700 clay-interactive"
                    >
                      <a
                        href="https://api.naturalproducts.net/latest/docs"
                        target="_blank"
                        rel="noopener noreferrer"
                      >
                        <ExternalLink className="mr-2 h-5 w-5" />
                        api.naturalproducts.net/latest/docs
                      </a>
                    </Button>
                  </motion.div>
                </div>
              </div>
            </div>
          </motion.div>

          {/* Citation Guidelines */}
          <motion.div
            className="glass-bold rounded-2xl shadow-xl overflow-hidden mb-12 transition-all duration-500 hover:shadow-2xl hover:border-blue-200/50 dark:hover:border-blue-700/30"
            variants={contentContainerVariant}
            initial="hidden"
            animate="visible"
            whileHover={{ y: -5 }}
            transition={{ duration: 0.3 }}
          >
            <div className="p-8 md:p-10">
              <div className="flex items-center mb-6">
                <FileText className="text-3xl text-cyan-500 dark:text-cyan-400 mr-4" />
                <h2 className="text-2xl font-bold text-slate-900 dark:text-white">
                  Citation Guidelines
                </h2>
              </div>
              <p className="text-slate-700 dark:text-slate-300 mb-6 leading-relaxed">
                It is strongly recommended that users cite both the software and the corresponding
                papers when utilizing this work. By providing appropriate citations for the
                Cheminformatics Microservice, users gain a convenient means to precisely track the
                original source of the utilized source code and data.
              </p>

              <div className="grid grid-cols-1 md:grid-cols-2 gap-8 mt-6">
                <CitationCard
                  title="Citing the paper (Application):"
                  citation="Rajan, K., Chandrasekhar, V., Sharma, N. et al. Cheminformatics Microservice V3: a web portal for chemical structure manipulation and analysis. J Cheminform 17, 142 (2025). https://doi.org/10.1186/s13321-025-01094-1"
                />
                <CitationCard
                  title="Citing the paper (Original - API):"
                  citation="Chandrasekhar, V., Sharma, N., Schaub, J. et al. Cheminformatics Microservice: unifying access to open cheminformatics toolkits. J Cheminform 15, 98 (2023). https://doi.org/10.1186/s13321-023-00762-4"
                />
                <div className="md:col-span-2 md:max-w-xl md:mx-auto md:w-full">
                  <CitationCard
                    title="Citing the software:"
                    citation="Venkata, C., Sharma, N., Schaub, J., Steinbeck, C., & Rajan, K. (2025). cheminformatics-microservice (Version v3.5.0) [Computer software]. https://doi.org/10.5281/zenodo.16890410"
                  />
                </div>
              </div>
            </div>
          </motion.div>

          {/* License Information */}
          <motion.div
            className="glass-bold rounded-2xl shadow-xl overflow-hidden mb-12 transition-all duration-500 hover:shadow-2xl hover:border-blue-200/50 dark:hover:border-blue-700/30"
            variants={contentContainerVariant}
            initial="hidden"
            animate="visible"
            whileHover={{ y: -5 }}
            transition={{ duration: 0.3 }}
          >
            <div className="p-8 md:p-10">
              <div className="flex items-center mb-6">
                <svg
                  className="w-8 h-8 text-blue-500 dark:text-blue-400 mr-4"
                  fill="currentColor"
                  viewBox="0 0 24 24"
                  xmlns="http://www.w3.org/2000/svg"
                >
                  <path d="M12 22C6.477 22 2 17.523 2 12S6.477 2 12 2s10 4.477 10 10-4.477 10-10 10zm0-2a8 8 0 1 0 0-16 8 8 0 0 0 0 16zM11 7h2v8h-2V7zm0 10h2v2h-2v-2z"></path>
                </svg>
                <h2 className="text-2xl font-bold text-slate-900 dark:text-white">
                  License Information
                </h2>
              </div>
              <div className="grid grid-cols-1 md:grid-cols-2 gap-6">
                <motion.div
                  className="bg-white/60 dark:bg-slate-800/60 p-6 rounded-xl shadow-md border border-slate-200/60 dark:border-slate-700/60"
                  whileHover={{
                    y: -5,
                    boxShadow:
                      "0 20px 25px -5px rgba(0, 0, 0, 0.1), 0 10px 10px -5px rgba(0, 0, 0, 0.04)",
                  }}
                  transition={{ duration: 0.2 }}
                >
                  <h3 className="font-semibold text-slate-900 dark:text-white mb-3">
                    Cheminformatics Microservice
                  </h3>
                  <p className="text-slate-700 dark:text-slate-300 mb-2">
                    Released under the{" "}
                    <a
                      href="https://opensource.org/licenses/MIT"
                      target="_blank"
                      rel="noopener noreferrer"
                      className="text-blue-600 dark:text-blue-400 hover:underline"
                    >
                      MIT license
                    </a>
                  </p>
                  <p className="text-slate-700 dark:text-slate-300">Copyright © {currentYear}</p>
                </motion.div>
                <motion.div
                  className="bg-white/60 dark:bg-slate-800/60 p-6 rounded-xl shadow-md border border-slate-200/60 dark:border-slate-700/60"
                  whileHover={{
                    y: -5,
                    boxShadow:
                      "0 20px 25px -5px rgba(0, 0, 0, 0.1), 0 10px 10px -5px rgba(0, 0, 0, 0.04)",
                  }}
                  transition={{ duration: 0.2 }}
                >
                  <h3 className="font-semibold text-slate-900 dark:text-white mb-3">
                    Integrated Components
                  </h3>
                  <p className="text-slate-700 dark:text-slate-300 mb-2">
                    <strong>Ketcher:</strong> Apache 2.0 License
                    <br />
                    Copyright © {currentYear} EPAM Systems, Inc.
                  </p>
                  <p className="text-slate-700 dark:text-slate-300">
                    Other components are used according to their respective licenses.
                  </p>
                </motion.div>
              </div>
            </div>
          </motion.div>

          {/* Acknowledgements */}
          <motion.div
            className="glass-bold rounded-2xl shadow-xl overflow-hidden mb-8"
            variants={contentContainerVariant}
            initial="hidden"
            animate="visible"
            whileHover={{ y: -5 }}
            transition={{ duration: 0.3 }}
          >
            <div className="p-8 md:p-10">
              <div className="flex items-center mb-6">
                <svg
                  className="w-8 h-8 text-purple-500 dark:text-purple-400 mr-4"
                  fill="currentColor"
                  viewBox="0 0 24 24"
                  xmlns="http://www.w3.org/2000/svg"
                >
                  <path d="M12 18.26l-7.053 3.948 1.575-7.928L.587 8.792l8.027-.952L12 .5l3.386 7.34 8.027.952-5.935 5.488 1.575 7.928L12 18.26zm0-2.292l4.247 2.377-.949-4.773 3.573-3.305-4.833-.573L12 5.275l-2.038 4.42-4.833.572 3.573 3.305-.949 4.773L12 15.968z"></path>
                </svg>
                <h2 className="text-2xl font-bold text-slate-900 dark:text-white">
                  Acknowledgments and Maintenance
                </h2>
              </div>
              <p className="text-slate-700 dark:text-slate-300 mb-6 leading-relaxed">
                Cheminformatics Microservice,{" "}
                <a
                  href="https://naturalproducts.net"
                  target="_blank"
                  rel="noopener noreferrer"
                  className="text-blue-600 dark:text-blue-400 hover:underline"
                >
                  Natural Products Online
                </a>
                , and{" "}
                <a
                  href="https://chemaudit.naturalproducts.net/"
                  target="_blank"
                  rel="noopener noreferrer"
                  className="text-blue-600 dark:text-blue-400 hover:underline"
                >
                  ChemAudit
                </a>{" "}
                are developed and maintained by the{" "}
                <a
                  href="https://cheminf.uni-jena.de"
                  target="_blank"
                  rel="noopener noreferrer"
                  className="text-blue-600 dark:text-blue-400 hover:underline"
                >
                  Steinbeck group
                </a>{" "}
                at the{" "}
                <a
                  href="https://www.uni-jena.de/en/"
                  target="_blank"
                  rel="noopener noreferrer"
                  className="text-blue-600 dark:text-blue-400 hover:underline"
                >
                  Friedrich Schiller University
                </a>{" "}
                Jena, Germany.
              </p>

              <div className="border-t border-slate-200/80 dark:border-slate-700/80 pt-6">
                <p className="text-slate-700 dark:text-slate-300 mb-8 leading-relaxed">
                  Funded by the{" "}
                  <a
                    href="https://www.dfg.de/"
                    target="_blank"
                    rel="noopener noreferrer"
                    className="text-blue-600 dark:text-blue-400 hover:underline"
                  >
                    DFG
                  </a>{" "}
                  (German Research Foundation) under{" "}
                  <a
                    href="https://nfdi4chem.de/"
                    target="_blank"
                    rel="noopener noreferrer"
                    className="text-blue-600 dark:text-blue-400 hover:underline"
                  >
                    NFDI4Chem
                  </a>{" "}
                  (Project: <strong>441958208</strong>) and{" "}
                  <a
                    href="https://www.chembiosys.de/en/"
                    target="_blank"
                    rel="noopener noreferrer"
                    className="text-blue-600 dark:text-blue-400 hover:underline"
                  >
                    ChemBioSys
                  </a>{" "}
                  (Project INF, SFB 1127, Project: <strong>239748522</strong>).
                </p>

                {/* Funding Logos */}
                <div className="flex flex-wrap justify-center items-center gap-6 md:gap-8">
                  {[
                    {
                      name: "DFG",
                      url: "https://www.dfg.de/",
                      logo: "https://github.com/Steinbeck-Lab/cheminformatics-microservice/blob/main/docs/public/dfg_logo_schriftzug_blau_foerderung_en.gif?raw=true",
                    },
                    {
                      name: "NFDI4Chem",
                      url: "https://nfdi4chem.de/",
                      logo: "https://www.nfdi4chem.de/wp-content/themes/wptheme/assets/img/logo.svg",
                    },
                    {
                      name: "Chembiosys",
                      url: "https://www.chembiosys.de/en/welcome.html",
                      logo: "https://github.com/Steinbeck-Lab/cheminformatics-microservice/assets/30716951/45c8e153-8322-4563-a51d-cbdbe4e08627",
                    },
                  ].map((funder) => (
                    <motion.div
                      key={funder.name}
                      whileHover={{ y: -5, scale: 1.05 }}
                      transition={{ duration: 0.2 }}
                    >
                      <Card className="glass-bold p-3 rounded-lg shadow-md hover:shadow-lg transition-all duration-300 py-0 gap-0">
                        <a
                          href={funder.url}
                          target="_blank"
                          rel="noopener noreferrer"
                          className="w-auto h-16 flex items-center justify-center"
                        >
                          <img
                            src={funder.logo}
                            alt={`${funder.name} Logo`}
                            className="max-h-12 object-contain"
                          />
                        </a>
                      </Card>
                    </motion.div>
                  ))}
                </div>
              </div>
            </div>
          </motion.div>
        </div>
      </div>

      {/* Fix for potential extra space */}
      <div style={{ height: 0, clear: "both", overflow: "hidden" }} />

      {/* GradientMesh + tailwind.css utilities provide glass/text-gradient -- no inline styles needed */}
    </motion.div>
  );
};
export default AboutPage;
