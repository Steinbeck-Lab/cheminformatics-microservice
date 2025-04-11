import React, { useEffect, useState } from "react";
import { motion, useScroll, useTransform } from "framer-motion";
import {
  HiOutlineExternalLink,
  HiOutlineCode,
  HiOutlineDocumentText,
  HiOutlineDatabase,
} from "react-icons/hi";
import { FaFlask, FaConnectdevelop, FaAtom } from "react-icons/fa";

// Animation Variants
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
    transition: { duration: 0.7, ease: [0.2, 0.65, 0.3, 0.9] },
  },
};

const contentContainerVariant = {
  hidden: { opacity: 0, y: 30, scale: 0.97 },
  visible: {
    opacity: 1,
    y: 0,
    scale: 1,
    transition: { duration: 0.7, delay: 0.3, ease: [0.2, 0.65, 0.3, 0.9] },
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
    boxShadow:
      "0 10px 25px -5px rgba(0, 0, 0, 0.1), 0 10px 10px -5px rgba(0, 0, 0, 0.04)",
    transition: { duration: 0.3 },
  },
};

// CitationCard component with copy functionality
const CitationCard = ({ title, citation }) => {
  const [copySuccess, setCopySuccess] = useState(false);

  const handleCopy = async () => {
    try {
      // Try to use the Clipboard API first (modern browsers)
      if (navigator.clipboard && window.isSecureContext) {
        await navigator.clipboard.writeText(citation);
        setCopySuccess(true);
        setTimeout(() => setCopySuccess(false), 2000);
        return;
      }

      // Fallback method for browsers that don't support clipboard API
      // or when not in a secure context (HTTPS)
      const textArea = document.createElement("textarea");
      textArea.value = citation;

      // Make the textarea out of viewport
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
      className="bg-white/60 dark:bg-slate-800/60 p-6 rounded-xl shadow-md border border-slate-200/60 dark:border-slate-700/60"
      whileHover={{
        y: -5,
        boxShadow:
          "0 20px 25px -5px rgba(0, 0, 0, 0.1), 0 10px 10px -5px rgba(0, 0, 0, 0.04)",
      }}
      transition={{ duration: 0.2 }}
    >
      <div className="flex justify-between items-center mb-3">
        <h3 className="font-semibold text-slate-900 dark:text-white">
          {title}
        </h3>
        <button
          onClick={handleCopy}
          className={`flex items-center justify-center p-2 rounded-lg text-white transition-all duration-300 ${
            copySuccess
              ? "bg-green-500 hover:bg-green-600"
              : "bg-blue-500 hover:bg-blue-600"
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
        </button>
      </div>
      <div className="text-sm text-slate-700 dark:text-slate-300 bg-white/80 dark:bg-slate-900/80 p-4 rounded-lg border border-slate-300/80 dark:border-slate-700/80 font-mono whitespace-pre-wrap">
        {citation}
      </div>
    </motion.div>
  );
};

const AboutPage = () => {
  const [currentYear] = useState(new Date().getFullYear());
  const { scrollYProgress } = useScroll();
  const backgroundY = useTransform(scrollYProgress, [0, 1], ["0%", "20%"]);
  const backgroundScale = useTransform(scrollYProgress, [0, 1], [1, 1.1]);
  const rotation = useTransform(scrollYProgress, [0, 0.5], [0, 5]);
  const [isLoaded, setIsLoaded] = useState(false);

  useEffect(() => {
    // Add a slight delay to ensure smooth animations after component mounts
    const timer = setTimeout(() => {
      setIsLoaded(true);
    }, 100);
    return () => clearTimeout(timer);
  }, []);

  return (
    <motion.div
      className="flex flex-col min-h-screen bg-gradient-to-br from-slate-50 to-slate-100 dark:from-slate-900 dark:to-gray-950 text-slate-900 dark:text-slate-100 font-sans overflow-x-hidden"
      style={{ maxHeight: "fit-content" }}
      variants={pageVariants}
      initial="hidden"
      animate="visible"
    >
      {/* Enhanced Background Effects */}
      <motion.div
        className="absolute inset-0 -z-20 overflow-hidden pointer-events-none"
        style={{ y: backgroundY, scale: backgroundScale, maxHeight: "100vh" }}
      >
        <div className="absolute inset-0 dark:bg-gradient-to-br dark:from-slate-900 dark:via-indigo-950/50 dark:to-slate-950 opacity-0 dark:opacity-100 transition-opacity duration-500"></div>
        <div className="absolute inset-0 bg-gradient-to-br from-blue-50/80 via-white to-indigo-50/80 opacity-100 dark:opacity-0 transition-opacity duration-500"></div>
        <div className="absolute inset-0 w-full h-full bg-[radial-gradient(#e5e7eb_1px,transparent_1px)] dark:bg-[radial-gradient(#1e293b30_1px,transparent_1px)] [background-size:20px_20px] opacity-50 dark:opacity-30"></div>
      </motion.div>

      {/* Enhanced Accent Patterns */}
      <div
        className="absolute -z-10 inset-0 overflow-hidden opacity-30 dark:opacity-20 pointer-events-none"
        style={{ maxHeight: "100vh" }}
      >
        <motion.svg
          className="absolute -top-[40%] -right-[30%] w-[80%] h-auto text-blue-500/10 dark:text-blue-300/10"
          viewBox="0 0 200 200"
          xmlns="http://www.w3.org/2000/svg"
          style={{ rotate: rotation }}
        >
          <path
            fill="currentColor"
            d="M39.5,-65.3C50.3,-56.7,57.6,-44.3,64.2,-31.1C70.8,-17.9,76.7,-3.8,74.6,9.4C72.6,22.5,62.6,34.8,51.5,43.9C40.5,53.1,28.3,59.1,14.9,63.5C1.5,68,-13.2,70.8,-24.5,65.8C-35.9,60.7,-43.9,47.7,-52.3,35.2C-60.6,22.7,-69.3,10.8,-70.8,-2.2C-72.3,-15.2,-66.5,-29.3,-57.1,-39.5C-47.7,-49.7,-34.5,-56,-21.5,-62.3C-8.6,-68.6,4.2,-74.9,17.4,-74.2C30.7,-73.5,44.5,-65.8,55.6,-55.2C66.7,-44.7,77.2,-31.3,77.1,-18.2C77.1,-5,66.7,7.9,59.1,19.5C51.6,31,47,41.1,38.7,47.8C30.5,54.5,18.6,57.7,6.5,59.1C-5.7,60.5,-18,60.1,-29.9,56.3C-41.8,52.5,-53.3,45.4,-59,35.1C-64.8,24.8,-64.8,11.3,-64.7,-2.2C-64.5,-15.7,-64.2,-29.4,-58.5,-40.1C-52.8,-50.9,-41.9,-58.8,-30,-64.3C-18.1,-69.9,-5.3,-73.1,7.7,-71.9C20.7,-70.8,34.6,-65.2,45.4,-56.4C56.1,-47.6,63.8,-35.5,67.4,-22.4C71,-9.3,70.5,4,64.9,14.2C59.3,24.5,48.5,31.9,39,40.3C29.4,48.8,21.1,58.4,10,63.9C-1.1,69.4,-15,71,-26.4,66.9C-37.9,62.9,-46.9,53.3,-54.6,42.8C-62.3,32.3,-68.7,20.8,-70.4,8.6C-72.1,-3.7,-69.2,-16.8,-63.9,-29C-58.6,-41.1,-51,-52.4,-40.7,-60.3C-30.4,-68.2,-17.5,-72.9,-3.5,-74.3C10.5,-75.7,22.6,-73.8,33.3,-68.9C44,-64,53.2,-56,58.5,-45.6C63.7,-35.2,65.1,-22.3,64.7,-9.9C64.3,2.4,62.1,14.3,57.1,24.7C52.1,35.1,44.4,44.2,34.9,50.6C25.4,57,14.1,60.7,2.4,62.9C-9.4,65.1,-21.5,65.7,-33.5,62.3C-45.5,58.9,-57.4,51.3,-63.3,40.5C-69.2,29.6,-69.2,15.6,-71.2,1.1C-73.3,-13.4,-77.6,-28.3,-73.1,-40.3C-68.6,-52.3,-55.4,-61.3,-41.6,-67.5C-27.9,-73.7,-13.5,-77.1,0.2,-77.2C13.9,-77.3,27.8,-74.2,39.5,-67.3C51.2,-60.4,60.7,-49.7,66.4,-37.5C72.1,-25.3,74,-11.6,72.5,-0.1C71,11.4,66.1,20.7,60.3,29.7C54.5,38.7,47.8,47.4,38.9,53C30,58.6,18.9,61.1,6.9,65.1C-5.1,69.2,-17.9,74.9,-29.3,73.5C-40.7,72.1,-50.6,63.5,-58.5,53.4C-66.3,43.3,-72,31.5,-74.8,18.9C-77.6,6.2,-77.4,-7.4,-73.7,-19.8C-70,-32.2,-62.9,-43.5,-52.9,-52.2C-42.9,-60.9,-30.1,-67,-17.5,-71.1C-4.9,-75.2,7.5,-77.3,19.7,-76.6C31.9,-75.9,43.9,-72.5,52.9,-65.3C61.9,-58.2,67.9,-47.3,72.2,-35.6C76.5,-23.9,79,-11.3,76.8,0.5C74.6,12.2,67.7,23.4,60.3,33.9C52.9,44.3,45,54,34.8,61.2C24.6,68.4,12.3,73.1,-0.6,74.6C-13.5,76.1,-27,74.4,-39.7,69.9C-52.4,65.4,-64.4,58,-68.6,46.9C-72.9,35.8,-69.5,20.9,-69.3,7.8C-69,-5.4,-71.8,-17.8,-68.4,-28.7C-65,-39.5,-55.3,-48.8,-44.4,-56.6C-33.4,-64.4,-21.2,-70.7,-8,-72.2C5.2,-73.7,19.4,-70.3,31.1,-64.2C42.8,-58,52,-49,57.9,-38.4C63.8,-27.8,66.4,-15.6,68.4,-2.9C70.4,9.8,71.9,23.1,66.7,33.1C61.5,43.1,49.7,49.9,38.3,56.3C26.9,62.8,15.8,68.9,3.4,71.8C-9,74.7,-22.7,74.2,-34.2,69.7C-45.7,65.1,-54.9,56.4,-61.9,45.9C-68.9,35.3,-73.6,22.9,-75.5,9.9C-77.4,-3.1,-76.5,-16.7,-71.6,-28.4C-66.7,-40.1,-57.9,-49.9,-46.9,-58.1C-35.9,-66.2,-22.8,-72.8,-8.8,-76.2C5.2,-79.6,20,-79.8,32.9,-74.8C45.7,-69.9,56.6,-59.7,65.2,-47.8C73.8,-35.9,80.1,-22.2,80,-8.8C79.8,4.7,73.2,18,66.5,29.9C59.8,41.7,53.1,52.2,43.2,59.4C33.3,66.6,20.3,70.6,6.5,73.7C-7.3,76.9,-21.8,79.2,-34.5,75.6C-47.1,72,-58,62.4,-66.2,50.8C-74.5,39.2,-80.2,25.4,-79.7,12.3C-79.2,-0.9,-72.6,-13.4,-66.3,-24.7C-59.9,-36,-53.8,-46.1,-44.5,-55.1C-35.2,-64.1,-22.8,-72.1,-8.8,-76.5C5.1,-80.9,20.6,-81.7,33.2,-76.3C45.8,-70.9,55.4,-59.2,64.3,-46.7C73.2,-34.2,81.3,-20.8,81.7,-7.1C82.1,6.6,74.9,20.7,67.2,33.2C59.5,45.7,51.3,56.7,40.6,64.2C29.8,71.8,16.5,75.9,2.9,77.9C-10.6,79.9,-24.5,79.7,-36.2,74.6C-47.9,69.5,-57.5,59.3,-65.6,47.8C-73.7,36.3,-80.4,23.3,-80.9,10.1C-81.5,-3.1,-75.9,-16.5,-69.5,-28.8C-63.1,-41.1,-55.9,-52.2,-45.4,-60.4C-34.9,-68.5,-21.2,-73.7,-6.7,-76.5C7.7,-79.3,23,-79.7,35.4,-74.3C47.8,-68.9,57.3,-57.7,65.9,-45.2C74.5,-32.8,82.3,-19.2,82.1,-5.9"
            transform="translate(100 100)"
          />
        </motion.svg>
        <motion.svg
          className="absolute top-[30%] -left-[20%] w-[70%] h-auto text-indigo-500/10 dark:text-indigo-300/10"
          viewBox="0 0 200 200"
          xmlns="http://www.w3.org/2000/svg"
          style={{ rotate: useTransform(scrollYProgress, [0, 0.5], [0, -5]) }}
        >
          <path
            fill="currentColor"
            d="M39.5,-65.3C50.3,-56.7,57.6,-44.3,64.2,-31.1C70.8,-17.9,76.7,-3.8,74.6,9.4C72.6,22.5,62.6,34.8,51.5,43.9C40.5,53.1,28.3,59.1,14.9,63.5C1.5,68,-13.2,70.8,-24.5,65.8C-35.9,60.7,-43.9,47.7,-52.3,35.2C-60.6,22.7,-69.3,10.8,-70.8,-2.2C-72.3,-15.2,-66.5,-29.3,-57.1,-39.5C-47.7,-49.7,-34.5,-56,-21.5,-62.3C-8.6,-68.6,4.2,-74.9,17.4,-74.2C30.7,-73.5,44.5,-65.8,55.6,-55.2C66.7,-44.7,77.2,-31.3,77.1,-18.2C77.1,-5,66.7,7.9,59.1,19.5C51.6,31,47,41.1,38.7,47.8C30.5,54.5,18.6,57.7,6.5,59.1C-5.7,60.5,-18,60.1,-29.9,56.3C-41.8,52.5,-53.3,45.4,-59,35.1C-64.8,24.8,-64.8,11.3,-64.7,-2.2C-64.5,-15.7,-64.2,-29.4,-58.5,-40.1C-52.8,-50.9,-41.9,-58.8,-30,-64.3C-18.1,-69.9,-5.3,-73.1,7.7,-71.9C20.7,-70.8,34.6,-65.2,45.4,-56.4C56.1,-47.6,63.8,-35.5,67.4,-22.4C71,-9.3,70.5,4,64.9,14.2C59.3,24.5,48.5,31.9,39,40.3C29.4,48.8,21.1,58.4,10,63.9C-1.1,69.4,-15,71,-26.4,66.9C-37.9,62.9,-46.9,53.3,-54.6,42.8C-62.3,32.3,-68.7,20.8,-70.4,8.6C-72.1,-3.7,-69.2,-16.8,-63.9,-29C-58.6,-41.1,-51,-52.4,-40.7,-60.3C-30.4,-68.2,-17.5,-72.9,-3.5,-74.3C10.5,-75.7,22.6,-73.8,33.3,-68.9C44,-64,53.2,-56,58.5,-45.6C63.7,-35.2,65.1,-22.3,64.7,-9.9C64.3,2.4,62.1,14.3,57.1,24.7C52.1,35.1,44.4,44.2,34.9,50.6C25.4,57,14.1,60.7,2.4,62.9C-9.4,65.1,-21.5,65.7,-33.5,62.3C-45.5,58.9,-57.4,51.3,-63.3,40.5C-69.2,29.6,-69.2,15.6,-71.2,1.1C-73.3,-13.4,-77.6,-28.3,-73.1,-40.3C-68.6,-52.3,-55.4,-61.3,-41.6,-67.5C-27.9,-73.7,-13.5,-77.1,0.2,-77.2C13.9,-77.3,27.8,-74.2,39.5,-67.3C51.2,-60.4,60.7,-49.7,66.4,-37.5C72.1,-25.3,74,-11.6,72.5,-0.1C71,11.4,66.1,20.7,60.3,29.7C54.5,38.7,47.8,47.4,38.9,53C30,58.6,18.9,61.1,6.9,65.1C-5.1,69.2,-17.9,74.9,-29.3,73.5C-40.7,72.1,-50.6,63.5,-58.5,53.4C-66.3,43.3,-72,31.5,-74.8,18.9C-77.6,6.2,-77.4,-7.4,-73.7,-19.8C-70,-32.2,-62.9,-43.5,-52.9,-52.2C-42.9,-60.9,-30.1,-67,-17.5,-71.1C-4.9,-75.2,7.5,-77.3,19.7,-76.6C31.9,-75.9,43.9,-72.5,52.9,-65.3C61.9,-58.2,67.9,-47.3,72.2,-35.6C76.5,-23.9,79,-11.3,76.8,0.5C74.6,12.2,67.7,23.4,60.3,33.9C52.9,44.3,45,54,34.8,61.2C24.6,68.4,12.3,73.1,-0.6,74.6C-13.5,76.1,-27,74.4,-39.7,69.9C-52.4,65.4,-64.4,58,-68.6,46.9C-72.9,35.8,-69.5,20.9,-69.3,7.8C-69,-5.4,-71.8,-17.8,-68.4,-28.7C-65,-39.5,-55.3,-48.8,-44.4,-56.6C-33.4,-64.4,-21.2,-70.7,-8,-72.2C5.2,-73.7,19.4,-70.3,31.1,-64.2C42.8,-58,52,-49,57.9,-38.4C63.8,-27.8,66.4,-15.6,68.4,-2.9C70.4,9.8,71.9,23.1,66.7,33.1C61.5,43.1,49.7,49.9,38.3,56.3C26.9,62.8,15.8,68.9,3.4,71.8C-9,74.7,-22.7,74.2,-34.2,69.7C-45.7,65.1,-54.9,56.4,-61.9,45.9C-68.9,35.3,-73.6,22.9,-75.5,9.9C-77.4,-3.1,-76.5,-16.7,-71.6,-28.4C-66.7,-40.1,-57.9,-49.9,-46.9,-58.1C-35.9,-66.2,-22.8,-72.8,-8.8,-76.2C5.2,-79.6,20,-79.8,32.9,-74.8C45.7,-69.9,56.6,-59.7,65.2,-47.8C73.8,-35.9,80.1,-22.2,80,-8.8C79.8,4.7,73.2,18,66.5,29.9C59.8,41.7,53.1,52.2,43.2,59.4C33.3,66.6,20.3,70.6,6.5,73.7C-7.3,76.9,-21.8,79.2,-34.5,75.6C-47.1,72,-58,62.4,-66.2,50.8C-74.5,39.2,-80.2,25.4,-79.7,12.3C-79.2,-0.9,-72.6,-13.4,-66.3,-24.7C-59.9,-36,-53.8,-46.1,-44.5,-55.1C-35.2,-64.1,-22.8,-72.1,-8.8,-76.5C5.1,-80.9,20.6,-81.7,33.2,-76.3C45.8,-70.9,55.4,-59.2,64.3,-46.7C73.2,-34.2,81.3,-20.8,81.7,-7.1C82.1,6.6,74.9,20.7,67.2,33.2C59.5,45.7,51.3,56.7,40.6,64.2C29.8,71.8,16.5,75.9,2.9,77.9C-10.6,79.9,-24.5,79.7,-36.2,74.6C-47.9,69.5,-57.5,59.3,-65.6,47.8C-73.7,36.3,-80.4,23.3,-80.9,10.1C-81.5,-3.1,-75.9,-16.5,-69.5,-28.8C-63.1,-41.1,-55.9,-52.2,-45.4,-60.4C-34.9,-68.5,-21.2,-73.7,-6.7,-76.5C7.7,-79.3,23,-79.7,35.4,-74.3C47.8,-68.9,57.3,-57.7,65.9,-45.2C74.5,-32.8,82.3,-19.2,82.1,-5.9"
            transform="translate(100 100)"
          />
        </motion.svg>
      </div>

      {/* Enhanced Floating Particles Effect - Limited height to prevent scrolling issues */}
      <div
        className="particle-container absolute inset-0 -z-10 overflow-hidden pointer-events-none opacity-40 dark:opacity-30"
        style={{ maxHeight: "100vh" }}
      >
        {Array(20)
          .fill()
          .map((_, i) => (
            <motion.div
              key={i}
              className="particle absolute rounded-full bg-gradient-to-br from-blue-400 to-indigo-400 dark:from-blue-500 dark:to-indigo-500"
              initial={{
                x: Math.random() * window.innerWidth,
                y: Math.random() * window.innerHeight,
                scale: Math.random() * 0.5 + 0.5,
              }}
              animate={{
                x: [
                  Math.random() * window.innerWidth,
                  Math.random() * window.innerWidth,
                  Math.random() * window.innerWidth,
                ],
                y: [
                  Math.random() * window.innerHeight,
                  Math.random() * window.innerHeight,
                  Math.random() * window.innerHeight,
                ],
                scale: [
                  Math.random() * 0.5 + 0.5,
                  Math.random() * 1 + 0.8,
                  Math.random() * 0.5 + 0.5,
                ],
              }}
              transition={{
                duration: Math.random() * 20 + 20,
                repeat: Infinity,
                ease: "linear",
              }}
              style={{
                width: `${Math.random() * 10 + 2}px`,
                height: `${Math.random() * 10 + 2}px`,
                opacity: Math.random() * 0.6 + 0.2,
                filter: `blur(${Math.random() * 2}px)`,
              }}
            />
          ))}
      </div>

      {/* Content Area - Outer padding container - reduced bottom padding */}
      <div className="relative w-full mx-auto px-4 sm:px-6 lg:px-8 py-12 md:py-16 pb-0 z-10 flex-grow">
        {/* Inner wrapper for width constraint */}
        <div className="w-full max-w-6xl mx-auto">
          {/* Page Header - Animated */}
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
              {/* Logo with dynamic color based on dark/light mode */}
              <img
                src="https://raw.githubusercontent.com/Steinbeck-Lab/cheminformatics-microservice/refs/heads/main/public/img/logo.svg"
                alt="Cheminformatics Microservice Logo"
                className="w-full filter dark:brightness-100 brightness-75 transition-all duration-500 transform hover:scale-105 drop-shadow-2xl"
              />
              <div className="absolute inset-0 bg-gradient-to-r from-blue-500/0 via-blue-500/10 to-blue-500/0 opacity-0 group-hover:opacity-100 animate-shimmer"></div>
            </motion.div>

            <motion.h1
              variants={headerItemVariants}
              className="text-4xl md:text-5xl lg:text-6xl font-bold mb-6"
            >
              <span className="text-transparent bg-clip-text bg-gradient-to-r from-blue-600 to-indigo-600 dark:from-blue-400 dark:to-indigo-400">
                Unifying Access to Open Cheminformatics Toolkits
              </span>
            </motion.h1>

            <motion.p
              variants={headerItemVariants}
              className="text-lg sm:text-xl text-slate-700 dark:text-slate-300 max-w-3xl mx-auto leading-relaxed"
            >
              A modern microservice platform providing seamless access to
              multiple open-source cheminformatics toolkits through a unified
              REST API interface.
            </motion.p>
          </motion.div>

          {/* Main Overview - Enhanced with better spacing and shadows */}
          <motion.div
            className="grid grid-cols-1 md:grid-cols-2 gap-8 mb-16"
            variants={contentContainerVariant}
            initial="hidden"
            animate="visible"
          >
            <motion.div
              className="glass rounded-xl shadow-xl border border-slate-200/50 dark:border-slate-700/30 overflow-hidden h-full transition-all duration-500 hover:shadow-2xl hover:border-blue-200/50 dark:hover:border-blue-700/30"
              whileHover={{ y: -10 }}
              transition={{ duration: 0.3 }}
            >
              <div className="p-8 md:p-10">
                <div className="flex items-center mb-6">
                  <FaConnectdevelop className="text-3xl text-blue-500 dark:text-blue-400 mr-4" />
                  <h2 className="text-2xl font-bold text-slate-900 dark:text-white">
                    Overview
                  </h2>
                </div>
                <p className="text-slate-700 dark:text-slate-300 mb-4 leading-relaxed">
                  The Cheminformatics Microservice offers a collection of
                  versatile functions accessible via REST endpoints that can
                  handle chemical data and perform various cheminformatics
                  tasks. These tasks include (but are not limited to) generating
                  chemical structure depictions, 3D conformers, descriptors,
                  IUPAC names, and converting between machine-readable formats.
                </p>
                <p className="text-slate-700 dark:text-slate-300 leading-relaxed">
                  Researchers and developers can effectively access open-source
                  cheminformatics toolkits such as{" "}
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
                  through this microservice and extend them easily to suit their
                  specific needs.
                </p>
              </div>
            </motion.div>

            <motion.div
              className="glass rounded-xl shadow-xl border border-slate-200/50 dark:border-slate-700/30 overflow-hidden h-full transition-all duration-500 hover:shadow-2xl hover:border-blue-200/50 dark:hover:border-blue-700/30"
              whileHover={{ y: -10 }}
              transition={{ duration: 0.3 }}
            >
              <div className="p-8 md:p-10">
                <div className="flex items-center mb-6">
                  <FaFlask className="text-3xl text-indigo-500 dark:text-indigo-400 mr-4" />
                  <h2 className="text-2xl font-bold text-slate-900 dark:text-white">
                    Key Features
                  </h2>
                </div>
                <ul className="space-y-3 text-slate-700 dark:text-slate-300">
                  <li className="flex items-start">
                    <span className="text-blue-600 dark:text-blue-400 mr-2 text-xl">
                      •
                    </span>
                    <span>
                      Unified access to multiple cheminformatics toolkits
                      (RDKit, CDK, OpenBabel)
                    </span>
                  </li>
                  <li className="flex items-start">
                    <span className="text-blue-600 dark:text-blue-400 mr-2 text-xl">
                      •
                    </span>
                    <span>High-quality 2D and 3D molecular visualizations</span>
                  </li>
                  <li className="flex items-start">
                    <span className="text-blue-600 dark:text-blue-400 mr-2 text-xl">
                      •
                    </span>
                    <span>Comprehensive molecular descriptor calculations</span>
                  </li>
                  <li className="flex items-start">
                    <span className="text-blue-600 dark:text-blue-400 mr-2 text-xl">
                      •
                    </span>
                    <span>
                      Format conversions (SMILES, InChI, InChIKey, SELFIES,
                      etc.)
                    </span>
                  </li>
                  <li className="flex items-start">
                    <span className="text-blue-600 dark:text-blue-400 mr-2 text-xl">
                      •
                    </span>
                    <span>Structure standardization and validation</span>
                  </li>
                  <li className="flex items-start">
                    <span className="text-blue-600 dark:text-blue-400 mr-2 text-xl">
                      •
                    </span>
                    <span>
                      Optical Chemical Structure Recognition (OCSR) using
                      DECIMER
                    </span>
                  </li>
                  <li className="flex items-start">
                    <span className="text-blue-600 dark:text-blue-400 mr-2 text-xl">
                      •
                    </span>
                    <span>Easy deployment via Docker containers</span>
                  </li>
                </ul>
              </div>
            </motion.div>
          </motion.div>

          {/* Integrated Tools Section - With logo integration */}
          <motion.div
            className="glass rounded-xl shadow-xl border border-slate-200/50 dark:border-slate-700/30 overflow-hidden mb-16 transition-all duration-500 hover:shadow-2xl hover:border-blue-200/50 dark:hover:border-blue-700/30"
            variants={contentContainerVariant}
            initial="hidden"
            animate="visible"
            whileHover={{ y: -5 }}
            transition={{ duration: 0.3 }}
          >
            <div className="p-8 md:p-10">
              <div className="flex items-center mb-8">
                <HiOutlineCode className="text-3xl text-purple-500 dark:text-purple-400 mr-4" />
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
                      <div className="flex-shrink-0 w-16 h-16 bg-white dark:bg-slate-700 rounded-lg flex items-center justify-center p-2 shadow-md group-hover:shadow-lg transition-shadow">
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
                          A collection of cheminformatics and machine learning
                          tools
                        </p>
                      </div>
                    </motion.div>

                    {/* CDK */}
                    <motion.div
                      className="flex items-center space-x-4 group"
                      whileHover={{ x: 5 }}
                      transition={{ duration: 0.2 }}
                    >
                      <div className="flex-shrink-0 w-16 h-16 bg-white dark:bg-slate-700 rounded-lg flex items-center justify-center p-2 shadow-md group-hover:shadow-lg transition-shadow">
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
                          Chemistry Development Kit, a Java library for
                          structural cheminformatics
                        </p>
                      </div>
                    </motion.div>

                    {/* OpenBabel */}
                    <motion.div
                      className="flex items-center space-x-4 group"
                      whileHover={{ x: 5 }}
                      transition={{ duration: 0.2 }}
                    >
                      <div className="flex-shrink-0 w-16 h-16 bg-white dark:bg-slate-700 rounded-lg flex items-center justify-center p-2 shadow-md group-hover:shadow-lg transition-shadow">
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
                          A chemical toolbox designed to search, convert,
                          analyze, or store data
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
                      <div className="flex-shrink-0 w-12 h-12 bg-white dark:bg-slate-700 rounded-lg flex items-center justify-center p-2 shadow-md mr-4">
                        <img
                          src="https://github.com/Steinbeck-Lab/cheminf-jena-logos/blob/main/DECIMER/DECIMER-1.png?raw=true"
                          alt="DECIMER Logo"
                          className="w-full h-full object-contain"
                        />
                      </div>
                      <div>
                        <span className="font-semibold">DECIMER</span>
                        <p className="text-sm">
                          Deep learning for Chemical Image Recognition
                        </p>
                      </div>
                    </motion.li>
                    <motion.li
                      className="flex items-start"
                      whileHover={{ x: 5 }}
                      transition={{ duration: 0.2 }}
                    >
                      <div className="flex-shrink-0 w-12 h-12 bg-white dark:bg-slate-700 rounded-lg flex items-center justify-center p-2 shadow-md mr-4">
                        <img
                          src="https://raw.githubusercontent.com/epam/ketcher/refs/heads/master/demo/public/favicon.ico"
                          alt="Ketcher Logo"
                          className="w-full h-full object-contain"
                        />
                      </div>
                      <div>
                        <span className="font-semibold">Ketcher</span>
                        <p className="text-sm">
                          Web-based chemical structure editor
                        </p>
                      </div>
                    </motion.li>
                    <motion.li
                      className="flex items-start"
                      whileHover={{ x: 5 }}
                      transition={{ duration: 0.2 }}
                    >
                      <div className="flex-shrink-0 w-12 h-12 bg-white dark:bg-slate-700 rounded-lg flex items-center justify-center p-2 shadow-md mr-4">
                        <img
                          src="https://upload.wikimedia.org/wikipedia/commons/b/b6/PubChem_logo.svg"
                          alt="PubChem Logo"
                          className="w-full h-full object-contain"
                        />
                      </div>
                      <div>
                        <span className="font-semibold">
                          PubChem PUG REST API
                        </span>
                        <p className="text-sm">
                          For compound search and retrieval
                        </p>
                      </div>
                    </motion.li>
                    <motion.li
                      className="flex items-start"
                      whileHover={{ x: 5 }}
                      transition={{ duration: 0.2 }}
                    >
                      <div className="flex-shrink-0 w-12 h-12 bg-white dark:bg-slate-700 rounded-lg flex items-center justify-center p-2 shadow-md mr-4">
                        <img
                          src="https://github.com/IUPAC-InChI/InChI-Web-Demo/blob/main/pages/img/InChI-logo2.png?raw=true"
                          alt="InChI Logo"
                          className="w-full h-full object-contain"
                        />
                      </div>
                      <div>
                        <span className="font-semibold">
                          InChI Web Application
                        </span>
                        <p className="text-sm">
                          InChI evaluation from chemical depictions
                        </p>
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
            <div className="text-center">
              <motion.a
                href="https://naturalproducts.net"
                target="_blank"
                rel="noopener noreferrer"
                className="inline-flex items-center justify-center px-8 py-4 bg-gradient-to-r from-blue-600 to-indigo-700 hover:from-blue-500 hover:to-indigo-600 text-white font-semibold rounded-lg shadow-lg"
                variants={buttonVariant}
                whileHover="hover"
              >
                <FaAtom className="mr-3 h-5 w-5" />
                <span className="text-lg">Natural Products Online</span>
                <HiOutlineExternalLink className="ml-3 h-5 w-5" />
              </motion.a>
            </div>
          </motion.div>

          {/* Public Instance */}
          <motion.div
            className="glass rounded-xl shadow-xl border border-slate-200/50 dark:border-slate-700/30 overflow-hidden mb-12 transition-all duration-500 hover:shadow-2xl hover:border-blue-200/50 dark:hover:border-blue-700/30"
            variants={contentContainerVariant}
            initial="hidden"
            animate="visible"
            whileHover={{ y: -5 }}
            transition={{ duration: 0.3 }}
          >
            <div className="p-8 md:p-10">
              <div className="flex items-center mb-6">
                <HiOutlineDatabase className="text-3xl text-green-500 dark:text-green-400 mr-4" />
                <h2 className="text-2xl font-bold text-slate-900 dark:text-white">
                  Public Instance
                </h2>
              </div>
              <div className="bg-blue-50/70 dark:bg-blue-900/20 rounded-xl p-6 border border-blue-200/70 dark:border-blue-800/50 mb-4">
                <p className="text-slate-700 dark:text-slate-300 mb-6">
                  A public instance of the Cheminformatics Microservice is
                  hosted at{" "}
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
                  <motion.a
                    href="https://api.naturalproducts.net/latest/docs"
                    target="_blank"
                    rel="noopener noreferrer"
                    className="inline-flex items-center px-8 py-4 rounded-lg bg-gradient-to-r from-blue-600 to-indigo-600 dark:from-blue-500 dark:to-indigo-500 text-white font-medium transition-all duration-300 hover:shadow-lg hover:from-blue-700 hover:to-indigo-700"
                    whileHover={{ y: -5, scale: 1.03 }}
                    transition={{ duration: 0.2 }}
                  >
                    <HiOutlineExternalLink className="mr-2 h-5 w-5" />
                    api.naturalproducts.net/latest/docs
                  </motion.a>
                </div>
              </div>
            </div>
          </motion.div>

          {/* Citation Guidelines */}
          <motion.div
            className="glass rounded-xl shadow-xl border border-slate-200/50 dark:border-slate-700/30 overflow-hidden mb-12 transition-all duration-500 hover:shadow-2xl hover:border-blue-200/50 dark:hover:border-blue-700/30"
            variants={contentContainerVariant}
            initial="hidden"
            animate="visible"
            whileHover={{ y: -5 }}
            transition={{ duration: 0.3 }}
          >
            <div className="p-8 md:p-10">
              <div className="flex items-center mb-6">
                <HiOutlineDocumentText className="text-3xl text-cyan-500 dark:text-cyan-400 mr-4" />
                <h2 className="text-2xl font-bold text-slate-900 dark:text-white">
                  Citation Guidelines
                </h2>
              </div>
              <p className="text-slate-700 dark:text-slate-300 mb-6 leading-relaxed">
                It is strongly recommended that users cite both the software and
                the corresponding paper when utilizing this work. By providing
                appropriate citations for the Cheminformatics Microservice,
                users gain a convenient means to precisely track the original
                source of the utilized source code and data.
              </p>

              <div className="grid grid-cols-1 md:grid-cols-2 gap-8 mt-6">
                <CitationCard
                  title="Citing the paper:"
                  citation="Chandrasekhar, V., Sharma, N., Schaub, J. et al. Cheminformatics Microservice: unifying access to open cheminformatics toolkits. J Cheminform 15, 98 (2023). https://doi.org/10.1186/s13321-023-00762-4"
                />
                <CitationCard
                  title="Citing the software:"
                  citation="Venkata, C., Sharma, N., Schaub, J., Steinbeck, C., & Rajan, K. (2023). cheminformatics-microservice (Version v2.6.0) [Computer software]. https://doi.org/10.5281/zenodo.13867839"
                />
              </div>
            </div>
          </motion.div>

          {/* License Information */}
          <motion.div
            className="glass rounded-xl shadow-xl border border-slate-200/50 dark:border-slate-700/30 overflow-hidden mb-12 transition-all duration-500 hover:shadow-2xl hover:border-blue-200/50 dark:hover:border-blue-700/30"
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
                  <p className="text-slate-700 dark:text-slate-300">
                    Copyright © {currentYear}
                  </p>
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
                    Other components are used according to their respective
                    licenses.
                  </p>
                </motion.div>
              </div>
            </div>
          </motion.div>

          {/* Acknowledgements */}
          <motion.div
            className="glass rounded-xl shadow-xl border border-slate-200/50 dark:border-slate-700/30 overflow-hidden mb-8"
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
                Cheminformatics Microservice and{" "}
                <a
                  href="https://naturalproducts.net"
                  target="_blank"
                  rel="noopener noreferrer"
                  className="text-blue-600 dark:text-blue-400 hover:underline"
                >
                  Natural Products Online
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
                    className="link-style"
                  >
                    DFG
                  </a>{" "}
                  (German Research Foundation) under{" "}
                  <a
                    href="https://nfdi4chem.de/"
                    target="_blank"
                    rel="noopener noreferrer"
                    className="link-style"
                  >
                    NFDI4Chem
                  </a>{" "}
                  (Project: <strong>441958208</strong>) and{" "}
                  <a
                    href="https://www.chembiosys.de/en/"
                    target="_blank"
                    rel="noopener noreferrer"
                    className="link-style"
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
                    <motion.a
                      key={funder.name}
                      href={funder.url}
                      target="_blank"
                      rel="noopener noreferrer"
                      className="block bg-white dark:bg-slate-800/60 p-3 rounded-lg shadow-md hover:shadow-lg transition-all duration-300 w-auto h-16 flex items-center justify-center" // Fixed height
                      whileHover={{ y: -5, scale: 1.05 }}
                      transition={{ duration: 0.2 }}
                    >
                      <img
                        src={funder.logo}
                        alt={`${funder.name} Logo`}
                        className="max-h-12 object-contain"
                      />{" "}
                      {/* Max height */}
                    </motion.a>
                  ))}
                </div>
              </div>
            </div>
          </motion.div>
        </div>
      </div>

      {/* Fix for potential extra space */}
      <div style={{ height: 0, clear: "both", overflow: "hidden" }} />

      {/* Custom CSS */}
      <style jsx global>{`
        .glass {
          @apply bg-white/80 dark:bg-slate-800/80 backdrop-blur-xl dark:backdrop-blur-2xl;
        }

        .text-gradient {
          @apply bg-clip-text text-transparent bg-gradient-to-r from-blue-600 to-indigo-600 dark:from-blue-400 dark:to-indigo-400;
        }

        /* Animation for shimmer effect */
        @keyframes shimmer {
          from {
            transform: translateX(-100%);
          }
          to {
            transform: translateX(100%);
          }
        }

        .animate-shimmer {
          animation: shimmer 2s infinite;
        }

        /* Custom scrollbar */
        ::-webkit-scrollbar {
          width: 10px;
        }

        ::-webkit-scrollbar-track {
          background: rgba(241, 245, 249, 0.1);
        }

        ::-webkit-scrollbar-thumb {
          background: #94a3b8;
          border-radius: 5px;
        }

        ::-webkit-scrollbar-thumb:hover {
          background: #64748b;
        }

        /* Smooth hover transitions for all links */
        a {
          transition: all 0.2s ease-in-out;
        }

        a:hover {
          transform: translateY(-2px);
        }

        /* Fix for extra whitespace */
        html,
        body {
          height: 100%;
          overflow-x: hidden;
        }

        #root {
          min-height: 100%;
          display: flex;
          flex-direction: column;
        }

        .particle-container {
          max-height: 100vh;
          overflow: hidden;
        }
      `}</style>
    </motion.div>
  );
};

export default AboutPage;
