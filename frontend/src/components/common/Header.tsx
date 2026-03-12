import React, { useState, useEffect } from "react";
import { Link, useLocation } from "react-router-dom";
import { useAppContext } from "../../context/AppContext";
import Navigation from "./Navigation";
import { HiOutlineMoon, HiOutlineSun, HiOutlineMenu, HiOutlineX } from "react-icons/hi";
import { motion, AnimatePresence, LayoutGroup } from "framer-motion";

// --- Enhanced Animation Variants ---
const headerVariants = {
  hidden: { opacity: 0, y: -20 },
  visible: {
    opacity: 1,
    y: 0,
    transition: {
      duration: 0.7,
      ease: [0.22, 1, 0.36, 1], // Custom cubic-bezier for smoother animation
      when: "beforeChildren",
      staggerChildren: 0.1,
    },
  },
};

const headerContentVariants = {
  hidden: { opacity: 0, y: -15, scale: 0.95 },
  visible: (i = 1) => ({
    opacity: 1,
    y: 0,
    scale: 1,
    transition: {
      duration: 0.6,
      delay: 0.05 * i,
      ease: [0.2, 0.65, 0.3, 0.9],
      scale: {
        type: "spring",
        stiffness: 400,
        damping: 25,
      },
    },
  }),
};

const mobileMenuVariants = {
  hidden: {
    opacity: 0,
    height: 0,
    y: -10,
    transformOrigin: "top center",
  },
  visible: {
    opacity: 1,
    height: "auto",
    y: 0,
    transition: {
      duration: 0.4,
      ease: [0.4, 0, 0.2, 1],
      staggerChildren: 0.05,
      when: "beforeChildren",
    },
  },
  exit: {
    opacity: 0,
    height: 0,
    y: -10,
    transition: {
      duration: 0.3,
      ease: [0.4, 0, 1, 1],
      when: "afterChildren",
    },
  },
};

const menuItemVariants = {
  hidden: { opacity: 0, x: -20 },
  visible: {
    opacity: 1,
    x: 0,
    transition: {
      type: "spring",
      stiffness: 400,
      damping: 30,
    },
  },
  exit: {
    opacity: 0,
    x: -20,
    transition: {
      duration: 0.2,
    },
  },
};

// Spring transition for the pill toggle
const pillSwitchTransition = {
  type: "spring",
  stiffness: 600,
  damping: 30,
};

// --- Component ---
const Header = () => {
  const { isDarkMode, toggleDarkMode } = useAppContext();
  const [isMenuOpen, setIsMenuOpen] = useState(false);
  const [isScrolled, setIsScrolled] = useState(false);
  const location = useLocation();

  // Close menu on route change
  useEffect(() => {
    setIsMenuOpen(false);
  }, [location]);

  // Check scroll position for glass effect intensity
  useEffect(() => {
    const handleScroll = () => {
      setIsScrolled(window.scrollY > 10);
    };

    window.addEventListener("scroll", handleScroll);
    return () => window.removeEventListener("scroll", handleScroll);
  }, []);

  const toggleMenu = () => {
    setIsMenuOpen(!isMenuOpen);
  };

  // Use higher quality logo URLs
  const logoDark =
    "https://raw.githubusercontent.com/Steinbeck-Lab/cheminformatics-microservice/main/public/img/logo_small_inverted.png";
  const logoLight =
    "https://raw.githubusercontent.com/Steinbeck-Lab/cheminformatics-microservice/main/public/img/logo_small.png";

  // Calculate header background classes based on scroll state - always white in light mode
  const headerBgClasses = isScrolled
    ? `backdrop-blur-lg ${isDarkMode ? "bg-slate-900/90" : "bg-white"} shadow-lg ${isDarkMode ? "shadow-slate-900/30" : "shadow-slate-200/50"}`
    : `backdrop-blur-sm ${isDarkMode ? "bg-slate-900/80" : "bg-white"} shadow-md ${isDarkMode ? "shadow-slate-900/20" : "shadow-slate-200/30"}`;

  return (
    <motion.header
      className={`sticky top-0 z-50 border-b transition-all duration-300 ${headerBgClasses} ${isDarkMode ? "border-slate-800/40" : "border-slate-200/70"}`}
      variants={headerVariants}
      initial="hidden"
      animate="visible"
    >
      <div className="w-full mx-auto px-4 sm:px-6 lg:px-8">
        <div className="flex items-center h-16 justify-between md:justify-center">
          {/* Logo and title */}
          <motion.div
            className="flex-shrink-0 flex items-center"
            variants={headerContentVariants}
            custom={0}
          >
            <Link to="/" className="flex items-center group" aria-label="Homepage">
              <motion.div
                className="overflow-hidden rounded-md mr-3"
                whileHover={{
                  scale: 1.08,
                  rotate: [0, 2, 0, -2, 0],
                  transition: { rotate: { repeat: Infinity, duration: 1.5 } },
                }}
              >
                <motion.img
                  src={isDarkMode ? logoDark : logoLight}
                  alt="Cheminformatics Microservice Logo"
                  className="h-9 w-auto p-0.5"
                  initial={{ scale: 0.9, opacity: 0 }}
                  animate={{ scale: 1, opacity: 1 }}
                  transition={{ duration: 0.5 }}
                />
              </motion.div>

              <div className="flex flex-col justify-center">
                <motion.span
                  className={`font-bold text-lg sm:text-xl leading-tight ${
                    isDarkMode
                      ? "bg-gradient-to-r from-sky-300 to-blue-400 text-transparent bg-clip-text"
                      : "bg-gradient-to-r from-sky-600 to-indigo-600 text-transparent bg-clip-text"
                  }`}
                  whileHover={{ scale: 1.03 }}
                >
                  Cheminformatics
                </motion.span>
                <motion.span
                  className={`hidden md:inline text-[11px] leading-tight -mt-0.5 transition-colors ${
                    isDarkMode
                      ? "text-sky-300 group-hover:text-sky-200"
                      : "text-indigo-500 group-hover:text-indigo-700"
                  }`}
                  initial={{ y: 5, opacity: 0 }}
                  animate={{ y: 0, opacity: 1 }}
                  transition={{ delay: 0.2, duration: 0.4 }}
                >
                  Microservices
                </motion.span>
              </div>
            </Link>
          </motion.div>

          {/* Desktop Navigation */}
          <motion.div
            className="hidden md:flex items-center ml-6 lg:ml-10"
            variants={headerContentVariants}
            custom={1}
          >
            <Navigation />
          </motion.div>

          {/* Theme Toggle */}
          <motion.div
            className="hidden md:flex items-center ml-6 lg:ml-10"
            variants={headerContentVariants}
            custom={2}
          >
            {/* Pill Theme Toggle */}
            <div
              className={`relative flex items-center w-[62px] h-8 rounded-full p-1 cursor-pointer transition-all duration-300 ease-in-out ${
                isDarkMode
                  ? "bg-gradient-to-r from-slate-800 to-slate-700 hover:shadow-inner hover:shadow-slate-900"
                  : "bg-gradient-to-r from-sky-100 to-indigo-100 hover:shadow-md hover:shadow-indigo-200/50"
              }`}
              onClick={toggleDarkMode}
              aria-label={isDarkMode ? "Switch to light mode" : "Switch to dark mode"}
              role="switch"
              aria-checked={isDarkMode}
            >
              <LayoutGroup>
                <motion.div
                  className={`absolute z-10 h-6 w-6 rounded-full shadow-lg ${
                    isDarkMode
                      ? "bg-gradient-to-br from-slate-800 to-slate-900"
                      : "bg-gradient-to-br from-white to-sky-50"
                  }`}
                  layout
                  transition={pillSwitchTransition}
                  style={{
                    left: isDarkMode ? "auto" : "4px",
                    right: isDarkMode ? "4px" : "auto",
                  }}
                />
              </LayoutGroup>
              <div className="relative z-0 flex justify-between w-full px-1">
                <motion.div whileHover={{ rotate: 45 }} transition={{ duration: 0.3 }}>
                  <HiOutlineSun
                    className={`h-4 w-4 transition-colors ${isDarkMode ? "text-slate-500" : "text-yellow-500"}`}
                  />
                </motion.div>
                <motion.div whileHover={{ rotate: -45 }} transition={{ duration: 0.3 }}>
                  <HiOutlineMoon
                    className={`h-4 w-4 transition-colors ${isDarkMode ? "text-yellow-300" : "text-slate-400"}`}
                  />
                </motion.div>
              </div>
            </div>
          </motion.div>

          {/* Mobile menu button */}
          <div className="flex items-center md:hidden">
            <motion.button
              onClick={toggleMenu}
              className={`p-2 rounded-full transition-colors duration-200 ${
                isDarkMode
                  ? "text-slate-300 hover:text-sky-300 hover:bg-slate-700/70"
                  : "text-slate-600 hover:text-sky-600 hover:bg-slate-200/70"
              } focus:outline-none focus:ring-2 focus:ring-offset-2 focus:ring-offset-transparent focus:ring-sky-500`}
              aria-label="Toggle mobile menu"
              aria-expanded={isMenuOpen}
              whileHover={{ scale: 1.05 }}
              whileTap={{ scale: 0.95 }}
            >
              <AnimatePresence mode="wait" initial={false}>
                {isMenuOpen ? (
                  <motion.div
                    key="close"
                    initial={{ rotate: -90, opacity: 0 }}
                    animate={{ rotate: 0, opacity: 1 }}
                    exit={{ rotate: 90, opacity: 0 }}
                    transition={{ duration: 0.2 }}
                  >
                    <HiOutlineX className="h-6 w-6" />
                  </motion.div>
                ) : (
                  <motion.div
                    key="open"
                    initial={{ rotate: 90, opacity: 0 }}
                    animate={{ rotate: 0, opacity: 1 }}
                    exit={{ rotate: -90, opacity: 0 }}
                    transition={{ duration: 0.2 }}
                  >
                    <HiOutlineMenu className="h-6 w-6" />
                  </motion.div>
                )}
              </AnimatePresence>
            </motion.button>
          </div>
        </div>
      </div>

      {/* Mobile menu */}
      <AnimatePresence>
        {isMenuOpen && (
          <motion.div
            className={`md:hidden absolute top-full left-0 right-0 overflow-hidden backdrop-blur-md ${
              isDarkMode
                ? "bg-slate-900/95 border-t border-slate-800/40 shadow-lg shadow-slate-900/50"
                : "bg-white border-t border-slate-200/70 shadow-lg shadow-slate-200/50"
            }`}
            variants={mobileMenuVariants}
            initial="hidden"
            animate="visible"
            exit="exit"
          >
            <motion.div className="px-3 pt-2 pb-4 space-y-1" variants={menuItemVariants}>
              <Navigation
                isMobile={true}
                closeMenu={() => setIsMenuOpen(false)}
                isAnimated={true}
                menuItemVariants={menuItemVariants}
              />

              {/* Add mobile theme toggle */}
              <div className="mt-4 pt-3 pb-1 border-t border-dashed border-slate-300/20 flex justify-center">
                <motion.div className="flex items-center gap-3" variants={menuItemVariants}>
                  <span className={`text-sm ${isDarkMode ? "text-slate-400" : "text-slate-500"}`}>
                    Theme:
                  </span>
                  <div
                    className={`relative flex items-center w-[62px] h-8 rounded-full p-1 cursor-pointer transition-all duration-300 ease-in-out ${
                      isDarkMode
                        ? "bg-gradient-to-r from-slate-800 to-slate-700"
                        : "bg-gradient-to-r from-sky-100 to-indigo-100"
                    }`}
                    onClick={toggleDarkMode}
                    aria-label={isDarkMode ? "Switch to light mode" : "Switch to dark mode"}
                    role="switch"
                    aria-checked={isDarkMode}
                  >
                    <LayoutGroup>
                      <motion.div
                        className={`absolute z-10 h-6 w-6 rounded-full shadow-lg ${
                          isDarkMode
                            ? "bg-gradient-to-br from-slate-800 to-slate-900"
                            : "bg-gradient-to-br from-white to-sky-50"
                        }`}
                        layout
                        transition={pillSwitchTransition}
                        style={{
                          left: isDarkMode ? "auto" : "4px",
                          right: isDarkMode ? "4px" : "auto",
                        }}
                      />
                    </LayoutGroup>
                    <div className="relative z-0 flex justify-between w-full px-1">
                      <HiOutlineSun
                        className={`h-4 w-4 transition-colors ${isDarkMode ? "text-slate-500" : "text-yellow-500"}`}
                      />
                      <HiOutlineMoon
                        className={`h-4 w-4 transition-colors ${isDarkMode ? "text-yellow-300" : "text-slate-400"}`}
                      />
                    </div>
                  </div>
                </motion.div>
              </div>
            </motion.div>
          </motion.div>
        )}
      </AnimatePresence>
    </motion.header>
  );
};

export default Header;
