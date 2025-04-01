// Description: Header component with logo, navigation, and theme toggle.
import React, { useState, useEffect } from 'react';
import { Link, useLocation } from 'react-router-dom';
import { useAppContext } from '../../context/AppContext'; // Assuming context exists
import Navigation from './Navigation'; // Assuming Navigation component exists
import { HiOutlineMoon, HiOutlineSun, HiOutlineMenu, HiOutlineX } from 'react-icons/hi';
import { motion, AnimatePresence, LayoutGroup } from 'framer-motion';

// --- Animation Variants ---
const headerVariants = {
  hidden: { opacity: 0, y: -10 },
  visible: { opacity: 1, y: 0, transition: { duration: 0.5, ease: "easeOut" } }
};

// Stagger children slightly less aggressively
const headerContentVariants = {
  hidden: { opacity: 0, y: -10 },
  visible: (i = 1) => ({
    opacity: 1,
    y: 0,
    transition: { duration: 0.6, delay: 0.1 + i * 0.08, ease: [0.2, 0.65, 0.3, 0.9] }
  })
};

const mobileMenuVariants = {
  hidden: { opacity: 0, height: 0, y: -10 },
  visible: { opacity: 1, height: 'auto', y: 0, transition: { duration: 0.3, ease: [0.4, 0, 0.2, 1] } },
  exit: { opacity: 0, height: 0, y: -10, transition: { duration: 0.25, ease: "easeOut" } }
};

// Spring transition for the pill toggle knob
const pillSwitchTransition = {
  type: "spring",
  stiffness: 500,
  damping: 40
};

// --- Component ---
const Header = () => {
  const { isDarkMode, toggleDarkMode } = useAppContext();
  const [isMenuOpen, setIsMenuOpen] = useState(false);
  const location = useLocation();

  useEffect(() => {
    setIsMenuOpen(false); // Close menu on route change
  }, [location]);

  const toggleMenu = () => {
    setIsMenuOpen(!isMenuOpen);
  };

  const logoDark = "https://raw.githubusercontent.com/Steinbeck-Lab/cheminformatics-microservice/main/public/img/logo_small_inverted.png";
  const logoLight = "https://raw.githubusercontent.com/Steinbeck-Lab/cheminformatics-microservice/main/public/img/logo_small.png";

  return (
    <motion.header
      className="glass sticky top-0 z-50 shadow-md dark:shadow-lg dark:shadow-slate-900/50 border-b border-slate-200/70 dark:border-slate-700/40"
      variants={headerVariants}
      initial="hidden"
      animate="visible"
    >
      <div className="w-full mx-auto px-4 sm:px-6 lg:px-8">
        {/* Main flex container:
            - Mobile (<md): justify-between (logo left, button right)
            - Desktop (md+): justify-center (all visible items centered)
        */}
        <div className="flex items-center h-16 justify-between md:justify-center">

          {/* Logo and title (Always visible, part of flex flow) */}
          <motion.div
            className="flex-shrink-0 flex items-center" // Keep shrink-0
            variants={headerContentVariants}
            custom={0} // Stagger index 0
          >
            <Link to="/" className="flex items-center group" aria-label="Homepage">
              <motion.img
                src={isDarkMode ? logoDark : logoLight}
                alt="Cheminformatics Microservice Logo"
                className="h-8 w-auto rounded-md p-0.5 mr-2 sm:mr-3"
                whileHover={{ scale: 1.05 }}
                transition={{ type: 'spring', stiffness: 300 }}
              />
              <div className="flex flex-col justify-center">
                <span className="font-bold text-lg sm:text-xl text-gradient leading-tight">
                  Cheminformatics
                </span>
                <span className="hidden md:inline text-[11px] text-slate-500 dark:text-blue-300 group-hover:text-slate-700 dark:group-hover:text-blue-200 transition-colors leading-tight -mt-0.5">
                  Microservices
                </span>
              </div>
            </Link>
          </motion.div>

          {/* Desktop Navigation (Hidden on mobile, part of centered group on desktop) */}
          <motion.div
            // Added margin left for spacing from logo ONLY on desktop
            className="hidden md:flex items-center ml-6 lg:ml-10"
            variants={headerContentVariants}
            custom={1} // Stagger index 1
          >
            <Navigation />
          </motion.div>

          {/* Theme Toggle (Hidden on mobile, part of centered group on desktop) */}
          {/* Use a wrapper div for layout and animation */}
          <motion.div
            // Added margin left for spacing from nav ONLY on desktop
            className="hidden md:flex items-center ml-6 lg:ml-10"
            variants={headerContentVariants}
            custom={2} // Stagger index 2
          >
            {/* --- Pill Theme Toggle --- */}
            <div
              className={`relative flex items-center w-[62px] h-8 rounded-full p-1 cursor-pointer transition-colors duration-300 ease-in-out ${isDarkMode ? 'bg-slate-700 hover:bg-slate-600' : 'bg-sky-100 hover:bg-sky-200'
                }`}
              onClick={toggleDarkMode}
              aria-label={isDarkMode ? 'Switch to light mode' : 'Switch to dark mode'}
              role="switch"
              aria-checked={isDarkMode}
            >
              <LayoutGroup>
                <motion.div
                  className={`absolute z-10 h-6 w-6 rounded-full shadow-md ${isDarkMode ? 'bg-slate-900' : 'bg-white'
                    }`}
                  layout
                  transition={pillSwitchTransition}
                  style={{ left: isDarkMode ? 'auto' : '4px', right: isDarkMode ? '4px' : 'auto' }}
                />
              </LayoutGroup>
              <div className="relative z-0 flex justify-between w-full px-1">
                <HiOutlineSun className={`h-4 w-4 transition-colors ${isDarkMode ? 'text-slate-500' : 'text-yellow-500'}`} />
                <HiOutlineMoon className={`h-4 w-4 transition-colors ${isDarkMode ? 'text-yellow-300' : 'text-slate-400'}`} />
              </div>
            </div>
            {/* --- End Pill Theme Toggle --- */}
          </motion.div>


          {/* Mobile menu button (Visible mobile only, pushed right by justify-between) */}
          <div className="flex items-center md:hidden"> {/* Wrapper div ensures it's treated as a single element for justify-between */}
            <button
              onClick={toggleMenu}
              className="p-2 rounded-full text-slate-600 dark:text-slate-300 hover:text-sky-600 dark:hover:text-sky-300 hover:bg-slate-200/70 dark:hover:bg-slate-700/70 focus:outline-none focus:ring-2 focus:ring-offset-2 focus:ring-offset-transparent dark:focus:ring-offset-transparent focus:ring-sky-500 dark:focus:ring-sky-400 transition-colors duration-200"
              aria-label="Toggle mobile menu"
              aria-expanded={isMenuOpen}
            >
              <AnimatePresence mode="wait" initial={false}>
                {isMenuOpen ? (
                  <motion.div key="close" initial={{ rotate: -90, opacity: 0 }} animate={{ rotate: 0, opacity: 1 }} exit={{ rotate: 90, opacity: 0 }} transition={{ duration: 0.2 }}>
                    <HiOutlineX className="h-6 w-6" />
                  </motion.div>
                ) : (
                  <motion.div key="open" initial={{ rotate: 90, opacity: 0 }} animate={{ rotate: 0, opacity: 1 }} exit={{ rotate: -90, opacity: 0 }} transition={{ duration: 0.2 }}>
                    <HiOutlineMenu className="h-6 w-6" />
                  </motion.div>
                )}
              </AnimatePresence>
            </button>
          </div>

        </div> {/* End Main Flex Container */}
      </div> {/* End Padding Container */}

      {/* Mobile menu */}
      <AnimatePresence>
        {isMenuOpen && (
          <motion.div
            className="md:hidden absolute top-full left-0 right-0 shadow-lg glass border-t border-slate-200/70 dark:border-slate-700/40"
            variants={mobileMenuVariants}
            initial="hidden"
            animate="visible"
            exit="exit"
          >
            <div className="px-2 pt-2 pb-3 space-y-1 sm:px-3">
              <Navigation isMobile={true} closeMenu={() => setIsMenuOpen(false)} />
            </div>
          </motion.div>
        )}
      </AnimatePresence>
    </motion.header>
  );
};

export default Header;
