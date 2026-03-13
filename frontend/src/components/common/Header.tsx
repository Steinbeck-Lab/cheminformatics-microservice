import React, { useState, useEffect } from "react";
import { Link, useLocation } from "react-router-dom";
import { useAppContext } from "../../context/AppContext";
import Navigation from "./Navigation";
import { motion, LayoutGroup } from "motion/react";
import { Menu, Moon, Sun } from "lucide-react";
import { Button } from "@/components/ui/button";
import { Sheet, SheetTrigger, SheetContent, SheetHeader, SheetTitle } from "@/components/ui/sheet";

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

// Spring transition for the pill toggle
const pillSwitchTransition = {
  type: "spring",
  stiffness: 600,
  damping: 30,
};

// Reusable theme toggle component
const ThemeToggle = ({
  isDarkMode,
  toggleDarkMode,
}: {
  isDarkMode: boolean;
  toggleDarkMode: () => void;
}) => (
  <div
    className={`relative flex items-center w-[62px] h-8 rounded-full p-1 cursor-pointer transition-all duration-300 ease-in-out ${
      isDarkMode
        ? "bg-linear-to-r from-slate-800 to-slate-700 hover:shadow-inner hover:shadow-slate-900"
        : "bg-linear-to-r from-sky-100 to-indigo-100 hover:shadow-md hover:shadow-indigo-200/50"
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
            ? "bg-linear-to-br from-slate-800 to-slate-900"
            : "bg-linear-to-br from-white to-sky-50"
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
      <Sun
        className={`h-4 w-4 transition-colors ${isDarkMode ? "text-slate-500" : "text-yellow-500"}`}
      />
      <Moon
        className={`h-4 w-4 transition-colors ${isDarkMode ? "text-yellow-300" : "text-slate-400"}`}
      />
    </div>
  </div>
);

// --- Component ---
const Header = () => {
  const { isDarkMode, toggleDarkMode } = useAppContext();
  const [isSheetOpen, setIsSheetOpen] = useState(false);
  const [isScrolled, setIsScrolled] = useState(false);
  const location = useLocation();

  // Close sheet on route change
  useEffect(() => {
    setIsSheetOpen(false);
  }, [location]);

  // Check scroll position for glass effect intensity
  useEffect(() => {
    const handleScroll = () => {
      setIsScrolled(window.scrollY > 10);
    };

    window.addEventListener("scroll", handleScroll);
    return () => window.removeEventListener("scroll", handleScroll);
  }, []);

  // Use higher quality logo URLs
  const logoDark =
    "https://raw.githubusercontent.com/Steinbeck-Lab/cheminformatics-microservice/main/public/img/logo_small_inverted.png";
  const logoLight =
    "https://raw.githubusercontent.com/Steinbeck-Lab/cheminformatics-microservice/main/public/img/logo_small.png";

  // Calculate header background classes based on scroll state using dark: variant
  const headerBgClasses = isScrolled
    ? "backdrop-blur-lg bg-white dark:bg-slate-900/90 shadow-lg shadow-slate-200/50 dark:shadow-slate-900/30"
    : "backdrop-blur-sm bg-white dark:bg-slate-900/80 shadow-md shadow-slate-200/30 dark:shadow-slate-900/20";

  return (
    <motion.header
      className={`sticky top-0 z-50 border-b transition-all duration-300 ${headerBgClasses} border-slate-200/70 dark:border-slate-800/40`}
      variants={headerVariants}
      initial="hidden"
      animate="visible"
    >
      <div className="w-full mx-auto px-4 sm:px-6 lg:px-8">
        <div className="flex items-center h-16 justify-between md:justify-center">
          {/* Logo and title */}
          <motion.div
            className="shrink-0 flex items-center"
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
                  className="font-bold text-lg sm:text-xl leading-tight bg-linear-to-r from-sky-600 to-indigo-600 dark:from-sky-300 dark:to-blue-400 text-transparent bg-clip-text"
                  whileHover={{ scale: 1.03 }}
                >
                  Cheminformatics
                </motion.span>
                <motion.span
                  className="hidden md:inline text-[11px] leading-tight -mt-0.5 transition-colors text-indigo-500 group-hover:text-indigo-700 dark:text-sky-300 dark:group-hover:text-sky-200"
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

          {/* Desktop Theme Toggle */}
          <motion.div
            className="hidden md:flex items-center ml-6 lg:ml-10"
            variants={headerContentVariants}
            custom={2}
          >
            <ThemeToggle isDarkMode={isDarkMode} toggleDarkMode={toggleDarkMode} />
          </motion.div>

          {/* Mobile menu button - Sheet trigger */}
          <div className="flex items-center md:hidden">
            <Sheet open={isSheetOpen} onOpenChange={setIsSheetOpen}>
              <SheetTrigger asChild>
                <motion.div whileHover={{ scale: 1.05 }} whileTap={{ scale: 0.95 }}>
                  <Button
                    variant="ghost"
                    size="icon"
                    className="rounded-full text-slate-600 hover:text-sky-600 hover:bg-slate-200/70 dark:text-slate-300 dark:hover:text-sky-300 dark:hover:bg-slate-700/70 focus:outline-hidden focus:ring-2 focus:ring-offset-2 focus:ring-offset-transparent focus:ring-sky-500"
                    aria-label="Toggle mobile menu"
                  >
                    <Menu className="h-6 w-6" />
                  </Button>
                </motion.div>
              </SheetTrigger>
              <SheetContent
                side="right"
                showCloseButton={true}
                className="w-[280px] bg-white dark:bg-slate-900 border-l border-slate-200 dark:border-slate-800"
              >
                <SheetHeader>
                  <SheetTitle className="text-foreground">Navigation</SheetTitle>
                </SheetHeader>
                <div className="px-2 py-4 space-y-1">
                  <Navigation isMobile={true} closeMenu={() => setIsSheetOpen(false)} />

                  {/* Mobile theme toggle */}
                  <div className="mt-6 pt-4 border-t border-dashed border-slate-300/30 dark:border-slate-700/50 flex justify-center">
                    <div className="flex items-center gap-3">
                      <span className="text-sm text-muted-foreground">Theme:</span>
                      <ThemeToggle isDarkMode={isDarkMode} toggleDarkMode={toggleDarkMode} />
                    </div>
                  </div>
                </div>
              </SheetContent>
            </Sheet>
          </div>
        </div>
      </div>
    </motion.header>
  );
};
export default Header;
