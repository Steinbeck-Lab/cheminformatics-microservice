import React, { useState, useEffect } from "react";
import { Link, useLocation } from "react-router-dom";
import { useAppContext } from "../../context/AppContext";
import Navigation from "./Navigation";
import { motion } from "motion/react";
import { Menu, Moon, Sun } from "lucide-react";
import { Button } from "@/components/ui/button";
import { Sheet, SheetTrigger, SheetContent, SheetHeader, SheetTitle } from "@/components/ui/sheet";

// --- Animation Variants ---
const headerVariants = {
  hidden: { opacity: 0, y: -20 },
  visible: {
    opacity: 1,
    y: 0,
    transition: {
      duration: 0.7,
      ease: [0.22, 1, 0.36, 1],
      when: "beforeChildren",
      staggerChildren: 0.1,
    },
  },
};

const headerContentVariants = {
  hidden: { opacity: 0, y: -10, scale: 0.97 },
  visible: (i = 1) => ({
    opacity: 1,
    y: 0,
    scale: 1,
    transition: {
      duration: 0.5,
      delay: 0.05 * i,
      ease: [0.2, 0.65, 0.3, 0.9],
    },
  }),
};

// --- Day/Night Theme Toggle ---
const ThemeToggle = ({
  isDarkMode,
  toggleDarkMode,
}: {
  isDarkMode: boolean;
  toggleDarkMode: () => void;
}) => (
  <motion.button
    onClick={toggleDarkMode}
    className={`relative flex items-center w-[52px] h-[28px] rounded-full p-[3px] cursor-pointer overflow-hidden ring-1 ring-inset transition-colors duration-500 ${
      isDarkMode
        ? "bg-gradient-to-r from-indigo-950 via-slate-900 to-slate-800 ring-indigo-400/20"
        : "bg-gradient-to-r from-sky-300 via-sky-200 to-amber-200 ring-amber-300/30"
    }`}
    whileHover={{ scale: 1.08 }}
    whileTap={{ scale: 0.94 }}
    transition={{ type: "spring", stiffness: 400, damping: 20 }}
    aria-label={isDarkMode ? "Switch to light mode" : "Switch to dark mode"}
    role="switch"
    aria-checked={isDarkMode}
  >
    {/* Outer glow */}
    <motion.div
      className="absolute -inset-1 rounded-full blur-md -z-10"
      animate={{
        opacity: 0.5,
        background: isDarkMode
          ? "radial-gradient(circle, rgba(99,102,241,0.4), transparent)"
          : "radial-gradient(circle, rgba(250,204,21,0.35), transparent)",
      }}
      transition={{ duration: 0.6 }}
    />

    {/* Stars (dark mode only) */}
    {isDarkMode && (
      <>
        <motion.div
          className="absolute top-[5px] right-[10px] w-[3px] h-[3px] bg-white rounded-full"
          animate={{ opacity: [0.2, 0.9, 0.2], scale: [0.8, 1, 0.8] }}
          transition={{
            repeat: Infinity,
            duration: 2,
            ease: "easeInOut",
          }}
        />
        <motion.div
          className="absolute bottom-[5px] right-[16px] w-[2px] h-[2px] bg-white/70 rounded-full"
          animate={{ opacity: [0.4, 1, 0.4] }}
          transition={{
            repeat: Infinity,
            duration: 1.6,
            delay: 0.5,
            ease: "easeInOut",
          }}
        />
      </>
    )}

    {/* Sliding knob with icon */}
    <motion.div
      className="relative z-10 h-[22px] w-[22px] rounded-full flex items-center justify-center"
      animate={{
        x: isDarkMode ? 24 : 0,
        backgroundColor: isDarkMode ? "#1e1b4b" : "#ffffff",
        boxShadow: isDarkMode
          ? "0 1px 8px rgba(129,140,248,0.45)"
          : "0 1px 8px rgba(250,204,21,0.4)",
      }}
      transition={{ type: "spring", stiffness: 500, damping: 28 }}
    >
      <motion.div
        key={isDarkMode ? "moon" : "sun"}
        initial={{ rotate: isDarkMode ? -90 : 90, scale: 0.5, opacity: 0 }}
        animate={{ rotate: 0, scale: 1, opacity: 1 }}
        transition={{ type: "spring", stiffness: 300, damping: 15 }}
      >
        {isDarkMode ? (
          <Moon className="h-3 w-3 text-indigo-300" />
        ) : (
          <Sun className="h-3 w-3 text-amber-500" />
        )}
      </motion.div>
    </motion.div>
  </motion.button>
);

// --- Component ---
const Header = () => {
  const { isDarkMode, toggleDarkMode } = useAppContext();
  const [isSheetOpen, setIsSheetOpen] = useState(false);
  const [isScrolled, setIsScrolled] = useState(false);
  const location = useLocation();

  useEffect(() => {
    setIsSheetOpen(false);
  }, [location]);

  useEffect(() => {
    const handleScroll = () => setIsScrolled(window.scrollY > 10);
    window.addEventListener("scroll", handleScroll);
    return () => window.removeEventListener("scroll", handleScroll);
  }, []);

  const logoDark =
    "https://raw.githubusercontent.com/Steinbeck-Lab/cheminformatics-microservice/main/public/img/logo_small_inverted.png";
  const logoLight =
    "https://raw.githubusercontent.com/Steinbeck-Lab/cheminformatics-microservice/main/public/img/logo_small.png";

  // Floating pill — glass intensifies on scroll
  const pillClasses = isScrolled
    ? "backdrop-blur-xl bg-white/85 dark:bg-slate-900/90 shadow-lg shadow-slate-900/5 dark:shadow-black/25 border-slate-200/80 dark:border-slate-700/60"
    : "backdrop-blur-md bg-white/70 dark:bg-slate-900/70 shadow-md shadow-slate-900/[0.03] dark:shadow-black/10 border-slate-200/50 dark:border-slate-700/40";

  return (
    <motion.header
      className="fixed top-0 left-0 right-0 z-50 px-3 sm:px-4 pt-2"
      variants={headerVariants}
      initial="hidden"
      animate="visible"
    >
      <div
        className={`max-w-7xl mx-auto rounded-full border transition-all duration-300 ${pillClasses}`}
      >
        <div className="flex items-center h-14 px-4 sm:px-5 lg:px-6">
          {/* Logo and title */}
          <motion.div
            className="shrink-0 flex items-center"
            variants={headerContentVariants}
            custom={0}
          >
            <Link to="/" className="flex items-center group" aria-label="Homepage">
              <motion.div
                className="overflow-hidden rounded-md mr-2"
                whileHover={{
                  scale: 1.1,
                  rotate: [0, 3, -3, 0],
                }}
                whileTap={{ scale: 0.92 }}
                transition={{ type: "spring", stiffness: 400, damping: 17 }}
              >
                <img
                  src={isDarkMode ? logoDark : logoLight}
                  alt="Cheminformatics Microservice Logo"
                  className="h-8 w-auto"
                />
              </motion.div>

              <div className="flex flex-col justify-center">
                <span className="font-bold text-base sm:text-lg leading-tight bg-linear-to-r from-sky-600 to-indigo-600 dark:from-sky-300 dark:to-blue-400 text-transparent bg-clip-text">
                  Cheminformatics
                </span>
                <span className="hidden xl:inline text-[10px] leading-tight -mt-0.5 text-indigo-500 dark:text-sky-300">
                  Microservices
                </span>
              </div>
            </Link>
          </motion.div>

          {/* Separator — logo / nav */}
          <div className="hidden lg:block self-stretch my-3 w-px bg-slate-300 dark:bg-slate-600 mx-4 xl:mx-5 shrink-0" />

          {/* Desktop Navigation — fills center */}
          <motion.div
            className="hidden lg:flex items-center flex-1 min-w-0 justify-center"
            variants={headerContentVariants}
            custom={1}
          >
            <Navigation />
          </motion.div>

          {/* Separator — nav / toggle */}
          <div className="hidden lg:block self-stretch my-3 w-px bg-slate-300 dark:bg-slate-600 mx-4 xl:mx-5 shrink-0" />

          {/* Desktop Theme Toggle */}
          <motion.div
            className="hidden lg:flex items-center shrink-0"
            variants={headerContentVariants}
            custom={2}
          >
            <ThemeToggle isDarkMode={isDarkMode} toggleDarkMode={toggleDarkMode} />
          </motion.div>

          {/* Mobile menu button */}
          <div className="flex items-center lg:hidden ml-auto">
            <Sheet open={isSheetOpen} onOpenChange={setIsSheetOpen}>
              <SheetTrigger asChild>
                <motion.div whileHover={{ scale: 1.08 }} whileTap={{ scale: 0.92 }}>
                  <Button
                    variant="ghost"
                    size="icon"
                    className="rounded-full h-8 w-8 text-slate-600 hover:text-sky-600 hover:bg-slate-200/70 dark:text-slate-300 dark:hover:text-sky-300 dark:hover:bg-slate-700/70 focus:outline-hidden"
                    aria-label="Toggle mobile menu"
                  >
                    <Menu className="h-5 w-5" />
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
