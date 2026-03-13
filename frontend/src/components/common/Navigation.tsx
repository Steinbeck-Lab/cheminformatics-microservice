import React from "react";
import { NavLink } from "react-router-dom";
import { motion, LayoutGroup } from "motion/react";
import {
  Camera,
  FlaskConical,
  House,
  LineChart,
  RefreshCw,
  SlidersHorizontal,
  Users,
} from "lucide-react";

// Navigation links — shortLabel used at lg, full label at xl+
const navLinks = [
  { to: "/home", label: "Home", shortLabel: "Home", icon: House, exact: true },
  { to: "/chem", label: "Chemical Analysis", shortLabel: "Analysis", icon: FlaskConical },
  { to: "/convert", label: "Format Conversion", shortLabel: "Convert", icon: RefreshCw },
  { to: "/depict", label: "Depiction", shortLabel: "Depict", icon: LineChart },
  { to: "/tools", label: "Tools", shortLabel: "Tools", icon: SlidersHorizontal },
  { to: "/ocsr", label: "OCSR", shortLabel: "OCSR", icon: Camera },
  { to: "/about", label: "About", shortLabel: "About", icon: Users },
];

// Smooth spring for icon hover/tap
const iconSpring = { type: "spring", stiffness: 400, damping: 17 };

// Icon animation variants — propagated from parent whileHover/whileTap
const iconVariants = {
  idle: { scale: 1, rotate: 0 },
  hover: { scale: 1.18, rotate: 8, transition: iconSpring },
  tap: { scale: 0.88, rotate: -4, transition: { duration: 0.1 } },
};

const Navigation = ({
  isMobile = false,
  closeMenu = () => {},
}: {
  isMobile?: boolean;
  closeMenu?: () => void;
}) => {
  // --- Mobile layout ---
  if (isMobile) {
    return (
      <nav className="flex flex-col space-y-1">
        {navLinks.map((link) => {
          const Icon = link.icon;
          return (
            <NavLink
              key={link.to}
              to={link.to}
              onClick={closeMenu}
              end={link.exact}
              className={({ isActive }) =>
                `flex items-center min-h-[44px] py-3 px-3 text-base font-medium rounded-lg transition-colors ${
                  isActive
                    ? "bg-slate-900 dark:bg-slate-100 text-white dark:text-slate-900 font-semibold border-l-4 border-primary"
                    : "text-slate-700 dark:text-slate-300 hover:bg-slate-100 dark:hover:bg-slate-800"
                }`
              }
            >
              <Icon className="shrink-0 h-5 w-5 mr-3" aria-hidden="true" />
              <span>{link.label}</span>
            </NavLink>
          );
        })}
      </nav>
    );
  }

  // --- Desktop layout with animated active pill ---
  return (
    <LayoutGroup id="desktopNav">
      <nav className="flex items-center gap-0.5 xl:gap-1">
        {navLinks.map((link) => {
          const Icon = link.icon;
          return (
            <NavLink key={link.to} to={link.to} end={link.exact} className="relative outline-none">
              {({ isActive }) => (
                <motion.div
                  className={`relative flex items-center px-2.5 xl:px-3 py-1.5 text-sm font-medium rounded-full cursor-pointer select-none ${
                    isActive
                      ? "text-white dark:text-slate-900"
                      : "text-slate-600 dark:text-slate-400"
                  }`}
                  whileHover="hover"
                  whileTap="tap"
                  initial="idle"
                  animate="idle"
                >
                  {/* Animated active pill background */}
                  {isActive && (
                    <motion.div
                      layoutId="activeNavPill"
                      className="absolute inset-0 bg-slate-900 dark:bg-slate-200 rounded-full shadow-sm"
                      transition={{
                        type: "spring",
                        stiffness: 350,
                        damping: 30,
                      }}
                    />
                  )}

                  {/* Hover highlight for inactive items */}
                  {!isActive && (
                    <motion.div
                      className="absolute inset-0 rounded-full bg-slate-200/0 dark:bg-slate-700/0"
                      variants={{
                        idle: {
                          backgroundColor: "rgba(0,0,0,0)",
                        },
                        hover: {
                          backgroundColor: "rgba(148,163,184,0.15)",
                        },
                      }}
                      transition={{ duration: 0.2 }}
                    />
                  )}

                  {/* Icon — hidden at lg, visible at xl */}
                  <motion.span
                    className="relative z-10 shrink-0 mr-1.5 hidden xl:inline-flex"
                    variants={iconVariants}
                  >
                    <Icon className="h-4 w-4" aria-hidden="true" />
                  </motion.span>

                  {/* Label — short at lg, full at xl */}
                  <span className="relative z-10 whitespace-nowrap">
                    <span className="xl:hidden">{link.shortLabel}</span>
                    <span className="hidden xl:inline">{link.label}</span>
                  </span>
                </motion.div>
              )}
            </NavLink>
          );
        })}
      </nav>
    </LayoutGroup>
  );
};
export default Navigation;
