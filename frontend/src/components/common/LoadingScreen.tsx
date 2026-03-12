// Description: A loading screen component that displays a molecule animation with a loading message.
import React from "react";
import { motion } from "framer-motion";

// --- Configuration (Moved outside component for clarity) ---

// Particle positions and properties for molecule animation
const particles = [
  // Center atom (larger, different color)
  {
    id: "c",
    x: 0,
    y: 0,
    size: 16,
    colorLight: "#4A90E2",
    colorDark: "#8ACDEA",
  }, // Blue tones
  // Surrounding atoms
  {
    id: "h1",
    x: -40,
    y: 0,
    size: 10,
    colorLight: "#F5A623",
    colorDark: "#FFD580",
  }, // Orange/Yellow tones
  {
    id: "h2",
    x: 40,
    y: 0,
    size: 10,
    colorLight: "#F5A623",
    colorDark: "#FFD580",
  }, // Orange/Yellow tones
  {
    id: "h3",
    x: 0,
    y: -40,
    size: 10,
    colorLight: "#7ED321",
    colorDark: "#B8E986",
  }, // Green tones
  {
    id: "h4",
    x: 0,
    y: 40,
    size: 10,
    colorLight: "#7ED321",
    colorDark: "#B8E986",
  }, // Green tones
];

// Bond connections (center to others)
const bonds = [
  { from: 0, to: 1 },
  { from: 0, to: 2 },
  { from: 0, to: 3 },
  { from: 0, to: 4 },
];

// --- Component ---

const LoadingScreen = ({ text = "Loading..." }) => {
  // Use a simple state or prop to determine if dark mode is active
  // For this example, we'll rely on Tailwind's 'dark:' prefix triggered by a parent class.
  // If you manage theme differently (e.g., context), adjust color selection accordingly.
  // NOTE: The document.documentElement check might not be reliable during server-side rendering or initial hydration.
  // A context-based or prop-based theme check is generally more robust in React.
  const isDarkMode =
    typeof window !== "undefined" && document.documentElement.classList.contains("dark");

  return (
    // Use AnimatePresence if this component might be conditionally rendered
    // and you want exit animations. Wrap the parent component usage if needed.
    <motion.div
      key="loading-screen" // Added key for potential AnimatePresence usage
      initial={{ opacity: 0 }}
      animate={{ opacity: 1 }}
      exit={{ opacity: 0 }}
      transition={{ duration: 0.3 }}
      // Overlay with backdrop blur, theme-aware background
      className="fixed inset-0 z-50 flex items-center justify-center backdrop-blur-md bg-gray-100/60 dark:bg-gray-900/70"
      aria-busy="true"
      aria-live="assertive" // Announce loading state changes
    >
      {/* Modal container with theme-aware background, border, shadow */}
      <motion.div
        initial={{ scale: 0.9, opacity: 0 }}
        animate={{ scale: 1, opacity: 1 }}
        transition={{ type: "spring", stiffness: 260, damping: 20, delay: 0.1 }}
        className="bg-gradient-to-br from-white via-gray-50 to-gray-100 dark:from-gray-800 dark:via-gray-800 dark:to-gray-900 p-8 rounded-2xl shadow-xl dark:shadow-2xl border border-gray-200 dark:border-gray-700/50 flex flex-col items-center w-64"
      >
        {/* Animation Container */}
        <motion.div
          className="w-32 h-32 mb-6 relative"
          animate={{ rotate: 360 }} // Rotate the whole molecule slowly
          transition={{ duration: 20, repeat: Infinity, ease: "linear" }}
        >
          {/* Bonds (SVG) */}
          <svg className="absolute inset-0 w-full h-full overflow-visible">
            <defs>{/* Define gradients for bonds if desired */}</defs>
            {bonds.map((bond, i) => {
              const p1 = particles[bond.from];
              const p2 = particles[bond.to];
              // Determine color based on theme
              const bondColor = isDarkMode ? p2.colorDark : p2.colorLight;

              return (
                <motion.line
                  key={`bond-${i}`}
                  // Centered coordinates within the 128x128 SVG canvas (w-32 h-32)
                  x1={p1.x + 64}
                  y1={p1.y + 64}
                  x2={p2.x + 64}
                  y2={p2.y + 64}
                  stroke={bondColor} // Use theme-aware color
                  strokeWidth="2.5" // Slightly thicker bond
                  strokeLinecap="round" // Rounded ends
                  initial={{ pathLength: 0, opacity: 0 }}
                  animate={{ pathLength: 1, opacity: 0.7 }}
                  transition={{
                    pathLength: {
                      duration: 1,
                      delay: 0.2 + i * 0.1,
                      ease: "easeInOut",
                    },
                    opacity: { duration: 0.5, delay: 0.2 + i * 0.1 },
                  }}
                />
              );
            })}
          </svg>

          {/* Atoms (Divs) */}
          {particles.map((particle, i) => {
            // Determine color based on theme
            const atomColor = isDarkMode ? particle.colorDark : particle.colorLight;
            return (
              <motion.div
                key={`atom-${particle.id}`}
                className="absolute rounded-full shadow-md" // Standard shadow
                style={{
                  width: particle.size,
                  height: particle.size,
                  backgroundColor: atomColor, // Use theme-aware color
                  // Position relative to the center (50%, 50%) of the container
                  left: `calc(50% + ${particle.x}px - ${particle.size / 2}px)`,
                  top: `calc(50% + ${particle.y}px - ${particle.size / 2}px)`,
                  // Add a subtle border that matches background for depth
                  border: `1px solid ${atomColor}50`, // Semi-transparent border
                }}
                // Individual atom animation (pulse/scale)
                animate={{
                  scale: [1, 1.15, 1], // Slightly larger pulse
                }}
                transition={{
                  duration: 1.8,
                  repeat: Infinity,
                  delay: i * 0.15, // Stagger animation slightly
                  ease: "easeInOut",
                }}
              />
            );
          })}
        </motion.div>

        {/* Loading Text */}
        <motion.div
          // Theme-aware text color
          className="text-lg font-medium text-gray-600 dark:text-gray-300"
          role="status" // Indicate this element describes the status
          aria-live="polite" // Announce changes politely
          animate={{
            opacity: [0.6, 1, 0.6], // Fade in/out effect
          }}
          transition={{
            duration: 1.8,
            repeat: Infinity,
            ease: "easeInOut",
          }}
        >
          {text} {/* Display the passed text prop */}
        </motion.div>
      </motion.div>
    </motion.div>
  );
};

export default LoadingScreen;
