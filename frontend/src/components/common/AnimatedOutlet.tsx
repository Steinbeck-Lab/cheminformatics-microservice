/**
 * AnimatedOutlet -- wraps React Router's Outlet with page transition animations.
 *
 * Uses AnimatePresence mode="wait" with a fast opacity crossfade.
 */
import { cloneElement } from "react";
import { useOutlet, useLocation } from "react-router-dom";
import { AnimatePresence, motion } from "motion/react";

const pageTransition = {
  duration: 0.25,
  ease: [0.25, 0.1, 0.25, 1] as const,
};

export function AnimatedOutlet() {
  const location = useLocation();
  const element = useOutlet();

  // Key on the base route segment so tab-internal navigation
  // (e.g. /chem/stereoisomers → /chem/descriptors) does NOT
  // trigger a page transition — only true page switches do.
  const segments = location.pathname.split("/").filter(Boolean);
  const pageKey = segments[0] || "/";

  return (
    <AnimatePresence mode="wait">
      {element && (
        <motion.div
          key={pageKey}
          initial={{ opacity: 0 }}
          animate={{ opacity: 1 }}
          exit={{ opacity: 0 }}
          transition={pageTransition}
        >
          {cloneElement(element, { key: pageKey })}
        </motion.div>
      )}
    </AnimatePresence>
  );
}
