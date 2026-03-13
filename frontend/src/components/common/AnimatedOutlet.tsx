/**
 * AnimatedOutlet -- wraps React Router's Outlet with page transition animations.
 *
 * Uses AnimatePresence mode="sync" for overlap transitions per user decision.
 * Exiting page uses absolute positioning to prevent vertical stacking during overlap.
 */
import React from "react";
import { useOutlet, useLocation } from "react-router-dom";
import { AnimatePresence, motion } from "motion/react";

const pageTransition = {
  type: "spring" as const,
  stiffness: 100,
  damping: 20,
};

export function AnimatedOutlet() {
  const location = useLocation();
  const element = useOutlet();

  return (
    <div style={{ position: "relative" }}>
      <AnimatePresence mode="sync">
        {element && (
          <motion.div
            key={location.pathname}
            initial={{ opacity: 0, y: 8 }}
            animate={{ opacity: 1, y: 0 }}
            exit={{
              opacity: 0,
              y: -8,
              position: "absolute" as never,
              width: "100%",
            }}
            transition={pageTransition}
            style={{ width: "100%" }}
          >
            {React.cloneElement(element, { key: location.pathname })}
          </motion.div>
        )}
      </AnimatePresence>
    </div>
  );
}
