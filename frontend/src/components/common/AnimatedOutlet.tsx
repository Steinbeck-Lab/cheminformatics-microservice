/**
 * AnimatedOutlet -- wraps React Router's Outlet with page transition animations.
 *
 * Uses AnimatePresence mode="wait": the old page fades out fast, the new one fades and settles in.
 */
import { cloneElement, Suspense } from "react";
import { useOutlet, useLocation } from "react-router-dom";
import { AnimatePresence, motion } from "motion/react";
import { RouteLoadingFallback } from "@/components/feedback/RouteLoadingFallback";

// exit faster than entrance; the new page settles with a confident deceleration
const enter = { duration: 0.32, ease: [0.16, 1, 0.3, 1] as const };
const leave = { duration: 0.14, ease: [0.4, 0, 1, 1] as const };

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
        <motion.div key={pageKey} exit={{ opacity: 0, transition: leave }}>
          {/* Suspense per page: the placeholder shows inside the transition, and the page fades in
              when its code arrives instead of popping in. */}
          <Suspense fallback={<RouteLoadingFallback />}>
            <motion.div
              initial={{ opacity: 0, y: 8 }}
              animate={{ opacity: 1, y: 0 }}
              transition={enter}
            >
              {cloneElement(element, { key: pageKey })}
            </motion.div>
          </Suspense>
        </motion.div>
      )}
    </AnimatePresence>
  );
}
