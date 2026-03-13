import { motion } from "motion/react";
import { GlassSkeleton } from "./GlassSkeleton";

export function RouteLoadingFallback() {
  return (
    <motion.div
      initial={{ opacity: 0 }}
      animate={{ opacity: 1 }}
      transition={{ duration: 0.3 }}
      data-testid="route-loading-fallback"
      className="w-full max-w-4xl mx-auto px-4 py-8 space-y-6"
    >
      <GlassSkeleton variant="card" />
      <GlassSkeleton variant="text" />
      <GlassSkeleton variant="text" className="w-3/4" />
    </motion.div>
  );
}
