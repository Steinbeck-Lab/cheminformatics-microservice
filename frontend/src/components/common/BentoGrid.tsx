import React, { useRef, useEffect } from "react";
import { motion, AnimatePresence, LayoutGroup } from "motion/react";
import { X } from "lucide-react";
import { Button } from "@/components/ui/button";
import { useMediaQuery } from "@/hooks/useMediaQuery";
import { ResizableToolPanel } from "@/components/layout/ResizableToolPanel";
import type { BentoSize } from "@/config/bentoLayouts";

/** Spring transition shared by all bento layout animations. */
const bentoSpring = { type: "spring" as const, stiffness: 300, damping: 25 };

// --- BentoGrid (container) ---

interface BentoGridProps {
  /** Page key used for unique LayoutGroup id. */
  pageKey: string;
  /** Currently expanded card id, or null if none. */
  expandedId: string | null;
  children: React.ReactNode;
}

/**
 * Reusable bento grid container.
 *
 * Wraps children in a LayoutGroup and renders a responsive CSS grid.
 * When a card is expanded, switches to single-column layout.
 */
export function BentoGrid({ pageKey, expandedId, children }: BentoGridProps) {
  return (
    <LayoutGroup id={`bento-${pageKey}`}>
      <div
        className={`grid gap-4 ${
          expandedId ? "grid-cols-1" : "grid-cols-1 sm:grid-cols-2 lg:grid-cols-3 xl:grid-cols-4"
        }`}
      >
        {children}
      </div>
    </LayoutGroup>
  );
}

// --- BentoCard (individual card) ---

interface BentoCardProps {
  /** Unique tool identifier. */
  id: string;
  /** Page key for layoutId namespacing. */
  pageKey: string;
  /** Display name of the tool. */
  name: string;
  /** Short description of the tool. */
  description: string;
  /** Icon component from lucide-react. */
  icon: React.ComponentType<{ className?: string }>;
  /** Size category from bento layout config. */
  size: BentoSize;
  /** Whether this card is currently expanded. */
  isExpanded: boolean;
  /** Whether another card is currently expanded. */
  isOtherExpanded: boolean;
  /** Toggle expand/collapse. */
  onToggle: () => void;
  /** Tool component rendered inside expanded card (backward-compatible). */
  children?: React.ReactNode;
  /** Left/top panel content for resizable split layout (optional). */
  inputContent?: React.ReactNode;
  /** Right/bottom panel content for resizable split layout (optional). */
  outputContent?: React.ReactNode;
  /** Extra content rendered in expanded card header (e.g., AddToCompareButton). */
  headerExtra?: React.ReactNode;
}

/**
 * Get grid span classes based on bento size category.
 *
 * Hero = 2x2 on lg+, medium = 2-col span on sm+, compact = 1x1.
 * On mobile (<640px), all cards are col-span-1.
 */
function getSizeClasses(size: BentoSize, isMobile: boolean): string {
  if (isMobile) return "col-span-1";

  switch (size) {
    case "hero":
      return "col-span-1 sm:col-span-2 lg:col-span-2 lg:row-span-2";
    case "medium":
      return "col-span-1 sm:col-span-2";
    case "compact":
    default:
      return "col-span-1";
  }
}

/**
 * Reusable bento card with expand-in-place animation.
 *
 * Uses framer-motion layoutId for FLIP animation with spring physics.
 * When collapsed: shows icon, name, description.
 * When expanded: full-width with tool component and close button.
 * When another card is expanded: fades to 30% opacity, 95% scale.
 */
export function BentoCard({
  id,
  pageKey,
  name,
  description,
  icon: Icon,
  size,
  isExpanded,
  isOtherExpanded,
  onToggle,
  children,
  inputContent,
  outputContent,
  headerExtra,
}: BentoCardProps) {
  const cardRef = useRef<HTMLDivElement>(null);
  const isMobile = useMediaQuery("(max-width: 639px)");

  // Scroll expanded card into view
  useEffect(() => {
    if (isExpanded && cardRef.current) {
      // Small delay to let layout animation start
      const timer = setTimeout(() => {
        cardRef.current?.scrollIntoView({ behavior: "smooth", block: "nearest" });
      }, 100);
      return () => clearTimeout(timer);
    }
  }, [isExpanded]);

  const sizeClasses = isExpanded ? "col-span-full z-10" : getSizeClasses(size, isMobile);

  return (
    <motion.div
      ref={cardRef}
      layoutId={`bento-${pageKey}-${id}`}
      transition={bentoSpring}
      animate={{
        opacity: isOtherExpanded ? 0.3 : 1,
        scale: isOtherExpanded ? 0.95 : 1,
      }}
      className={`${sizeClasses} glass-bold rounded-2xl overflow-hidden ${
        isOtherExpanded ? "pointer-events-none" : "cursor-pointer"
      }`}
      onClick={isExpanded ? undefined : onToggle}
      role="button"
      tabIndex={isExpanded ? undefined : 0}
      onKeyDown={
        isExpanded
          ? undefined
          : (e: React.KeyboardEvent) => {
              if (e.key === "Enter" || e.key === " ") {
                e.preventDefault();
                onToggle();
              }
            }
      }
    >
      {/* Card Header */}
      <div className={`flex items-start gap-3 p-5 ${isExpanded ? "border-b border-white/10" : ""}`}>
        <Icon className="h-6 w-6 shrink-0 text-primary mt-0.5" />
        <div className="min-w-0 flex-1">
          <h3 className="text-base font-semibold text-foreground truncate">{name}</h3>
          <p className="text-sm text-muted-foreground mt-0.5 line-clamp-2">{description}</p>
        </div>
        {isExpanded && headerExtra}
        {isExpanded && (
          <Button
            variant="ghost"
            size="icon"
            className="shrink-0 -mt-1 -mr-2"
            onClick={(e: React.MouseEvent) => {
              e.stopPropagation();
              onToggle();
            }}
            aria-label="Close expanded tool"
          >
            <X className="h-5 w-5" />
          </Button>
        )}
      </div>

      {/* Expanded content */}
      <AnimatePresence>
        {isExpanded && (inputContent || outputContent || children) && (
          <motion.div
            initial={{ opacity: 0 }}
            animate={{ opacity: 1, transition: { delay: 0.2 } }}
            exit={{ opacity: 0 }}
            className="p-5 pt-0"
          >
            {inputContent && outputContent ? (
              <ResizableToolPanel
                toolId={id}
                inputContent={inputContent}
                outputContent={outputContent}
              />
            ) : (
              children
            )}
          </motion.div>
        )}
      </AnimatePresence>
    </motion.div>
  );
}
