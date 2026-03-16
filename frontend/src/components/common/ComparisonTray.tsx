/**
 * ComparisonTray -- fixed bottom slide-up tray for molecule comparison.
 *
 * Renders a persistent glass-styled bar at the bottom of the viewport.
 * Slides up when isOpen is true (first molecule added).
 * Shows molecule cards with remove buttons, "Compare Side-by-Side" action,
 * and "Clear All" link. Visible across all tool pages.
 */
import { motion, AnimatePresence } from "motion/react";
import { X, GitCompareArrows, Trash2 } from "lucide-react";
import { Button } from "@/components/ui/button";
import { useComparison } from "@/context/ComparisonContext";

const traySpring = { type: "spring" as const, stiffness: 300, damping: 25 };

export function ComparisonTray() {
  const { molecules, isOpen, canAdd, clearAll, removeMolecule, startCompare } = useComparison();

  return (
    <AnimatePresence>
      {isOpen && molecules.length > 0 && (
        <motion.div
          initial={{ y: "100%" }}
          animate={{ y: 0 }}
          exit={{ y: "100%" }}
          transition={traySpring}
          className="fixed bottom-0 left-0 right-0 z-50 glass-bold rounded-t-2xl"
        >
          <div className="max-w-4xl mx-auto px-4 py-4">
            {/* Header row */}
            <div className="flex items-center justify-between mb-3">
              <h3 className="text-sm font-semibold text-foreground">
                Compare Molecules ({molecules.length}/2)
              </h3>
              <button
                onClick={clearAll}
                className="text-xs text-muted-foreground hover:text-destructive transition-colors flex items-center gap-1"
              >
                <Trash2 className="h-3 w-3" />
                Clear All
              </button>
            </div>

            {/* Molecule cards row */}
            <div className="flex items-center gap-3">
              {molecules.map((mol) => (
                <div
                  key={mol.smiles}
                  className="flex items-center gap-2 glass-bold rounded-xl px-3 py-2 min-w-0 flex-1"
                >
                  <div className="min-w-0 flex-1">
                    <p className="text-sm font-medium text-foreground truncate">{mol.title}</p>
                    <p className="text-xs text-muted-foreground font-mono truncate">{mol.smiles}</p>
                  </div>
                  <button
                    onClick={() => removeMolecule(mol.smiles)}
                    className="shrink-0 p-1 rounded-full hover:bg-destructive/10 transition-colors"
                    aria-label={`Remove ${mol.title} from comparison`}
                  >
                    <X className="h-4 w-4 text-muted-foreground hover:text-destructive" />
                  </button>
                </div>
              ))}

              {/* Placeholder for second molecule slot */}
              {canAdd && (
                <div className="flex-1 rounded-xl border-2 border-dashed border-white/20 px-3 py-2 text-center">
                  <p className="text-xs text-muted-foreground">Add another molecule</p>
                </div>
              )}

              {/* Compare button */}
              <Button
                onClick={startCompare}
                disabled={molecules.length < 2}
                className="shrink-0 clay-interactive"
                size="sm"
              >
                <GitCompareArrows className="h-4 w-4" />
                <span className="hidden sm:inline">Compare Side-by-Side</span>
              </Button>
            </div>
          </div>
        </motion.div>
      )}
    </AnimatePresence>
  );
}
