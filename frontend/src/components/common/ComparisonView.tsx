/**
 * ComparisonView -- side-by-side molecule comparison overlay.
 *
 * Renders a full-screen glass overlay when isComparing is true.
 * Shows two molecules in parallel columns with their structures,
 * properties, descriptors, and SMILES strings.
 */
import { motion, AnimatePresence } from "motion/react";
import { X, ArrowLeft, Copy, Check } from "lucide-react";
import { useState } from "react";
import { Button } from "@/components/ui/button";
import { useComparison } from "@/context/ComparisonContext";
import type { ComparisonMolecule } from "@/context/ComparisonContext";

/** Copy-to-clipboard button for SMILES strings. */
function CopySmiles({ smiles }: { smiles: string }) {
  const [copied, setCopied] = useState(false);

  const handleCopy = async () => {
    try {
      await navigator.clipboard.writeText(smiles);
      setCopied(true);
      setTimeout(() => setCopied(false), 2000);
    } catch {
      // Clipboard API not available in some contexts
    }
  };

  return (
    <button
      onClick={handleCopy}
      className="inline-flex items-center gap-1 text-xs text-muted-foreground hover:text-foreground transition-colors"
      title="Copy SMILES"
    >
      {copied ? <Check className="h-3 w-3 text-green-500" /> : <Copy className="h-3 w-3" />}
      {copied ? "Copied" : "Copy"}
    </button>
  );
}

/** Single molecule column in the comparison view. */
function MoleculeColumn({ mol }: { mol: ComparisonMolecule }) {
  return (
    <div className="glass-bold rounded-2xl p-5 flex flex-col gap-4">
      {/* Title */}
      <h3 className="text-lg font-semibold text-foreground">{mol.title}</h3>

      {/* Structure image or SMILES display */}
      {mol.imageUrl ? (
        <div className="flex justify-center">
          <img
            src={mol.imageUrl}
            alt={`Structure of ${mol.title}`}
            className="max-w-full max-h-48 object-contain rounded-lg"
          />
        </div>
      ) : (
        <div className="glass-bold rounded-xl p-3 text-center">
          <p className="text-sm text-muted-foreground font-mono break-all">{mol.smiles}</p>
        </div>
      )}

      {/* Descriptors table */}
      {mol.descriptors && Object.keys(mol.descriptors).length > 0 && (
        <div className="space-y-1">
          <h4 className="text-sm font-medium text-foreground">Properties</h4>
          <div className="rounded-xl overflow-hidden">
            {Object.entries(mol.descriptors).map(([key, value]) => (
              <div
                key={key}
                className="flex justify-between px-3 py-1.5 text-sm even:bg-white/5 odd:bg-transparent"
              >
                <span className="text-muted-foreground">{key}</span>
                <span className="text-foreground font-medium">{String(value)}</span>
              </div>
            ))}
          </div>
        </div>
      )}

      {/* SMILES string */}
      <div className="space-y-1">
        <div className="flex items-center justify-between">
          <h4 className="text-sm font-medium text-foreground">SMILES</h4>
          <CopySmiles smiles={mol.smiles} />
        </div>
        <p className="text-xs text-muted-foreground font-mono bg-white/5 rounded-lg p-2 break-all">
          {mol.smiles}
        </p>
      </div>

      {/* Source tool */}
      {mol.sourceToolId && (
        <p className="text-xs text-muted-foreground">
          Source: <span className="capitalize">{mol.sourceToolId}</span>
        </p>
      )}
    </div>
  );
}

/** Highlight differences between two molecules' descriptors. */
function DifferenceHighlight({
  molecules,
}: {
  molecules: [ComparisonMolecule, ComparisonMolecule];
}) {
  const [a, b] = molecules;
  if (!a.descriptors || !b.descriptors) return null;

  const allKeys = new Set([...Object.keys(a.descriptors), ...Object.keys(b.descriptors)]);
  const diffs: { key: string; valA: string; valB: string }[] = [];

  allKeys.forEach((key) => {
    const valA = a.descriptors?.[key];
    const valB = b.descriptors?.[key];
    if (valA !== valB) {
      diffs.push({
        key,
        valA: valA != null ? String(valA) : "-",
        valB: valB != null ? String(valB) : "-",
      });
    }
  });

  if (diffs.length === 0) return null;

  return (
    <div className="glass-bold rounded-2xl p-5">
      <h3 className="text-sm font-semibold text-foreground mb-3">Differences</h3>
      <div className="rounded-xl overflow-hidden">
        <div className="grid grid-cols-3 gap-2 px-3 py-1.5 text-xs font-medium text-muted-foreground border-b border-white/10">
          <span>Property</span>
          <span className="text-center">{a.title}</span>
          <span className="text-center">{b.title}</span>
        </div>
        {diffs.map(({ key, valA, valB }) => (
          <div key={key} className="grid grid-cols-3 gap-2 px-3 py-1.5 text-sm even:bg-white/5">
            <span className="text-muted-foreground">{key}</span>
            <span className="text-center text-foreground font-medium">{valA}</span>
            <span className="text-center text-foreground font-medium">{valB}</span>
          </div>
        ))}
      </div>
    </div>
  );
}

export function ComparisonView() {
  const { molecules, isComparing, stopCompare } = useComparison();

  return (
    <AnimatePresence>
      {isComparing && molecules.length === 2 && (
        <motion.div
          initial={{ opacity: 0 }}
          animate={{ opacity: 1 }}
          exit={{ opacity: 0 }}
          transition={{ duration: 0.3 }}
          className="fixed inset-0 z-[60] glass-bold overflow-y-auto"
        >
          <div className="max-w-6xl mx-auto px-4 py-8">
            {/* Header */}
            <div className="flex items-center justify-between mb-6">
              <h2 className="text-2xl font-bold text-foreground">Comparing Molecules</h2>
              <Button
                variant="ghost"
                size="icon"
                onClick={stopCompare}
                aria-label="Close comparison"
              >
                <X className="h-5 w-5" />
              </Button>
            </div>

            {/* Two-column layout */}
            <div className="grid grid-cols-1 md:grid-cols-2 gap-6 mb-6">
              {molecules.map((mol) => (
                <MoleculeColumn key={mol.smiles} mol={mol} />
              ))}
            </div>

            {/* Differences section */}
            <DifferenceHighlight
              molecules={molecules as [ComparisonMolecule, ComparisonMolecule]}
            />

            {/* Back to tray button */}
            <div className="mt-6 flex justify-center">
              <Button variant="outline" onClick={stopCompare} className="clay-interactive">
                <ArrowLeft className="h-4 w-4" />
                Back to Tray
              </Button>
            </div>
          </div>
        </motion.div>
      )}
    </AnimatePresence>
  );
}
