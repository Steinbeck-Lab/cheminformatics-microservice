/**
 * AddToCompareButton -- reusable "Add to Compare" button for tool result views.
 *
 * Self-contained button that any tool result view can drop in by providing
 * a SMILES string. Reads comparison state from context to show appropriate
 * state: ready to add, already added, or compare full.
 */
import { Plus, Check, GitCompareArrows } from "lucide-react";
import { Button } from "@/components/ui/button";
import { useComparison } from "@/context/ComparisonContext";
import { cn } from "@/lib/utils";

interface AddToCompareButtonProps {
  smiles: string;
  title?: string;
  imageUrl?: string;
  descriptors?: Record<string, string | number>;
  sourceToolId?: string;
  className?: string;
}

export function AddToCompareButton({
  smiles,
  title = "Molecule",
  imageUrl,
  descriptors,
  sourceToolId,
  className,
}: AddToCompareButtonProps) {
  const { addMolecule, canAdd, molecules } = useComparison();

  const isAlreadyAdded = molecules.some((m) => m.smiles === smiles);
  const isDisabled = !smiles || isAlreadyAdded || !canAdd;

  const handleClick = (e: React.MouseEvent) => {
    e.stopPropagation();
    if (!smiles || isAlreadyAdded || !canAdd) return;
    addMolecule({ smiles, title, imageUrl, descriptors, sourceToolId });
  };

  // Already added state
  if (isAlreadyAdded) {
    return (
      <Button
        variant="ghost"
        size="sm"
        disabled
        className={cn("gap-1.5 text-green-500", className)}
        title="Already in comparison"
      >
        <Check className="h-3.5 w-3.5" />
        <span className="hidden sm:inline text-xs">Added</span>
      </Button>
    );
  }

  // Compare full state
  if (!canAdd) {
    return (
      <Button
        variant="ghost"
        size="sm"
        disabled
        className={cn("gap-1.5", className)}
        title="Compare is full (2 molecules)"
      >
        <GitCompareArrows className="h-3.5 w-3.5" />
        <span className="hidden sm:inline text-xs">Full</span>
      </Button>
    );
  }

  // Default: ready to add
  return (
    <Button
      variant="outline"
      size="sm"
      onClick={handleClick}
      disabled={isDisabled}
      className={cn("gap-1.5 clay-interactive", className)}
      title="Add to comparison"
    >
      <Plus className="h-3.5 w-3.5" />
      <span className="hidden sm:inline text-xs">Compare</span>
    </Button>
  );
}
