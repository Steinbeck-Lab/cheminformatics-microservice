import { useState, useRef, useMemo, useEffect, type KeyboardEvent, type MouseEvent } from "react";
import { useAppContext } from "../../context/AppContext";
import { Check, ChevronDown, Copy, Download, X } from "lucide-react";
import { toast } from "sonner";
import { Button } from "@/components/ui/button";
import { Dialog, DialogContent, DialogTitle } from "@/components/ui/dialog";
import { cn } from "@/lib/utils";
import { Input } from "@/components/ui/input";

// Fallback shown on the white stage when the depiction fails: "Error loading structure" in #4a4650
const FALLBACK_SVG_BASE64 =
  "data:image/svg+xml;base64,PHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHZpZXdCb3g9IjAgMCAxMDAgMTAwIj48dGV4dCB4PSI1MCIgeT0iNTAiIGRvbWluYW50LWJhc2VsaW5lPSJtaWRkbGUiIHRleHQtYW5jaG9yPSJtaWRkbGUiIGZvbnQtc2l6ZT0iMTAiIGZpbGw9IiM0YTQ2NTAiPkVycm9yIGxvYWRpbmcgc3RydWN0dXJlPC90ZXh0Pjwvc3ZnPg==";

// An Ertl functional group as the API returns it; a bare string carries no highlight data.
type FunctionalGroup =
  | string
  | {
      None?: unknown;
      description?: string;
      type?: string;
      atoms?: string;
      atomIds?: number[];
    };

const HighlightedMoleculeCard = ({
  smiles,
  title = "Molecule",
  description,
  functionalGroups = [],
  highlightedGroupIndex = null, // Index of the functional group to highlight
  showActions = true,
  size = "xl",
  onClick = null,
  isExpandable = false, // Enable expand/collapse for SMILES details
}: {
  smiles: string;
  title?: string;
  description?: string;
  functionalGroups?: FunctionalGroup[];
  highlightedGroupIndex?: number | null;
  showActions?: boolean;
  size?: "sm" | "md" | "lg" | "xl";
  onClick?: ((smiles: string) => void) | null;
  isExpandable?: boolean;
}) => {
  const [copied, setCopied] = useState(false);
  const [showFullSmiles, setShowFullSmiles] = useState(false);
  const [showCopyModal, setShowCopyModal] = useState(false);
  const [expanded, setExpanded] = useState(false);
  const copyTextRef = useRef<HTMLInputElement>(null);

  const { addRecentMolecule } = useAppContext();

  // Height of the depiction stage for each size
  const sizeClasses = {
    sm: "h-32",
    md: "h-48",
    lg: "h-64",
    xl: "h-80",
  };

  const baseUrl = (
    import.meta.env.VITE_API_URL || "https://dev.api.naturalproducts.net/latest"
  ).replace(/\/latest$/, "");

  const rawGroup =
    highlightedGroupIndex !== null ? functionalGroups[highlightedGroupIndex] : undefined;
  const selectedGroup = rawGroup && typeof rawGroup === "object" ? rawGroup : null;
  // A bare string group still counts as selected (it has no highlight data)
  const hasSelection = !!rawGroup && !selectedGroup?.None;

  // Generate highlighting information from functional groups
  const highlightInfo = useMemo((): { atomIds: number[] | null; pattern: string } => {
    const group = selectedGroup;
    if (!group || group.None) return { atomIds: null, pattern: "" };

    // Priority 1: Use the type field as SMARTS pattern - this contains the actual substructure pattern
    if (group.type) {
      const typeStr = group.type.toString().trim();
      if (typeStr) {
        // If we also have atomIds, return both for combined highlighting
        if (group.atomIds && Array.isArray(group.atomIds) && group.atomIds.length > 0) {
          return {
            atomIds: group.atomIds,
            pattern: typeStr,
          };
        }
        // Otherwise just use the pattern
        return { atomIds: null, pattern: typeStr };
      }
    }

    // Priority 2: Use atom indices alone if no type pattern available
    if (group.atomIds && Array.isArray(group.atomIds) && group.atomIds.length > 0) {
      return {
        atomIds: group.atomIds,
        pattern: "",
      };
    }

    // Priority 2: Generate SMARTS pattern based on atom type and functional group context
    if (group.atoms) {
      const atomStr = group.atoms.toLowerCase();
      let pattern = "";

      // Create more specific patterns based on the atom type and context
      switch (atomStr) {
        case "o":
          // Check if it's part of a carbonyl, hydroxyl, etc.
          if (group.type && group.type.toString().toLowerCase().includes("=o")) {
            pattern = "[C]=[O]"; // Carbonyl oxygen
          } else if (group.type && group.type.toString().toLowerCase().includes("oh")) {
            pattern = "[OH]"; // Hydroxyl oxygen
          } else {
            pattern = "[O]"; // Generic oxygen
          }
          break;
        case "n":
          // Check for different nitrogen contexts
          if (group.type && group.type.toString().toLowerCase().includes("nh2")) {
            pattern = "[NH2]"; // Primary amine
          } else if (group.type && group.type.toString().toLowerCase().includes("nh")) {
            pattern = "[NH]"; // Secondary amine
          } else {
            pattern = "[N]"; // Generic nitrogen
          }
          break;
        case "s":
          pattern = "[S]";
          break;
        case "p":
          pattern = "[P]";
          break;
        case "c":
          // For carbon, try to be more specific based on context
          if (group.type && group.type.toString().toLowerCase().includes("=o")) {
            pattern = "[C]=O"; // Carbonyl carbon
          } else if (group.type && group.type.toString().toLowerCase().includes("aromatic")) {
            pattern = "[c]"; // Aromatic carbon
          } else {
            pattern = "[C]"; // Aliphatic carbon
          }
          break;
        default:
          // Try to use the atom string directly as an element symbol
          pattern = `[${atomStr.toUpperCase()}]`;
      }

      return { atomIds: null, pattern };
    }

    // Priority 3: Store atom indices for future use when backend supports it
    if (group.atomIds && Array.isArray(group.atomIds) && group.atomIds.length > 0) {
      return {
        atomIds: group.atomIds,
        pattern: "", // Will fallback to basic pattern if atomIds not supported
      };
    }

    // Priority 4: Final fallback - try to extract meaningful patterns from description
    if (group.description) {
      const desc = group.description.toLowerCase();

      // Common functional group patterns
      const patterns: Record<string, string> = {
        hydroxyl: "[OH]",
        alcohol: "[OH]",
        carbonyl: "[C]=[O]",
        ketone: "[C](=[O])[C]",
        aldehyde: "[CHO]",
        carboxyl: "[C](=[O])[OH]",
        carboxylic: "[C](=[O])[OH]",
        amine: "[NH2,NH,N]",
        amino: "[NH2]",
        amide: "[C](=[O])[NH]",
        ester: "[C](=[O])[O]",
        ether: "[O]([C])[C]",
        phenyl: "c1ccccc1",
        benzene: "c1ccccc1",
        aromatic: "c",
        nitro: "[N+](=O)[O-]",
        sulfur: "[S]",
        phosphorus: "[P]",
      };

      for (const [name, pattern] of Object.entries(patterns)) {
        if (desc.includes(name)) {
          return { atomIds: null, pattern };
        }
      }
    }

    return { atomIds: null, pattern: "" };
  }, [selectedGroup]);

  // Generate image URL with highlighting
  const imageUrl = useMemo(() => {
    if (!smiles) return FALLBACK_SVG_BASE64;

    const params = new URLSearchParams({
      smiles: smiles,
      width: "512",
      height: "512",
      toolkit: "rdkit", // Use RDKit for better highlighting support
    });

    // Add highlighting based on available data (use anchor atoms + SMARTS pattern for precision)
    if (highlightInfo.atomIds && highlightInfo.atomIds.length > 0 && highlightInfo.pattern) {
      params.set("atomIds", highlightInfo.atomIds.join(","));
      params.set("highlight", highlightInfo.pattern);
    } else if (highlightInfo.atomIds && highlightInfo.atomIds.length > 0) {
      params.set("atomIds", highlightInfo.atomIds.join(","));
    } else if (highlightInfo.pattern) {
      params.set("highlight", highlightInfo.pattern);
    }

    return `${baseUrl}/latest/depict/2D?${params.toString()}`;
  }, [smiles, highlightInfo, baseUrl]);

  // Copy functionality
  const handleCopy = (e: MouseEvent) => {
    e.stopPropagation();

    const smilesToCopy = typeof smiles === "string" ? smiles : "";
    if (!smilesToCopy) {
      console.error("No SMILES to copy");
      return;
    }

    try {
      if (navigator.clipboard && navigator.clipboard.writeText) {
        navigator.clipboard
          .writeText(smilesToCopy)
          .then(() => {
            setCopied(true);
            setTimeout(() => setCopied(false), 2000);
          })
          .catch(() => {
            fallbackCopy();
          });
        return;
      }

      fallbackCopy();
    } catch (err) {
      console.error("Copy failed with all methods:", err);
      setShowCopyModal(true);
    }
  };

  const fallbackCopy = () => {
    try {
      const textArea = document.createElement("textarea");
      textArea.value = smiles;
      textArea.style.position = "fixed";
      textArea.style.left = "-999999px";
      textArea.style.top = "-999999px";
      document.body.appendChild(textArea);
      textArea.select();

      const successful = document.execCommand("copy");
      document.body.removeChild(textArea);

      if (successful) {
        setCopied(true);
        setTimeout(() => setCopied(false), 2000);
        return;
      } else {
        throw new Error("execCommand copy failed");
      }
    } catch (err) {
      console.debug("Fallback copy failed:", err);
      setShowCopyModal(true);
    }
  };

  // Auto-select text in copy modal
  useEffect(() => {
    if (showCopyModal) {
      const timer = setTimeout(() => copyTextRef.current?.select(), 100);
      return () => clearTimeout(timer);
    }
  }, [showCopyModal]);

  // Handle downloading the structure
  const handleDownload = (e: MouseEvent) => {
    e.stopPropagation();
    if (typeof smiles !== "string" || !smiles.trim()) {
      console.error("Cannot download: Invalid SMILES string provided.");
      return;
    }

    const molUrl = `${baseUrl}/latest/convert/mol2D?smiles=${encodeURIComponent(smiles)}&toolkit=cdk`;

    fetch(molUrl)
      .then((response) => {
        if (!response.ok) {
          throw new Error(`Failed to download structure: ${response.statusText}`);
        }
        return response.text();
      })
      .then((molData) => {
        const blob = new Blob([molData], { type: "chemical/x-mdl-molfile" });
        const url = URL.createObjectURL(blob);

        const a = document.createElement("a");
        a.href = url;
        const filename = `${title.replace(/[^a-z0-9]/gi, "_").toLowerCase() || "molecule"}.sdf`;
        a.download = filename;
        document.body.appendChild(a);
        a.click();
        document.body.removeChild(a);
        URL.revokeObjectURL(url);
      })
      .catch((err) => {
        console.error("Structure download failed:", err);
        toast.error("Failed to download structure. Please try again later.");
      });
  };

  // Handle card click
  const handleCardClick = () => {
    if (onClick) {
      onClick(smiles);
    }
    addRecentMolecule({
      smiles,
      name: title !== "Molecule" ? title : undefined,
      timestamp: new Date().toISOString(),
    });
  };

  // Only keys pressed on the card itself, not on the buttons inside it
  const handleCardKey = (e: KeyboardEvent<HTMLDivElement>) => {
    if (e.target !== e.currentTarget) return;
    if (e.key === "Enter" || e.key === " ") {
      e.preventDefault();
      handleCardClick();
    }
  };

  // Format SMILES for display
  const MAX_SMILES_LENGTH = 35;
  const isLong = !!smiles && smiles.length > MAX_SMILES_LENGTH;
  const formattedSmiles =
    !showFullSmiles && isLong ? `${smiles.substring(0, MAX_SMILES_LENGTH)}...` : smiles;

  return (
    <>
      <div
        className={cn(
          "grid content-start gap-3 rounded-2xl border border-[var(--line)] bg-[var(--inset)] p-3 transition-colors",
          onClick &&
            "cursor-pointer hover:border-primary/45 focus-visible:outline-2 focus-visible:outline-offset-2 focus-visible:outline-ring"
        )}
        onClick={onClick ? handleCardClick : undefined}
        onKeyDown={onClick ? handleCardKey : undefined}
        role={onClick ? "button" : "figure"}
        tabIndex={onClick ? 0 : undefined}
        aria-label={title}
      >
        {/* Title row */}
        <div className="flex min-h-8 items-center justify-between gap-2">
          <h3 className="m-0 truncate text-sm font-semibold text-[var(--head)]" title={title}>
            {title}
          </h3>
          <div className="flex items-center gap-2">
            {hasSelection && <span className="pc-tag">Highlighted</span>}
            {isExpandable && (
              <Button
                variant="ghost"
                size="icon-sm"
                onClick={(e) => {
                  e.stopPropagation();
                  setExpanded((prev) => !prev);
                }}
                aria-label={expanded ? "Collapse details" : "Expand details"}
                aria-expanded={expanded}
              >
                <ChevronDown
                  aria-hidden="true"
                  className={cn("transition-transform", expanded && "rotate-180")}
                />
              </Button>
            )}
          </div>
        </div>

        {/* Depiction on the white stage */}
        <div className={cn("pc-stage p-2", sizeClasses[size] || "h-60")}>
          {smiles ? (
            <img
              key={imageUrl}
              src={imageUrl}
              alt={`Chemical structure of ${title}`}
              onError={(e) => {
                console.error(`Error loading image for ${title} from ${imageUrl}:`, e);
                const img = e.currentTarget;
                if (img.src === FALLBACK_SVG_BASE64) return;
                img.src = FALLBACK_SVG_BASE64;
                img.classList.add("p-4");
              }}
            />
          ) : (
            <p className="pc-stage-empty">No structure to display.</p>
          )}
        </div>

        {/* SMILES line */}
        <div className="grid gap-1 text-sm">
          <div className="flex items-start gap-2">
            <span className="shrink-0 text-xs font-semibold text-muted-foreground">SMILES</span>
            {isLong ? (
              <button
                type="button"
                className="min-w-0 flex-1 cursor-pointer border-0 bg-transparent p-0 text-left text-foreground hover:text-primary"
                onClick={(e) => {
                  e.stopPropagation();
                  setShowFullSmiles(!showFullSmiles);
                }}
                aria-expanded={showFullSmiles}
                title={showFullSmiles ? "Click to shorten" : "Click to expand"}
              >
                <span className="pc-mono font-mono text-xs break-all">{formattedSmiles}</span>
                <span className="ml-1 text-xs font-medium text-primary">
                  {showFullSmiles ? "less" : "more"}
                </span>
              </button>
            ) : (
              <span className="pc-mono min-w-0 flex-1 font-mono text-xs break-all text-foreground">
                {formattedSmiles || "-"}
              </span>
            )}
          </div>
          {description && <p className="pc-note">{description}</p>}
        </div>

        {/* Expandable details */}
        {isExpandable && expanded && (
          <div className="grid gap-2 border-t border-[var(--line)] pt-3 text-sm">
            <div>
              <span className="text-xs text-muted-foreground">Full SMILES</span>
              <div className="pc-mono mt-0.5 rounded-lg border border-[var(--line)] p-1.5 font-mono text-xs break-all text-foreground">
                {smiles || "N/A"}
              </div>
            </div>
            {rawGroup && (
              <div>
                <span className="text-xs text-muted-foreground">Highlighted group</span>
                <div className="mt-0.5 text-xs text-foreground">
                  {selectedGroup?.description || selectedGroup?.type || "Unknown group"}
                </div>
              </div>
            )}
          </div>
        )}

        {/* Actions */}
        {showActions && (
          <div className="flex justify-end gap-1 border-t border-[var(--line)] pt-2">
            <Button
              variant="ghost"
              size="icon-sm"
              onClick={handleCopy}
              disabled={!smiles}
              title={copied ? "Copied!" : "Copy SMILES to clipboard"}
              aria-label={copied ? "SMILES copied" : "Copy SMILES"}
            >
              {copied ? <Check aria-hidden="true" /> : <Copy aria-hidden="true" />}
            </Button>
            <Button
              variant="ghost"
              size="icon-sm"
              onClick={handleDownload}
              disabled={!smiles}
              title="Download structure as SDF"
              aria-label="Download structure as SDF"
            >
              <Download aria-hidden="true" />
            </Button>
          </div>
        )}
      </div>

      {/* Manual copy dialog (shown when automatic copy fails) */}
      <Dialog open={showCopyModal} onOpenChange={setShowCopyModal}>
        <DialogContent
          showCloseButton={false}
          className="block w-full max-w-lg gap-0 rounded-xl border-0 bg-card p-6"
        >
          <div className="mb-4 flex items-center justify-between">
            <DialogTitle className="text-lg font-semibold text-[var(--head)]">
              Copy SMILES
            </DialogTitle>
            <Button
              variant="ghost"
              size="icon-sm"
              onClick={() => setShowCopyModal(false)}
              aria-label="Close"
            >
              <X aria-hidden="true" />
            </Button>
          </div>
          <p className="mb-4 text-foreground">
            Automatic copying failed. Please select and copy this text manually:
          </p>
          <div className="mb-4">
            <Input
              type="text"
              ref={copyTextRef}
              value={smiles || ""}
              readOnly
              className="w-full font-mono text-sm"
              onClick={(e) => (e.target as HTMLInputElement).select()}
            />
          </div>
          <div className="flex justify-end gap-3">
            <Button variant="outline" onClick={() => setShowCopyModal(false)}>
              Close
            </Button>
          </div>
        </DialogContent>
      </Dialog>
    </>
  );
};
export default HighlightedMoleculeCard;
