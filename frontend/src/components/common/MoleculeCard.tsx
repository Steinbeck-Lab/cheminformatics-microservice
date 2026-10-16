import {
  useState,
  useRef,
  useEffect,
  type KeyboardEvent,
  type MouseEvent,
  type SyntheticEvent,
} from "react";
import { useAppContext } from "../../context/AppContext";
import { AlertCircle, Check, ChevronDown, Copy, Download, Info, RefreshCw, X } from "lucide-react";
import { toast } from "sonner";
import { Button } from "@/components/ui/button";
import { Dialog, DialogContent, DialogTitle } from "@/components/ui/dialog";
import { cn } from "@/lib/utils";
import { errorText } from "@/lib/error-messages";
import { Input } from "@/components/ui/input";

// Fallback shown on the white stage when the depiction fails: "Error loading structure" in #4a4650
const FALLBACK_SVG_BASE64 =
  "data:image/svg+xml;base64,PHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHZpZXdCb3g9IjAgMCAxMDAgMTAwIj48dGV4dCB4PSI1MCIgeT0iNTAiIGRvbWluYW50LWJhc2VsaW5lPSJtaWRkbGUiIHRleHQtYW5jaG9yPSJtaWRkbGUiIGZvbnQtc2l6ZT0iMTAiIGZpbGw9IiM0YTQ2NTAiPkVycm9yIGxvYWRpbmcgc3RydWN0dXJlPC90ZXh0Pjwvc3ZnPg==";

type PropertyMap = Record<string, string | number | boolean | null | undefined>;

interface MoleculeDetails {
  smiles: string;
  physicalProperties?: PropertyMap;
  atomCounts?: PropertyMap;
  drugLikeProperties?: PropertyMap;
  molecularFormula?: string;
  inchi?: string;
  inchiKey?: string;
  estimatedProperties?: Record<string, number>;
  error?: string;
}

// Swap a broken depiction for the fallback once; a second failure is left alone.
const showFallback = (e: SyntheticEvent<HTMLImageElement>, pad = false) => {
  const img = e.currentTarget;
  if (img.src === FALLBACK_SVG_BASE64) return;
  img.src = FALLBACK_SVG_BASE64;
  if (pad) img.classList.add("p-4");
};

const panel = "rounded-xl border border-[var(--line)] bg-[var(--inset)] p-3";
const sectionTitle = "m-0 mb-2 text-sm font-semibold text-[var(--head)]";

const MoleculeCard = ({
  smiles,
  title = "Molecule",
  description,
  showActions = true,
  size = "xl", // sm, md, lg, xl
  onClick = null, // Optional click handler for the whole card
  isExpandable = false, // Enable expand/collapse for additional molecule details
}: {
  smiles: string;
  title?: string;
  description?: string;
  showActions?: boolean;
  size?: "sm" | "md" | "lg" | "xl";
  onClick?: ((smiles: string) => void) | null;
  isExpandable?: boolean;
}) => {
  const [copied, setCopied] = useState(false);
  const [showFullSmiles, setShowFullSmiles] = useState(false);
  const [showCopyModal, setShowCopyModal] = useState(false);
  const [showDetailsModal, setShowDetailsModal] = useState(false);
  const [moleculeDetails, setMoleculeDetails] = useState<MoleculeDetails | null>(null);
  const [isLoadingDetails, setIsLoadingDetails] = useState(false);
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
  const imageUrl = smiles
    ? `${baseUrl}/latest/depict/2D?smiles=${encodeURIComponent(smiles)}&width=512&height=512&toolkit=cdk`
    : FALLBACK_SVG_BASE64;

  // Copy with the Clipboard API, then execCommand, then a manual-copy dialog
  const handleCopy = (e: MouseEvent) => {
    e.stopPropagation(); // Prevent card click if clicking button

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
          .catch((err) => {
            console.debug("Clipboard API failed, trying alternative method", err);
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

  // Fallback copy method using execCommand
  const fallbackCopy = () => {
    try {
      const textArea = document.createElement("textarea");
      textArea.value = typeof smiles === "string" ? smiles : "";
      textArea.style.position = "fixed";
      textArea.style.left = "-999999px";
      textArea.style.top = "-999999px";
      document.body.appendChild(textArea);
      textArea.focus();
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

  // Auto-select text in copy modal when it appears
  useEffect(() => {
    if (showCopyModal) {
      const timer = setTimeout(() => copyTextRef.current?.select(), 100);
      return () => clearTimeout(timer);
    }
  }, [showCopyModal]);

  // Download the structure as an SDF (molfile) from the mol2D endpoint
  const handleDownload = (e: MouseEvent) => {
    e.stopPropagation(); // Prevent card click
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

  // Fetch molecule details from the descriptor endpoint
  const fetchMoleculeDetails = async (smiles: string) => {
    if (!smiles) return null;

    setIsLoadingDetails(true);

    try {
      const detailsUrl = `${baseUrl}/latest/chem/descriptors`;

      const params = new URLSearchParams();
      params.append("smiles", smiles);
      params.append("toolkit", "rdkit");
      params.append("format", "json");

      const response = await fetch(`${detailsUrl}?${params.toString()}`);

      if (!response.ok) {
        throw new Error(`Failed to fetch molecule details: ${response.statusText}`);
      }

      const data = await response.json();

      const processedData: MoleculeDetails = {
        smiles,
        physicalProperties: {
          molecularWeight: data.molecular_weight,
          exactMolecularWeight: data.exact_molecular_weight,
          alogp: data.alogp,
          rotatableBondCount: data.rotatable_bond_count,
          polarSurfaceArea: data.topological_polar_surface_area,
          formalCharge: data.formal_charge,
          fractionCsp3: data.fractioncsp3,
          vanDerWaalsVolume: data.van_der_waals_volume,
        },
        atomCounts: {
          atomCount: data.atom_count,
          heavyAtomCount: data.heavy_atom_count,
          aromaticRingsCount: data.aromatic_rings_count,
          minimumRingCount: data.number_of_minimal_rings,
        },
        drugLikeProperties: {
          hydrogenBondAcceptors: data.hydrogen_bond_acceptors,
          hydrogenBondDonors: data.hydrogen_bond_donors,
          hydrogenBondAcceptorsLipinski: data.hydrogen_bond_acceptors_lipinski,
          hydrogenBondDonorsLipinski: data.hydrogen_bond_donors_lipinski,
          lipinskiViolations: data.lipinski_rule_of_five_violations,
          qedDrugLikeliness: data.qed_drug_likeliness,
          npLikeness: data.nplikeness,
          hasLinearSugars: data.linear_sugars,
          hasCircularSugars: data.circular_sugars,
          murckoFramework: data.murcko_framework,
        },
        molecularFormula: data.molecular_formula,
      };

      setMoleculeDetails(processedData);
      return processedData;
    } catch (error) {
      console.error("Error fetching molecule descriptors:", error);

      // Rough estimates parsed from the SMILES string
      const fallbackDetails: MoleculeDetails = {
        smiles,
        estimatedProperties: {
          heavyAtomCount: (smiles.match(/[A-Z][a-z]?/g) || []).length,
          carbonCount: (smiles.match(/C/g) || []).length,
          oxygenCount: (smiles.match(/O/g) || []).length,
          nitrogenCount: (smiles.match(/N/g) || []).length,
          ringCount: (smiles.match(/[0-9]/g) || []).length / 2,
        },
        error: errorText(error),
      };

      setMoleculeDetails(fallbackDetails);
      return fallbackDetails;
    } finally {
      setIsLoadingDetails(false);
    }
  };

  const handleShowDetails = async (e: MouseEvent) => {
    e.stopPropagation(); // Prevent card click

    if (!smiles) {
      console.warn("No SMILES available to show details");
      return;
    }

    setShowDetailsModal(true);

    if (!moleculeDetails) {
      await fetchMoleculeDetails(smiles);
    }
  };

  // Click on the card itself (only wired when onClick is provided)
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

  // Format SMILES for display (shorten if too long)
  const MAX_SMILES_LENGTH = 35;
  const isLong = !!smiles && smiles.length > MAX_SMILES_LENGTH;
  const formattedSmiles =
    !showFullSmiles && isLong ? `${smiles.substring(0, MAX_SMILES_LENGTH)}...` : smiles;

  const renderPropertySection = (heading: string, properties: PropertyMap) => {
    if (!properties || Object.keys(properties).length === 0) return null;

    const formatValue = (key: string, value: PropertyMap[string]) => {
      if (value === null || value === undefined) return "N/A";

      if (typeof value === "number") {
        if (key.toLowerCase().includes("weight") || key.toLowerCase().includes("mass")) {
          return value.toFixed(2) + " g/mol";
        } else if (key.toLowerCase().includes("logp")) {
          return value.toFixed(2);
        } else if (key.toLowerCase().includes("surface") || key.toLowerCase().includes("area")) {
          return value.toFixed(2) + " Å²";
        } else if (Number.isInteger(value)) {
          return value.toString();
        } else {
          return value.toFixed(2);
        }
      }

      return value.toString();
    };

    const formatPropertyName = (key: string) => {
      const specialCases: Record<string, string> = {
        logP: "LogP",
        tpsa: "TPSA",
        hBondDonors: "H-Bond Donors",
        hBondAcceptors: "H-Bond Acceptors",
        molecularWeight: "Molecular Weight",
        polarSurfaceArea: "Polar Surface Area",
        rotatableBonds: "Rotatable Bonds",
        heavyAtomCount: "Heavy Atom Count",
        aromaticRingCount: "Aromatic Ring Count",
      };

      if (specialCases[key]) return specialCases[key];

      // camelCase to Title Case with spaces
      return key
        .replace(/([A-Z])/g, " $1")
        .replace(/^./, (str) => str.toUpperCase())
        .trim();
    };

    return (
      <div>
        <h4 className={sectionTitle}>{heading}</h4>
        <div className={panel}>
          {Object.entries(properties)
            .filter(([, value]) => value !== null)
            .map(([key, value]) => (
              <div
                key={key}
                className="flex justify-between gap-3 border-b border-[var(--line)] py-1.5 text-sm last:border-0"
              >
                <span className="text-muted-foreground">{formatPropertyName(key)}</span>
                <span className="font-mono text-foreground">{formatValue(key, value)}</span>
              </div>
            ))}
        </div>
      </div>
    );
  };

  const detailRow = (label: string, value: string | number) => (
    <div className="flex justify-between gap-3">
      <span className="text-muted-foreground">{label}</span>
      <span className="font-mono text-foreground">{value}</span>
    </div>
  );

  const phys = moleculeDetails?.physicalProperties;
  const counts = moleculeDetails?.atomCounts;
  const drug = moleculeDetails?.drugLikeProperties;
  const qed = drug?.qedDrugLikeliness;

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
          {isExpandable && (
            <Button
              variant="ghost"
              size="icon-sm"
              onClick={(e) => {
                e.stopPropagation();
                setExpanded((prev) => !prev);
                // Fetch details when expanding if not already loaded
                if (!expanded && !moleculeDetails) {
                  fetchMoleculeDetails(smiles);
                }
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

        {/* Depiction on the white stage */}
        <div className={cn("pc-stage p-2", sizeClasses[size] || "h-60")}>
          {smiles ? (
            <img
              key={imageUrl}
              src={imageUrl}
              alt={`Chemical structure of ${title}`}
              onError={(e) => {
                console.error(`Error loading image for ${title} from ${imageUrl}:`, e);
                showFallback(e, true);
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
            {isLoadingDetails && <p className="pc-note">Loading details...</p>}
            {!isLoadingDetails && moleculeDetails && (
              <>
                {moleculeDetails.molecularFormula &&
                  detailRow("Formula", moleculeDetails.molecularFormula)}
                {phys?.molecularWeight != null &&
                  detailRow("MW", `${Number(phys.molecularWeight).toFixed(2)} g/mol`)}
                {phys?.alogp != null && detailRow("ALogP", Number(phys.alogp).toFixed(2))}
                {counts?.heavyAtomCount != null &&
                  detailRow("Heavy atoms", String(counts.heavyAtomCount))}
                {drug?.hydrogenBondDonors != null &&
                  detailRow(
                    "HBD / HBA",
                    `${drug.hydrogenBondDonors} / ${drug.hydrogenBondAcceptors}`
                  )}
                <div className="pt-1">
                  <span className="text-xs text-muted-foreground">SMILES</span>
                  <div className="pc-mono mt-0.5 rounded-lg border border-[var(--line)] p-1.5 font-mono text-xs break-all text-foreground">
                    {smiles || "N/A"}
                  </div>
                </div>
              </>
            )}
            {!isLoadingDetails && moleculeDetails?.error && (
              <p className="pc-err on m-0">Could not load full details.</p>
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
            <Button
              variant="ghost"
              size="icon-sm"
              onClick={handleShowDetails}
              disabled={!smiles}
              title="View molecule details"
              aria-label="View molecule details"
            >
              <Info aria-hidden="true" />
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

      {/* Molecule details dialog */}
      <Dialog open={showDetailsModal} onOpenChange={setShowDetailsModal}>
        <DialogContent
          showCloseButton={false}
          className="block max-h-[90vh] w-full max-w-2xl gap-0 overflow-y-auto rounded-xl border-0 bg-card p-6"
        >
          <div className="mb-4 flex items-center justify-between">
            <DialogTitle className="text-lg font-semibold text-[var(--head)]">
              {title !== "Molecule" ? title : "Molecule Details"}
            </DialogTitle>
            <Button
              variant="ghost"
              size="icon-sm"
              onClick={() => setShowDetailsModal(false)}
              aria-label="Close"
            >
              <X aria-hidden="true" />
            </Button>
          </div>

          {isLoadingDetails && (
            <p className="pc-note py-8 text-center" role="status">
              Loading molecule information...
            </p>
          )}

          {!isLoadingDetails && moleculeDetails?.error && (
            <p className="pc-err on mb-4" role="alert">
              <AlertCircle className="pc-i size-4 shrink-0" aria-hidden="true" />
              <span>Could not fetch complete molecule details. Showing limited information.</span>
            </p>
          )}

          {!isLoadingDetails && moleculeDetails && (
            <div className="grid gap-4">
              <div className="pc-stage h-64 p-2">
                <img
                  src={imageUrl}
                  alt={`Chemical structure of ${title}`}
                  onError={(e) => showFallback(e)}
                />
              </div>

              <dl className="pc-ids">
                <div className="pc-id">
                  <dt>SMILES</dt>
                  <dd>{smiles || "N/A"}</dd>
                </div>
                {moleculeDetails.inchi && (
                  <div className="pc-id">
                    <dt>InChI</dt>
                    <dd>{moleculeDetails.inchi}</dd>
                  </div>
                )}
                {moleculeDetails.inchiKey && (
                  <div className="pc-id">
                    <dt>InChI Key</dt>
                    <dd>{moleculeDetails.inchiKey}</dd>
                  </div>
                )}
                {moleculeDetails.molecularFormula && (
                  <div className="pc-id">
                    <dt>Molecular Formula</dt>
                    <dd className="pc-formula">{moleculeDetails.molecularFormula}</dd>
                  </div>
                )}
              </dl>

              {renderPropertySection("Physical Properties", phys || {})}
              {renderPropertySection("Atom Count", counts || {})}
              {drug &&
                renderPropertySection("Drug-likeness Properties", {
                  hydrogenBondAcceptors: drug.hydrogenBondAcceptors,
                  hydrogenBondDonors: drug.hydrogenBondDonors,
                  lipinskiViolations: drug.lipinskiViolations,
                  qedDrugLikeliness: drug.qedDrugLikeliness,
                  npLikeness: drug.npLikeness,
                })}

              {typeof qed === "number" && (
                <div>
                  <h4 className={sectionTitle}>QED Score</h4>
                  <div className={panel}>
                    <p className="pc-note mb-2">
                      Quantitative Estimate of Drug-likeness (0-1, higher is better)
                    </p>
                    <div
                      className="h-3 w-full overflow-hidden rounded-full bg-[var(--inset-2)]"
                      role="meter"
                      aria-label="QED score"
                      aria-valuemin={0}
                      aria-valuemax={1}
                      aria-valuenow={qed}
                    >
                      <div
                        className="h-full rounded-full bg-primary"
                        style={{ width: `${(qed * 100).toFixed(0)}%` }}
                      />
                    </div>
                    <div className="mt-1 flex justify-between text-xs text-muted-foreground">
                      <span>0.0</span>
                      <span className="font-mono font-medium text-foreground">
                        {qed.toFixed(2)}
                      </span>
                      <span>1.0</span>
                    </div>
                  </div>
                </div>
              )}

              {moleculeDetails.error && (
                <div className={panel}>
                  <div className="mb-2 flex items-center gap-2">
                    <AlertCircle className="size-5 text-saffron-ink" aria-hidden="true" />
                    <h4 className="m-0 text-sm font-semibold text-[var(--head)]">
                      Estimated Information
                    </h4>
                  </div>
                  <p className="pc-note mb-2">
                    Unable to fetch complete molecule data. These values are rough estimates parsed
                    from the SMILES string:
                  </p>
                  <ul className="list-disc pl-5 text-sm text-foreground">
                    {Object.entries(moleculeDetails.estimatedProperties || {}).map(
                      ([key, value]) => (
                        <li key={key}>
                          <span className="font-medium">
                            {key.replace(/([A-Z])/g, " $1").trim()}:
                          </span>{" "}
                          {value}
                        </li>
                      )
                    )}
                  </ul>
                  <div className="mt-3 border-t border-[var(--line)] pt-3">
                    <Button
                      variant="outline"
                      size="sm"
                      onClick={() => fetchMoleculeDetails(smiles)}
                    >
                      <RefreshCw aria-hidden="true" />
                      Retry fetching details
                    </Button>
                  </div>
                </div>
              )}
            </div>
          )}

          <div className="mt-6 flex justify-end gap-3">
            <Button onClick={handleDownload}>
              <Download aria-hidden="true" />
              Download Structure
            </Button>
            <Button variant="outline" onClick={() => setShowDetailsModal(false)}>
              Close
            </Button>
          </div>
        </DialogContent>
      </Dialog>
    </>
  );
};
export default MoleculeCard;
