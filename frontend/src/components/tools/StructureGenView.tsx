import { useState, type FormEvent } from "react";
import { Check, CheckSquare, Download, FlaskConical, Loader2, XCircle } from "lucide-react";
import MoleculeCard from "../common/MoleculeCard";
import { ToolSkeleton } from "@/components/feedback/ToolSkeleton";
import { GlassErrorCard } from "@/components/feedback/GlassErrorCard";
import { getErrorMessage } from "@/lib/error-messages";
import toolsService from "../../services/toolsService";
import { Button } from "@/components/ui/button";
import { Cell } from "@/components/common/porcelain-ui";
import { useRevealOnResult } from "@/hooks/useRevealOnResult";

/** Metadata the generator returns next to the structures (newer API format). */
interface GenerationResult {
  structures: string[];
  total_count?: number;
  generated_count?: number;
  limit_applied?: boolean;
  settings?: Record<string, string>;
}

/** Structure generation: every constitutional isomer of a molecular formula (SURGE). */
const StructureGenView = () => {
  const [formula, setFormula] = useState("");
  const [isLoading, setIsLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [structures, setStructures] = useState<string[]>([]);
  const [selectedStructures, setSelectedStructures] = useState<string[]>([]);
  const [generationResult, setGenerationResult] = useState<GenerationResult | null>(null);

  const handleGenerateStructures = async (e: FormEvent) => {
    e.preventDefault();
    if (!formula) {
      setError("Please enter a molecular formula");
      return;
    }
    setIsLoading(true);
    setError(null);
    setStructures([]);
    setSelectedStructures([]);
    setGenerationResult(null);

    try {
      // The service is typed as string[], but the API has answered in several shapes.
      const result: unknown = await toolsService.generateStructures(formula);
      const output = (result as { output?: unknown } | null)?.output;
      const generation = output as GenerationResult | undefined;
      if (typeof result === "string") {
        setError(result);
        setStructures([]);
        setGenerationResult(null);
      } else if (generation?.structures && Array.isArray(generation.structures)) {
        setStructures(generation.structures);
        setGenerationResult(generation);
      } else if (Array.isArray(output)) {
        // Legacy format: just the structures array
        setStructures(output as string[]);
        setGenerationResult(null);
      } else if (Array.isArray(result)) {
        // Direct array format
        setStructures(result as string[]);
        setGenerationResult(null);
      } else {
        console.error("Unexpected API response format:", result);
        setError("Unexpected response format from API");
        setStructures([]);
        setGenerationResult(null);
      }
    } catch (err) {
      console.error("Error generating structures:", err);
      setError(getErrorMessage("tools", err));
      setStructures([]);
      setGenerationResult(null);
    } finally {
      setIsLoading(false);
    }
  };

  const toggleSelectStructure = (smiles: string) => {
    setSelectedStructures((prev) =>
      prev.includes(smiles) ? prev.filter((s) => s !== smiles) : [...prev, smiles]
    );
  };

  const downloadSelectedStructures = () => {
    if (selectedStructures.length === 0) return;
    const content = selectedStructures.join("\n");
    const blob = new Blob([content], { type: "text/plain" });
    const url = URL.createObjectURL(blob);
    const a = document.createElement("a");
    a.href = url;
    a.download = `${formula || "generated"}-structures.smi`;
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
    URL.revokeObjectURL(url);
  };

  const done = structures.length > 0 && !isLoading;
  const title = done
    ? `${structures.length} structure${structures.length === 1 ? "" : "s"}`
    : "Results";

  useRevealOnResult("sg-results", isLoading, structures.length > 0, !!generationResult, !!error);

  return (
    <div className="pc-bento">
      <Cell id="sg-input" title="Formula" className="pc-t-input">
        <form onSubmit={handleGenerateStructures} className="pc-tool-form">
          <div className="grid gap-2">
            <label htmlFor="formula-input" className="pc-lbl">
              Molecular formula
            </label>
            <input
              id="formula-input"
              type="text"
              className="pc-well"
              value={formula}
              onChange={(e) => setFormula(e.target.value)}
              placeholder="e.g., C6H6"
              spellCheck={false}
              autoComplete="off"
              required
            />
            <p className="pc-note">
              For example C6H6 or C8H10, with at most 10 heavy atoms. Results are limited to the
              first 1000 structures when the total exceeds this limit.
            </p>
          </div>
          <div className="pc-row pc-actions">
            <Button type="submit" disabled={!formula || isLoading}>
              {isLoading ? <Loader2 className="animate-spin" /> : <FlaskConical />}
              {isLoading ? "Generating…" : "Generate structures"}
            </Button>
          </div>
        </form>
      </Cell>

      <Cell id="sg-about" title="About structure generation" className="pc-t-side">
        <p className="pc-tool-note">
          Generates the chemical structures that match a molecular formula with the SURGE algorithm.
          Because of the computational cost, formulas are limited to 10 heavy atoms.
        </p>
        <dl className="pc-defs">
          <div>
            <dt>Reference</dt>
            <dd>
              McKay, B.D., Yirik, M.A. &amp; Steinbeck, C. Surge: a fast open-source chemical graph
              generator. J Cheminform 14, 24 (2022).
            </dd>
          </div>
        </dl>
      </Cell>

      <Cell
        id="sg-results"
        title={title}
        className="pc-t-wide"
        meta={done && formula ? <span className="pc-tag">{formula}</span> : undefined}
      >
        {isLoading && <ToolSkeleton variant="molecule" />}
        {error && !isLoading && (
          <GlassErrorCard
            message={error}
            onRetry={() => {
              setError(null);
              document.getElementById("formula-input")?.focus();
            }}
          />
        )}
        {!done && !isLoading && !error && (
          <p className="pc-empty">Enter a molecular formula to generate its structures.</p>
        )}
        {done && (
          <div className="pc-tool-form">
            {generationResult && (
              <>
                <dl className="pc-defs grid-cols-2 md:grid-cols-4">
                  <div>
                    <dt>Total possible</dt>
                    <dd>{generationResult.total_count?.toLocaleString()}</dd>
                  </div>
                  <div>
                    <dt>Generated</dt>
                    <dd>{generationResult.generated_count?.toLocaleString()}</dd>
                  </div>
                  <div>
                    <dt>Selected</dt>
                    <dd>{selectedStructures.length}</dd>
                  </div>
                  <div>
                    <dt>Results</dt>
                    <dd className={generationResult.limit_applied ? "text-saffron-ink" : ""}>
                      {generationResult.limit_applied ? "Limited" : "Complete"}
                    </dd>
                  </div>
                </dl>

                {generationResult.settings && (
                  <div className="grid gap-2">
                    <span className="pc-lbl">Generation settings</span>
                    <ul className="grid gap-1 text-sm">
                      {Object.entries(generationResult.settings).map(([flag, description]) => (
                        <li key={flag} className="flex flex-wrap items-center gap-2">
                          <code className="pc-tag">{flag}</code>
                          <span className="text-muted-foreground">{description}</span>
                        </li>
                      ))}
                    </ul>
                  </div>
                )}

                {generationResult.limit_applied && (
                  <p className="pc-note text-saffron-ink">
                    Results limited to first {generationResult.generated_count} structures out of{" "}
                    {generationResult.total_count?.toLocaleString()} total possible isomers.
                  </p>
                )}
              </>
            )}

            <div className="pc-row">
              <Button
                variant="outline"
                type="button"
                onClick={downloadSelectedStructures}
                disabled={selectedStructures.length === 0}
              >
                <Download />
                Download selected ({selectedStructures.length})
              </Button>
              <Button
                variant="outline"
                type="button"
                onClick={() => setSelectedStructures([...structures])}
              >
                <CheckSquare />
                Select all
              </Button>
              <Button
                variant="outline"
                type="button"
                onClick={() => setSelectedStructures([])}
                disabled={selectedStructures.length === 0}
              >
                <XCircle />
                Deselect all
              </Button>
            </div>
            <p className="pc-note">Select structures to download them as a SMILES file.</p>

            <div className="pc-result-grid">
              {structures.map((smiles, index) => {
                const selected = selectedStructures.includes(smiles);
                return (
                  <div
                    key={smiles + index}
                    className={`relative rounded-2xl ${selected ? "ring-2 ring-primary ring-offset-2 ring-offset-background" : ""}`}
                  >
                    <MoleculeCard
                      smiles={smiles}
                      title={`Isomer ${index + 1}`}
                      showActions={false}
                    />
                    {/* The whole card toggles selection, as before. */}
                    <button
                      type="button"
                      className="absolute inset-0 z-10 cursor-pointer rounded-2xl focus-visible:outline-2 focus-visible:outline-primary"
                      aria-pressed={selected}
                      aria-label={`Select Isomer ${index + 1}`}
                      onClick={() => toggleSelectStructure(smiles)}
                    >
                      {selected && (
                        <span className="absolute top-2 right-2 rounded-full bg-primary p-1">
                          <Check className="size-4 text-primary-foreground" aria-hidden="true" />
                        </span>
                      )}
                    </button>
                  </div>
                );
              })}
            </div>
          </div>
        )}
      </Cell>
    </div>
  );
};

export default StructureGenView;
