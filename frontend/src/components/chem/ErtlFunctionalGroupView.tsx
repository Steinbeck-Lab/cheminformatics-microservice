import { useState, type FormEvent } from "react";
import { Loader2, Search } from "lucide-react";
import SMILESInput from "../common/SMILESInput";
import HighlightedMoleculeCard from "../common/HighlightedMoleculeCard";
import { generateFunctionalGroups } from "../../services/chemService";
import { ToolSkeleton } from "@/components/feedback/ToolSkeleton";
import { GlassErrorCard } from "@/components/feedback/GlassErrorCard";
import { getErrorMessage } from "@/lib/error-messages";
import { Button } from "@/components/ui/button";
import { Cell } from "@/components/common/porcelain-ui";
import { useRevealOnResult } from "@/hooks/useRevealOnResult";

/** One group as the /chem/ertlfunctionalgroup endpoint may return it. */
type Group =
  | string
  | {
      None?: unknown;
      description?: string;
      type?: string;
      atoms?: string;
      atomIds?: number[];
      [key: string]: unknown;
    };

const isNone = (fg: Group) => typeof fg === "object" && !!fg?.None;

/** A readable label for one group, whatever shape the API sent. */
const formatFunctionalGroup = (fg: Group): string => {
  if (typeof fg === "string") return fg;
  if (!fg || typeof fg !== "object") return "Invalid group data";
  if (fg.None) return "No functional groups found";
  if (fg.description) return fg.description;
  if (fg.type && fg.atoms) return `${fg.type} (atoms: ${fg.atoms})`;
  if (fg.type && fg.atomIds) return `${fg.type} (atom IDs: ${fg.atomIds.join(", ")})`;
  if (fg.type) return fg.type;
  const props = Object.entries(fg)
    .filter(([key]) => key !== "None")
    .map(([key, value]) => `${key}: ${value}`)
    .join(", ");
  return props || JSON.stringify(fg);
};

/** Functional groups: Ertl's algorithm, with each group highlighted on the structure on click. */
const ErtlFunctionalGroupView = () => {
  const [smiles, setSmiles] = useState("");
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [functionalGroups, setFunctionalGroups] = useState<Group[]>([]);
  const [selectedGroupIndex, setSelectedGroupIndex] = useState<number | null>(null);

  const handleSubmit = async (e: FormEvent) => {
    e.preventDefault();
    const trimmedSmiles = smiles.trim();
    if (!trimmedSmiles) {
      setError("Please enter a SMILES string");
      setFunctionalGroups([]);
      return;
    }
    setLoading(true);
    setError(null);
    setFunctionalGroups([]);
    try {
      const result = await generateFunctionalGroups(trimmedSmiles);
      setFunctionalGroups(Array.isArray(result) ? (result as Group[]) : []);
    } catch (err) {
      console.error("Functional group detection error:", err);
      setError(getErrorMessage("chem", err));
      setFunctionalGroups([]);
    } finally {
      setLoading(false);
    }
  };

  const first = functionalGroups[0];
  const noGroupsFound =
    functionalGroups.length === 0 ||
    (functionalGroups.length === 1 && first !== undefined && isNone(first));
  const done = functionalGroups.length > 0 && !loading && !error;
  const selected = selectedGroupIndex !== null ? functionalGroups[selectedGroupIndex] : undefined;
  const title = !done
    ? "Results"
    : noGroupsFound
      ? "No functional groups found"
      : `${functionalGroups.length} functional group${functionalGroups.length !== 1 ? "s" : ""}`;

  useRevealOnResult("fg-results", loading, functionalGroups.length > 0, !!error);

  return (
    <div className="pc-bento">
      <Cell id="fg-input" title="Structure input" className="pc-t-input">
        <form onSubmit={handleSubmit} className="pc-tool-form">
          <SMILESInput value={smiles} onChange={setSmiles} required />
          <div className="pc-row pc-actions">
            <Button type="submit" disabled={!smiles.trim() || loading}>
              {loading ? <Loader2 className="animate-spin" /> : <Search />}
              {loading ? "Detecting…" : "Detect functional groups"}
            </Button>
          </div>
        </form>
      </Cell>

      <Cell id="fg-about" title="About functional groups" className="pc-t-side">
        <p className="pc-tool-note">
          Functional groups are arrangements of atoms that give a molecule its characteristic
          chemical properties. This tool finds them with the algorithm proposed by Peter Ertl.
        </p>
        <dl className="pc-defs">
          <div>
            <dt>Examples</dt>
            <dd>Hydroxyl (-OH), carbonyl (C=O), carboxyl (-COOH), amino (-NH₂) and many others.</dd>
          </div>
        </dl>
        <p className="pc-note">Ertl, P. J. Cheminform. 9, 36 (2017).</p>
      </Cell>

      <Cell id="fg-results" title={title} className="pc-t-wide">
        {loading && <ToolSkeleton variant="general" />}
        {error && !loading && (
          <GlassErrorCard
            message={error}
            onRetry={() => {
              setError(null);
              document.getElementById("smiles-input")?.focus();
            }}
          />
        )}
        {!functionalGroups.length && !loading && !error && (
          <p className="pc-empty">Enter a SMILES string to detect its functional groups.</p>
        )}
        {done && (
          <div className="grid grid-cols-1 gap-6 md:grid-cols-2">
            <div className="pc-tool-form">
              <HighlightedMoleculeCard
                smiles={smiles}
                title="Input Molecule"
                size="md"
                functionalGroups={functionalGroups}
                highlightedGroupIndex={selectedGroupIndex}
                showActions={true}
              />
              {selected !== undefined && !isNone(selected) && (
                <div className="pc-row justify-center">
                  <Button variant="outline" onClick={() => setSelectedGroupIndex(null)}>
                    Clear highlighting
                  </Button>
                </div>
              )}
            </div>
            <div className="pc-tool-form content-start">
              {noGroupsFound ? (
                <p className="pc-empty">No functional groups detected.</p>
              ) : (
                <>
                  <p className="pc-note">Select a group to highlight it in the structure.</p>
                  <ul className="grid gap-2">
                    {functionalGroups.map((fg, index) => {
                      const on = selectedGroupIndex === index;
                      const obj = typeof fg === "object" ? fg : null;
                      return (
                        <li key={index}>
                          <button
                            type="button"
                            aria-pressed={on}
                            disabled={isNone(fg)}
                            onClick={() => setSelectedGroupIndex(on ? null : index)}
                            className={`w-full rounded-xl border p-2 text-left text-sm transition-colors ${
                              on
                                ? "border-primary bg-accent text-primary"
                                : "border-border bg-card text-foreground hover:bg-muted"
                            }`}
                          >
                            <span className="flex items-center justify-between gap-2">
                              <span>{formatFunctionalGroup(fg)}</span>
                              <span className="text-xs text-muted-foreground">
                                {obj?.atomIds && obj.atomIds.length > 0 && (
                                  <>
                                    {obj.atomIds.length} atom{obj.atomIds.length !== 1 ? "s" : ""}
                                  </>
                                )}
                                {on && <span className="ml-2 font-medium">Highlighted</span>}
                              </span>
                            </span>
                            {obj?.atoms && obj.type && !obj.None && (
                              <span className="mt-1 block text-xs text-muted-foreground">
                                Type: {obj.type} | Atoms: {obj.atoms}
                              </span>
                            )}
                          </button>
                        </li>
                      );
                    })}
                  </ul>
                </>
              )}
            </div>
          </div>
        )}
      </Cell>
    </div>
  );
};

export default ErtlFunctionalGroupView;
