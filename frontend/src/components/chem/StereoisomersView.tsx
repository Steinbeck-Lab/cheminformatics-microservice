import { useState, type FormEvent } from "react";
import { Box, Loader2 } from "lucide-react";
import SMILESInput from "../common/SMILESInput";
import MoleculeCard from "../common/MoleculeCard";
import { generateStereoisomers } from "../../services/chemService";
import { ToolSkeleton } from "@/components/feedback/ToolSkeleton";
import { GlassErrorCard } from "@/components/feedback/GlassErrorCard";
import { getErrorMessage } from "@/lib/error-messages";
import { Button } from "@/components/ui/button";
import { Cell } from "@/components/common/porcelain-ui";
import { useRevealOnResult } from "@/hooks/useRevealOnResult";

/** Stereoisomers: every enantiomer and diastereomer of one structure. */
const StereoisomersView = () => {
  const [smiles, setSmiles] = useState("");
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [isomers, setIsomers] = useState<string[] | null>(null);

  const submit = async (e: FormEvent) => {
    e.preventDefault();
    const input = smiles.trim();
    if (!input) return setError("Enter a SMILES string.");
    setLoading(true);
    setError(null);
    setIsomers(null);
    try {
      const result = await generateStereoisomers(input);
      setIsomers(Array.isArray(result) ? result : []);
    } catch (err) {
      setError(getErrorMessage("chem", err));
    } finally {
      setLoading(false);
    }
  };

  // The service returns the input alone when there is nothing to enumerate.
  const none =
    isomers !== null &&
    (isomers.length === 0 || (isomers.length === 1 && isomers[0] === smiles.trim()));
  const title =
    isomers && !none
      ? `${isomers.length} stereoisomer${isomers.length === 1 ? "" : "s"}`
      : "Results";

  useRevealOnResult("st-results", loading, !!isomers, !!error);

  return (
    <div className="pc-bento">
      <Cell id="st-input" title="Structure input" className="pc-t-input">
        <form onSubmit={submit} className="pc-tool-form">
          <SMILESInput value={smiles} onChange={setSmiles} required />
          <div className="pc-row pc-actions">
            <Button type="submit" disabled={!smiles.trim() || loading}>
              {loading ? <Loader2 className="animate-spin" /> : <Box />}
              {loading ? "Generating…" : "Generate stereoisomers"}
            </Button>
          </div>
        </form>
      </Cell>

      <Cell id="st-about" title="About stereoisomers" className="pc-t-side">
        <p className="pc-tool-note">
          Same formula, same connections, different arrangement in space.
        </p>
        <dl className="pc-defs">
          <div>
            <dt>Enantiomers</dt>
            <dd>Mirror images that cannot be superimposed.</dd>
          </div>
          <div>
            <dt>Diastereomers</dt>
            <dd>Stereoisomers that are not mirror images.</dd>
          </div>
        </dl>
      </Cell>

      <Cell id="st-results" title={title} className="pc-t-wide">
        {loading && <ToolSkeleton variant="molecule" />}
        {error && !loading && <GlassErrorCard message={error} onRetry={() => setError(null)} />}
        {!loading && !error && isomers === null && (
          <p className="pc-empty">Enter a SMILES string to list its stereoisomers.</p>
        )}
        {!loading && none && (
          <p className="pc-empty">This structure has no distinct stereoisomers.</p>
        )}
        {!loading && isomers && !none && (
          <div className="pc-result-grid">
            {isomers.map((s, i) => (
              <MoleculeCard key={s} smiles={s} title={`Stereoisomer ${i + 1}`} size="sm" />
            ))}
          </div>
        )}
      </Cell>
    </div>
  );
};

export default StereoisomersView;
