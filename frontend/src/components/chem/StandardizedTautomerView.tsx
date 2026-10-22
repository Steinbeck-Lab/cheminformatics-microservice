import { useState, type FormEvent } from "react";
import { Loader2, RefreshCw } from "lucide-react";
import SMILESInput from "../common/SMILESInput";
import MoleculeCard from "../common/MoleculeCard";
import { generateStandardizedTautomer } from "../../services/chemService";
import { ToolSkeleton } from "@/components/feedback/ToolSkeleton";
import { GlassErrorCard } from "@/components/feedback/GlassErrorCard";
import { getErrorMessage } from "@/lib/error-messages";
import { Button } from "@/components/ui/button";
import { Cell, CopyButton } from "@/components/common/porcelain-ui";
import { useRevealOnResult } from "@/hooks/useRevealOnResult";

/** Standardized tautomer: the canonical tautomeric form of one structure (RDKit). */
const StandardizedTautomerView = () => {
  const [smiles, setSmiles] = useState("");
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [standardizedSmiles, setStandardizedSmiles] = useState<string | null>(null);

  const handleSubmit = async (e: FormEvent) => {
    e.preventDefault();
    const trimmedSmiles = smiles.trim();
    if (!trimmedSmiles) {
      setError("Please enter a SMILES string");
      setStandardizedSmiles(null);
      return;
    }
    setLoading(true);
    setError(null);
    setStandardizedSmiles(null);
    try {
      const result = await generateStandardizedTautomer(trimmedSmiles);
      setStandardizedSmiles(result);
      // Not a failure: the input is already canonical, so show it as the result with a note.
      if (!result || result === trimmedSmiles) {
        setError("Input is already the canonical tautomer or no different tautomers were found.");
        setStandardizedSmiles(trimmedSmiles);
      }
    } catch (err) {
      console.error("Tautomer standardization error:", err);
      setError(getErrorMessage("chem", err));
      setStandardizedSmiles(null);
    } finally {
      setLoading(false);
    }
  };

  const done = !!standardizedSmiles && !loading;

  useRevealOnResult("taut-results", loading, !!standardizedSmiles, !!error);

  return (
    <div className="pc-bento">
      <Cell id="taut-input" title="Structure input" className="pc-t-input">
        <form onSubmit={handleSubmit} className="pc-tool-form">
          <SMILESInput value={smiles} onChange={setSmiles} required />
          <div className="pc-row pc-actions">
            <Button type="submit" disabled={!smiles.trim() || loading}>
              {loading ? <Loader2 className="animate-spin" /> : <RefreshCw />}
              {loading ? "Standardizing…" : "Standardize tautomer"}
            </Button>
          </div>
        </form>
      </Cell>

      <Cell id="taut-about" title="About tautomers" className="pc-t-side">
        <p className="pc-tool-note">
          Tautomers are isomers that readily interconvert, usually by moving a proton and shifting a
          double bond. This tool returns the canonical form chosen by RDKit&apos;s
          TautomerEnumerator.
        </p>
        <dl className="pc-defs">
          <div>
            <dt>Keto form</dt>
            <dd className="font-mono">CC(=O)CC</dd>
          </div>
          <div>
            <dt>Enol form</dt>
            <dd className="font-mono">CC(O)=CC</dd>
          </div>
        </dl>
        <p className="pc-note">Both forms typically give the same standardized tautomer.</p>
      </Cell>

      <Cell
        id="taut-results"
        title={done ? "Standardized tautomer" : "Results"}
        className="pc-t-wide"
      >
        {loading && <ToolSkeleton variant="molecule" />}
        {error && !loading && !standardizedSmiles && (
          <GlassErrorCard
            message={error}
            onRetry={() => {
              setError(null);
              document.getElementById("smiles-input")?.focus();
            }}
          />
        )}
        {!loading && !error && !standardizedSmiles && (
          <p className="pc-empty">Enter a SMILES string to find its standardized tautomer.</p>
        )}
        {done && (
          <>
            {error && <p className="pc-note mb-3">{error}</p>}
            <div className="pc-result-grid">
              <MoleculeCard smiles={smiles} title="Original" size="md" />
              <MoleculeCard smiles={standardizedSmiles} title="Standardized" size="md" />
            </div>
            <dl className="pc-ids mt-4">
              <div className="pc-id">
                <dt>Original SMILES</dt>
                <dd>
                  {smiles}
                  <CopyButton icon text={smiles} label="Copy original SMILES" />
                </dd>
              </div>
              <div className="pc-id">
                <dt>Standardized SMILES</dt>
                <dd>
                  {standardizedSmiles}
                  <CopyButton icon text={standardizedSmiles} label="Copy standardized SMILES" />
                </dd>
              </div>
            </dl>
          </>
        )}
      </Cell>
    </div>
  );
};

export default StandardizedTautomerView;
