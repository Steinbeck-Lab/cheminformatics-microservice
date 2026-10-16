import { useState, type FormEvent } from "react";
import { AlertCircle, Check, Info, Loader2, ShieldCheck } from "lucide-react";
import SMILESInput from "../common/SMILESInput";
import MoleculeCard from "../common/MoleculeCard";
import { checkStructureErrors } from "../../services/chemService";
import { ToolSkeleton } from "@/components/feedback/ToolSkeleton";
import { GlassErrorCard } from "@/components/feedback/GlassErrorCard";
import { getErrorMessage } from "@/lib/error-messages";
import { Button } from "@/components/ui/button";
import { Cell, CopyButton, Switch } from "@/components/common/porcelain-ui";
import { useRevealOnResult } from "@/hooks/useRevealOnResult";

/** What the /chem/errors endpoint returns (the shared StructureErrorResult type is incomplete). */
interface CheckResult {
  original?: { smi?: string; messages?: string[] };
  standardized?: { smi: string; messages?: string[] };
}

/** Validation messages: none, the all-clear, or a list of issues. */
function Messages({ messages }: { messages?: string[] }) {
  if (!messages || messages.length === 0) {
    return (
      <p className="pc-row text-sm text-muted-foreground">
        <Info className="size-4 shrink-0" aria-hidden="true" />
        No validation messages reported.
      </p>
    );
  }
  if (messages.length === 1 && messages[0] === "No Errors Found") {
    return (
      <p className="pc-row text-sm font-medium text-primary">
        <Check className="size-4 shrink-0" aria-hidden="true" />
        {messages[0]}
      </p>
    );
  }
  return (
    <div className="text-sm text-saffron-ink">
      <p className="pc-row font-medium">
        <AlertCircle className="size-4 shrink-0" aria-hidden="true" />
        Structure issues found
      </p>
      <ul className="mt-1 list-disc space-y-1 pl-6">
        {messages.map((msg, index) => (
          <li key={index}>{msg}</li>
        ))}
      </ul>
    </div>
  );
}

/** Structure check: ChEMBL curation-pipeline validation, with optional standardization. */
const StructureErrorView = () => {
  const [smiles, setSmiles] = useState("");
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [result, setResult] = useState<CheckResult | null>(null);
  const [fix, setFix] = useState(true);

  const handleSubmit = async (e: FormEvent) => {
    e.preventDefault();
    const trimmedSmiles = smiles.trim();
    if (!trimmedSmiles) {
      setError("Please enter a SMILES string.");
      setResult(null);
      return;
    }
    setLoading(true);
    setError(null);
    setResult(null);
    try {
      const validationResult = await checkStructureErrors(trimmedSmiles, fix);
      setResult(validationResult as unknown as CheckResult);
    } catch (err) {
      console.error("Structure check error:", err);
      setError(getErrorMessage("chem", err));
      setResult(null);
    } finally {
      setLoading(false);
    }
  };

  const done = !!result && !loading;
  const originalSmiles = result?.original?.smi || smiles;
  const standardized = fix && done ? result?.standardized : undefined;

  useRevealOnResult("se-results", loading, !!result, !!error);

  return (
    <div className="pc-bento">
      <Cell id="se-input" title="Structure input" className="pc-t-input">
        <form onSubmit={handleSubmit} className="pc-tool-form">
          <SMILESInput value={smiles} onChange={setSmiles} required />
          <Switch
            checked={fix}
            onChange={setFix}
            label="Fix issues by standardizing"
            hint="Uses the ChEMBL curation pipeline"
          />
          <div className="pc-row pc-actions">
            <Button type="submit" disabled={!smiles.trim() || loading}>
              {loading ? <Loader2 className="animate-spin" /> : <ShieldCheck />}
              {loading ? "Checking…" : "Check structure"}
            </Button>
          </div>
        </form>
      </Cell>

      <Cell id="se-about" title="About structure checks" className="pc-t-side">
        <p className="pc-tool-note">
          Checks a structure for issues such as incorrect valences or invalid aromaticity, using the
          ChEMBL structure curation pipeline. With fixing on, standardization rules (normalization,
          aromaticity perception and others) resolve common problems.
        </p>
      </Cell>

      <Cell id="se-results" title={done ? "Original structure" : "Results"} className="pc-t-wide">
        {loading && <ToolSkeleton variant="molecule" />}
        {error && !loading && (
          <GlassErrorCard
            message={error}
            onRetry={() => {
              setError(null);
              document.getElementById("smiles-input")?.focus();
            }}
          />
        )}
        {!loading && !error && !result && (
          <p className="pc-empty">Enter a SMILES string to check it for errors.</p>
        )}
        {done && (
          <div className="pc-tool-form">
            <dl className="pc-ids">
              <div className="pc-id">
                <dt>SMILES</dt>
                <dd>
                  {originalSmiles}
                  <CopyButton icon text={originalSmiles} label="Copy original SMILES" />
                </dd>
              </div>
            </dl>
            <Messages messages={result?.original?.messages} />
          </div>
        )}
      </Cell>

      {standardized && (
        <Cell id="se-fixed" title="Standardized structure" className="pc-t-wide">
          <div className="pc-tool-form">
            <dl className="pc-ids">
              <div className="pc-id">
                <dt>SMILES</dt>
                <dd>
                  {standardized.smi}
                  <CopyButton icon text={standardized.smi} label="Copy standardized SMILES" />
                </dd>
              </div>
            </dl>
            <Messages messages={standardized.messages} />
            <div className="pc-result-grid">
              <MoleculeCard smiles={originalSmiles} title="Original" size="sm" />
              <MoleculeCard smiles={standardized.smi} title="Standardized / Fixed" size="sm" />
            </div>
          </div>
        </Cell>
      )}
    </div>
  );
};

export default StructureErrorView;
