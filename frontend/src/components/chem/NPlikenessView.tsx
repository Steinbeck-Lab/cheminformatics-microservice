import { useState, type FormEvent } from "react";
import { Calculator, Info, Loader2 } from "lucide-react";
import SMILESInput from "../common/SMILESInput";
import { calculateNPLikeness } from "../../services/chemService";
import { ToolSkeleton } from "@/components/feedback/ToolSkeleton";
import { GlassErrorCard } from "@/components/feedback/GlassErrorCard";
import { getErrorMessage } from "@/lib/error-messages";
import { Button } from "@/components/ui/button";
import { Cell } from "@/components/common/porcelain-ui";
import { useRevealOnResult } from "@/hooks/useRevealOnResult";

const getNpScoreColor = (score: number | null) => {
  if (score === null || isNaN(score)) return "text-muted-foreground";
  if (score >= 0) return "text-primary";
  if (score >= -2) return "text-saffron-ink";
  return "text-madder-ink";
};

const getNpScoreInterpretation = (score: number | null) => {
  if (score === null || isNaN(score)) return "Enter a SMILES string to calculate the score.";
  if (score >= 1) return "Very high natural product likeness - likely to be a natural product";
  if (score >= 0) return "High natural product likeness";
  if (score >= -1) return "Moderate natural product likeness";
  if (score >= -2) return "Low natural product likeness";
  return "Very low natural product likeness - likely to be synthetic";
};

/** NP-likeness: how close a structure is to known natural products, on a −5 to +5 scale. */
const NPlikenessView = () => {
  const [smiles, setSmiles] = useState("");
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [npScore, setNpScore] = useState<number | null>(null);
  const [showInfo, setShowInfo] = useState(false);

  const handleSubmit = async (e: FormEvent) => {
    e.preventDefault();
    const trimmedSmiles = smiles.trim();
    if (!trimmedSmiles) {
      setError("Please enter a SMILES string.");
      setNpScore(null);
      return;
    }
    setLoading(true);
    setError(null);
    setNpScore(null);
    try {
      const score = await calculateNPLikeness(trimmedSmiles);
      if (typeof score === "number" && !isNaN(score)) {
        setNpScore(score);
      } else {
        console.warn("Received non-numeric score:", score);
        throw new Error("Received an invalid score format from the server.");
      }
    } catch (err) {
      console.error("NP-likeness calculation error:", err);
      setError(getErrorMessage("chem", err));
      setNpScore(null);
    } finally {
      setLoading(false);
    }
  };

  const hasScore = typeof npScore === "number" && !isNaN(npScore) && !loading && !error;

  useRevealOnResult("np-results", loading, npScore !== null, !!error);

  return (
    <div className="pc-bento">
      <Cell id="np-input" title="Structure input" className="pc-t-input">
        <form onSubmit={handleSubmit} className="pc-tool-form">
          <SMILESInput value={smiles} onChange={setSmiles} required />
          <div className="pc-row pc-actions">
            <Button type="submit" disabled={!smiles.trim() || loading}>
              {loading ? <Loader2 className="animate-spin" /> : <Calculator />}
              {loading ? "Calculating…" : "Calculate NP score"}
            </Button>
          </div>
        </form>
      </Cell>

      <Cell id="np-about" title="About NP-likeness" className="pc-t-side">
        <p className="pc-tool-note">
          Measures how similar a molecule is to known natural products. Scores run from −5
          (synthetic-like) to +5 (natural-product-like).
        </p>
        <p className="pc-note">
          RDKit implementation of the method by Ertl et al., J. Cheminform. 1:8 (2008).
        </p>
      </Cell>

      <Cell
        id="np-results"
        title={hasScore ? "NP-likeness score" : "Results"}
        className="pc-t-wide"
        meta={
          hasScore && (
            <Button
              variant="ghost"
              size="icon-sm"
              aria-label="Score ranges"
              aria-expanded={showInfo}
              onClick={() => setShowInfo(!showInfo)}
            >
              <Info />
            </Button>
          )
        }
      >
        {loading && <ToolSkeleton variant="descriptors" />}
        {error && !loading && (
          <GlassErrorCard
            message={error}
            onRetry={() => {
              setError(null);
              document.getElementById("smiles-input")?.focus();
            }}
          />
        )}
        {npScore === null && !loading && !error && (
          <p className="pc-empty">Enter a SMILES string to calculate its NP-likeness score.</p>
        )}
        {hasScore && npScore !== null && (
          <div className="pc-tool-form">
            {showInfo && (
              <dl className="pc-defs">
                <div>
                  <dt>Above 1</dt>
                  <dd>Very high NP-likeness</dd>
                </div>
                <div>
                  <dt>0 to 1</dt>
                  <dd>High NP-likeness</dd>
                </div>
                <div>
                  <dt>−1 to 0</dt>
                  <dd>Moderate NP-likeness</dd>
                </div>
                <div>
                  <dt>−2 to −1</dt>
                  <dd>Low NP-likeness</dd>
                </div>
                <div>
                  <dt>Below −2</dt>
                  <dd>Very low NP-likeness (synthetic)</dd>
                </div>
              </dl>
            )}
            <div className="text-center">
              <div className={`mb-2 text-5xl font-bold ${getNpScoreColor(npScore)}`}>
                {npScore.toFixed(2)}
              </div>
              <p className="text-sm text-muted-foreground">{getNpScoreInterpretation(npScore)}</p>
            </div>
            <div className="relative px-2 pt-6">
              <div className="mb-4 flex h-2.5 overflow-hidden rounded-full bg-accent">
                <div className="h-full w-[30%] bg-linear-to-r from-destructive to-saffron" />
                <div className="h-full w-[10%] bg-saffron" />
                <div className="h-full w-[10%] bg-linear-to-r from-saffron to-primary" />
                <div className="h-full w-[50%] bg-primary" />
              </div>
              <div
                className="absolute bottom-4 h-4 w-4 -translate-x-1/2 rounded-full bg-card ring-2 ring-border transition-all duration-500 ease-out"
                style={{ left: `${Math.max(0, Math.min(100, ((npScore + 5) / 10) * 100))}%` }}
                title={`Score: ${npScore.toFixed(2)}`}
              />
              <div className="mt-1 flex justify-between text-xs">
                <span className="text-madder-ink">−5</span>
                <span className="text-saffron-ink">−2</span>
                <span className="text-saffron-ink">−1</span>
                <span className="text-primary">0</span>
                <span className="text-primary">1</span>
                <span className="text-primary">5</span>
              </div>
            </div>
          </div>
        )}
      </Cell>
    </div>
  );
};

export default NPlikenessView;
