import { useState, useEffect, type FormEvent } from "react";
import DOMPurify from "dompurify";
import { Check, CheckCircle, Clipboard, FlaskConical, Info, Loader2 } from "lucide-react";
import SMILESInput from "../common/SMILESInput";
import { fixRadicals } from "../../services/chemService";
import { generate2DDepictionEnhanced, generate2DDepiction } from "../../services/depictService";
import { ToolSkeleton } from "@/components/feedback/ToolSkeleton";
import { GlassErrorCard } from "@/components/feedback/GlassErrorCard";
import { getErrorMessage } from "@/lib/error-messages";
import { Button } from "@/components/ui/button";
import { Cell } from "@/components/common/porcelain-ui";
import { useRevealOnResult } from "@/hooks/useRevealOnResult";

/** What /chem/fixRadicals returns: counts, not the booleans the shared type declares. */
interface FixResult {
  fixed_smiles: string;
  radicals_detected: number;
  radicals_fixed: number;
}

const cleanSvg = (svg: string) => svg.replace(/ns0:/g, "").replace("xmlns:ns0=", "xmlns=");

function Depiction({ svg, label }: { svg: string | null; label: string }) {
  return svg ? (
    <div
      role="img"
      aria-label={label}
      className="flex h-80 w-full items-center justify-center rounded-xl border border-border bg-card p-4 [&>svg]:h-auto [&>svg]:max-h-full [&>svg]:w-auto [&>svg]:max-w-full"
      dangerouslySetInnerHTML={{
        __html: DOMPurify.sanitize(svg, { USE_PROFILES: { svg: true, svgFilters: true } }),
      }}
    />
  ) : (
    <div className="flex h-80 w-full items-center justify-center rounded-xl border border-border bg-muted">
      <p className="text-sm text-muted-foreground">Depiction unavailable</p>
    </div>
  );
}

/** Fix radicals: unpaired electrons on C, N and O are filled with implicit hydrogens (CDK). */
const FixRadicalsView = () => {
  const [smiles, setSmiles] = useState("");
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [result, setResult] = useState<FixResult | null>(null);
  const [copied, setCopied] = useState(false);
  const [originalDepiction, setOriginalDepiction] = useState<string | null>(null);
  const [fixedDepiction, setFixedDepiction] = useState<string | null>(null);
  const [depictionLoading, setDepictionLoading] = useState(false);

  // Draw both structures once there is a result.
  useEffect(() => {
    const generateDepictions = async () => {
      if (!result) return;
      setDepictionLoading(true);
      setOriginalDepiction(null);
      setFixedDepiction(null);
      try {
        // The original keeps its radicals visible (perceive_radicals).
        const originalSvg = await generate2DDepictionEnhanced(smiles, {
          width: 400,
          height: 400,
          perceive_radicals: true,
          annotate: "none",
          style: "cow",
        });
        setOriginalDepiction(cleanSvg(originalSvg));
        const fixedSvg = await generate2DDepiction(result.fixed_smiles, {
          width: 400,
          height: 400,
          toolkit: "cdk",
        });
        setFixedDepiction(cleanSvg(fixedSvg));
      } catch (err) {
        // The SMILES still show without a picture.
        console.error("Error generating depictions:", err);
      } finally {
        setDepictionLoading(false);
      }
    };
    generateDepictions();
  }, [result, smiles]);

  const handleCopySmiles = async () => {
    if (!result?.fixed_smiles) return;
    try {
      await navigator.clipboard.writeText(result.fixed_smiles);
      setCopied(true);
      setTimeout(() => setCopied(false), 2000);
    } catch (err) {
      console.error("Failed to copy SMILES:", err);
      try {
        const textArea = document.createElement("textarea");
        textArea.value = result.fixed_smiles;
        textArea.style.position = "fixed";
        textArea.style.left = "-999999px";
        document.body.appendChild(textArea);
        textArea.select();
        document.execCommand("copy");
        document.body.removeChild(textArea);
        setCopied(true);
        setTimeout(() => setCopied(false), 2000);
      } catch (fallbackErr) {
        console.error("Fallback copy also failed:", fallbackErr);
      }
    }
  };

  const handleSubmit = async (e: FormEvent) => {
    e.preventDefault();
    const trimmedSmiles = smiles.trim();
    if (!trimmedSmiles) {
      setError("Please enter a SMILES string");
      setResult(null);
      setOriginalDepiction(null);
      setFixedDepiction(null);
      return;
    }
    setLoading(true);
    setError(null);
    setResult(null);
    setOriginalDepiction(null);
    setFixedDepiction(null);
    try {
      const data = (await fixRadicals(trimmedSmiles)) as unknown as FixResult;
      setResult(data);
      // Not a failure: shown as a note above the result.
      if (data.radicals_detected === 0) {
        setError("No radicals detected in the input molecule.");
      }
    } catch (err) {
      console.error("Radical fixing error:", err);
      setError(getErrorMessage("chem", err));
      setResult(null);
    } finally {
      setLoading(false);
    }
  };

  const done = !!result && !loading;
  const rate =
    result && result.radicals_detected > 0
      ? Math.round((result.radicals_fixed / result.radicals_detected) * 100)
      : 100;
  const title = done
    ? `${result.radicals_fixed} radical${result.radicals_fixed === 1 ? "" : "s"} fixed`
    : "Results";

  useRevealOnResult("fr-results", loading, !!result, !!error);

  return (
    <div className="pc-bento">
      <Cell id="fr-input" title="Structure input" className="pc-t-input">
        <form onSubmit={handleSubmit} className="pc-tool-form">
          <SMILESInput
            value={smiles}
            onChange={setSmiles}
            required
            placeholder="Enter SMILES with radicals (e.g., [CH3], C[CH2], [OH])"
          />
          <div className="pc-row pc-actions">
            <Button type="submit" disabled={!smiles.trim() || loading}>
              {loading ? <Loader2 className="animate-spin" /> : <FlaskConical />}
              {loading ? "Fixing radicals…" : "Fix radicals"}
            </Button>
          </div>
        </form>
      </Cell>

      <Cell id="fr-about" title="About radical fixing" className="pc-t-side">
        <p className="pc-tool-note">
          Detects radical (unpaired) electrons on atoms and fixes them. Radicals on carbon, nitrogen
          and oxygen are supported.
        </p>
        <dl className="pc-defs">
          <div>
            <dt>Examples</dt>
            <dd>
              <span className="font-mono">[CH3]</span> methyl,{" "}
              <span className="font-mono">C[CH2]</span> ethyl,{" "}
              <span className="font-mono">[OH]</span> hydroxyl
            </dd>
          </div>
        </dl>
      </Cell>

      <Cell id="fr-results" title={title} className="pc-t-wide">
        {loading && <ToolSkeleton variant="molecule" />}
        {error && !loading && !result && (
          <GlassErrorCard
            message={error}
            onRetry={() => {
              setError(null);
              document.getElementById("smiles-input")?.focus();
            }}
          />
        )}
        {!result && !loading && !error && (
          <p className="pc-empty">Enter a SMILES string with radicals to fix them.</p>
        )}
        {done && (
          <div className="pc-tool-form">
            {error && <p className="pc-note">{error}</p>}
            <dl className="pc-defs sm:grid-cols-3">
              <div>
                <dt>Radicals detected</dt>
                <dd className="text-2xl font-bold text-foreground">{result.radicals_detected}</dd>
              </div>
              <div>
                <dt>Radicals fixed</dt>
                <dd className="text-2xl font-bold text-foreground">{result.radicals_fixed}</dd>
              </div>
              <div>
                <dt>Success rate</dt>
                <dd className="text-2xl font-bold text-foreground">{rate}%</dd>
              </div>
            </dl>
            {result.radicals_detected > result.radicals_fixed && (
              <p className="pc-row text-sm text-saffron-ink">
                <Info className="size-4 shrink-0" aria-hidden="true" />
                Some radicals were detected but not fixed. Only radicals on C, N and O are
                supported.
              </p>
            )}
          </div>
        )}
      </Cell>

      {done && (
        <Cell
          id="fr-compare"
          title="Structure comparison"
          className="pc-t-wide"
          meta={
            <Button
              variant="outline"
              onClick={handleCopySmiles}
              disabled={!result.fixed_smiles}
              title="Copy fixed SMILES to clipboard"
            >
              {copied ? <Check /> : <Clipboard />}
              {copied ? "Copied" : "Copy fixed SMILES"}
            </Button>
          }
        >
          {depictionLoading ? (
            <p className="pc-row justify-center py-8 text-muted-foreground">
              <Loader2 className="size-5 animate-spin" aria-hidden="true" />
              Generating depictions…
            </p>
          ) : (
            <div className="grid grid-cols-1 gap-6 md:grid-cols-2">
              <div className="pc-tool-form gap-3">
                <h3 className="text-sm font-semibold text-foreground">
                  Original structure (with radicals)
                </h3>
                <Depiction svg={originalDepiction} label="Original structure, with radicals" />
                <p className="rounded-xl bg-muted p-3 font-mono text-xs break-all text-foreground">
                  {smiles}
                </p>
              </div>
              <div className="pc-tool-form gap-3">
                <h3 className="text-sm font-semibold text-foreground">Fixed structure</h3>
                <Depiction svg={fixedDepiction} label="Fixed structure" />
                <p className="rounded-xl bg-muted p-3 font-mono text-xs break-all text-foreground">
                  {result.fixed_smiles}
                </p>
              </div>
            </div>
          )}
          {smiles === result.fixed_smiles && result.radicals_detected === 0 && (
            <p className="pc-note mt-4 text-center">
              The structures are identical: no radicals were detected in the input.
            </p>
          )}
          {smiles !== result.fixed_smiles && result.radicals_fixed > 0 && (
            <p className="pc-row mt-4 text-sm text-primary">
              <CheckCircle className="size-4 shrink-0" aria-hidden="true" />
              Fixed {result.radicals_fixed} radical{result.radicals_fixed !== 1 ? "s" : ""} by
              adding implicit hydrogens.
            </p>
          )}
        </Cell>
      )}
    </div>
  );
};

export default FixRadicalsView;
