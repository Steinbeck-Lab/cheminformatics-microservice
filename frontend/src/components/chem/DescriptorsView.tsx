import { useState, type FormEvent } from "react";
import DOMPurify from "dompurify";
import { Calculator, Loader2 } from "lucide-react";
import SMILESInput from "../common/SMILESInput";
import { getDescriptors } from "../../services/chemService";
import { ToolSkeleton } from "@/components/feedback/ToolSkeleton";
import { GlassErrorCard } from "@/components/feedback/GlassErrorCard";
import { getErrorMessage } from "@/lib/error-messages";
import { Button } from "@/components/ui/button";
import { Cell, CopyButton, Seg } from "@/components/common/porcelain-ui";
import { useRevealOnResult } from "@/hooks/useRevealOnResult";

type Toolkit = "rdkit" | "cdk" | "all";
type Format = "json" | "html";

/** Descriptors: molecular descriptors from RDKit, CDK or both, as copyable JSON or an HTML table. */
const DescriptorsView = () => {
  const [smiles, setSmiles] = useState("");
  const [descriptors, setDescriptors] = useState<unknown>(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [toolkit, setToolkit] = useState<Toolkit>("rdkit");
  const [format, setFormat] = useState<Format>("json");
  // the settings the shown result was fetched with: switching the controls later must not
  // re-render an old result in the wrong format or label it with the wrong toolkit
  const [shownAs, setShownAs] = useState<{ format: Format; toolkit: Toolkit }>({
    format,
    toolkit,
  });

  const handleCalculateDescriptors = async (e: FormEvent) => {
    e.preventDefault();
    const trimmedSmiles = smiles.trim();
    if (!trimmedSmiles) {
      setError("Please enter a SMILES string.");
      return;
    }
    setLoading(true);
    setError(null);
    setDescriptors(null);
    try {
      const result = await getDescriptors(trimmedSmiles, toolkit, format);
      setShownAs({ format, toolkit });
      setDescriptors(result);
    } catch (err) {
      console.error("Descriptor calculation error:", err);
      setError(getErrorMessage("chem", err));
    } finally {
      setLoading(false);
    }
  };

  const isHtml = shownAs.format === "html" && typeof descriptors === "string";
  const entries =
    shownAs.format === "json" && typeof descriptors === "object" && descriptors !== null
      ? Object.entries(descriptors as Record<string, unknown>)
      : null;
  const json = entries ? JSON.stringify(descriptors, null, 2) : "";
  const done = !!descriptors && !loading;
  const title = !done
    ? "Results"
    : entries && entries.length > 0
      ? `${entries.length} descriptor${entries.length === 1 ? "" : "s"}`
      : "Descriptors";

  const renderDescriptors = () => {
    if (isHtml) {
      return (
        <div
          className="prose prose-sm max-w-none overflow-x-auto dark:prose-invert"
          dangerouslySetInnerHTML={{ __html: DOMPurify.sanitize(descriptors as string) }}
        />
      );
    }
    if (entries) {
      if (entries.length === 0) {
        return <p className="pc-empty">No descriptors returned.</p>;
      }
      return (
        <pre className="pc-well pc-mono m-0 max-h-[560px] overflow-auto text-sm whitespace-pre">
          {json}
        </pre>
      );
    }
    return (
      <p className="text-madder-ink">
        Could not render descriptors. Unexpected data format received.
      </p>
    );
  };

  useRevealOnResult("desc-results", loading, !!descriptors, !!error);

  return (
    <div className="pc-bento">
      <Cell id="desc-input" title="Structure input" className="pc-t-input">
        <form onSubmit={handleCalculateDescriptors} className="pc-tool-form">
          <SMILESInput value={smiles} onChange={setSmiles} required />
          <div className="pc-row pc-actions">
            <Button type="submit" disabled={!smiles.trim() || loading}>
              {loading ? <Loader2 className="animate-spin" /> : <Calculator />}
              {loading ? "Calculating…" : "Calculate descriptors"}
            </Button>
          </div>
        </form>
      </Cell>

      <Cell id="desc-opts" title="Options" className="pc-t-side">
        <div className="pc-tool-form">
          <div className="pc-opt-group">
            <span className="pc-lbl" id="desc-toolkit">
              Toolkit
            </span>
            <Seg
              labelId="desc-toolkit"
              full
              value={toolkit}
              onChange={setToolkit}
              options={[
                { value: "rdkit", label: "RDKit" },
                { value: "cdk", label: "CDK" },
                { value: "all", label: "Both" },
              ]}
            />
          </div>
          <div className="pc-opt-group">
            <span className="pc-lbl" id="desc-format">
              Output format
            </span>
            <Seg
              labelId="desc-format"
              full
              value={format}
              onChange={setFormat}
              options={[
                { value: "json", label: "JSON" },
                { value: "html", label: "HTML" },
              ]}
            />
          </div>
        </div>
      </Cell>

      <Cell
        id="desc-results"
        title={title}
        className="pc-t-wide"
        meta={
          done && (
            <>
              <span className="pc-tag">
                {shownAs.toolkit === "all" ? "RDKit + CDK" : shownAs.toolkit.toUpperCase()} ·{" "}
                {shownAs.format.toUpperCase()}
              </span>
              {json && <CopyButton text={json} label="Copy JSON" />}
            </>
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
        {!descriptors && !loading && !error && (
          <p className="pc-empty">Enter a SMILES string to calculate its molecular descriptors.</p>
        )}
        {done && renderDescriptors()}
      </Cell>
    </div>
  );
};

export default DescriptorsView;
