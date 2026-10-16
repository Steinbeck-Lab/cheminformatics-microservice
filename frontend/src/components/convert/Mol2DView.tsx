import { useState, type FormEvent } from "react";
import { Download, FileText, Loader2 } from "lucide-react";
import SMILESInput from "../common/SMILESInput";
import { generate2DCoordinates } from "../../services/convertService";
import { ToolSkeleton } from "@/components/feedback/ToolSkeleton";
import { GlassErrorCard } from "@/components/feedback/GlassErrorCard";
import { getErrorMessage } from "@/lib/error-messages";
import { Button } from "@/components/ui/button";
import { Cell, CopyButton, Seg } from "@/components/common/porcelain-ui";
import { useRevealOnResult } from "@/hooks/useRevealOnResult";

type Toolkit = "cdk" | "rdkit" | "openbabel";

const TOOLKIT_OPTIONS: { value: Toolkit; label: string }[] = [
  { value: "cdk", label: "CDK" },
  { value: "rdkit", label: "RDKit" },
  { value: "openbabel", label: "Open Babel" },
];

/** 2D coordinates: a SMILES string to an MDL molfile (V2000) with 2D coordinates. */
const Mol2DView = () => {
  const [smiles, setSmiles] = useState("");
  const [toolkit, setToolkit] = useState<Toolkit>("cdk");
  const [isLoading, setIsLoading] = useState(false);
  const [molblock, setMolblock] = useState<string | null>(null);
  // The toolkit that produced the shown molblock
  const [usedToolkit, setUsedToolkit] = useState<Toolkit>("cdk");
  const [error, setError] = useState<string | null>(null);

  const handleSubmit = async (e: FormEvent) => {
    e.preventDefault();
    const trimmedSmiles = smiles.trim();
    if (!trimmedSmiles) {
      setError("Please enter a SMILES string.");
      setMolblock(null);
      return;
    }

    setIsLoading(true);
    setError(null);
    setMolblock(null);

    try {
      const data = await generate2DCoordinates(trimmedSmiles, toolkit);
      if (typeof data === "string" && data.trim()) {
        setMolblock(data);
        setUsedToolkit(toolkit);
      } else {
        throw new Error("Received empty or invalid Molblock data from the server.");
      }
    } catch (err) {
      console.error("Molblock generation error:", err);
      setError(getErrorMessage("convert", err));
      setMolblock(null);
    } finally {
      setIsLoading(false);
    }
  };

  // Download the molblock as a .mol file named after the SMILES
  const downloadMolblock = () => {
    if (!molblock) return;
    try {
      const blob = new Blob([molblock], { type: "chemical/x-mdl-molfile;charset=utf-8" });
      const url = URL.createObjectURL(blob);
      const a = document.createElement("a");
      a.href = url;
      const filenameBase = smiles
        ? smiles.replace(/[^a-z0-9]/gi, "_").substring(0, 30)
        : "molecule";
      a.download = `${filenameBase}_2d.mol`;
      document.body.appendChild(a);
      a.click();
      document.body.removeChild(a);
      URL.revokeObjectURL(url);
    } catch (err) {
      console.error("Error creating download link:", err);
      setError("Could not create file for download.");
    }
  };

  const done = !!molblock && !isLoading;

  useRevealOnResult("m2d-results", isLoading, !!molblock, !!error);

  return (
    <div className="pc-bento">
      <Cell id="m2d-input" title="Structure input" className="pc-t-input">
        <form onSubmit={handleSubmit} className="pc-tool-form">
          <SMILESInput
            value={smiles}
            onChange={setSmiles}
            placeholder="Enter SMILES string (e.g., CN1C=NC2=C1C(=O)N(C(=O)N2C)C)"
            required
          />
          <div className="pc-row pc-actions">
            <Button type="submit" disabled={!smiles.trim() || isLoading}>
              {isLoading ? <Loader2 className="animate-spin" /> : <FileText />}
              {isLoading ? "Generating…" : "Generate 2D coordinates"}
            </Button>
          </div>
        </form>
      </Cell>

      <Cell id="m2d-opts" title="Options" className="pc-t-side">
        <div className="pc-tool-form">
          <div className="pc-opt-group">
            <span className="pc-lbl" id="m2d-toolkit">
              Toolkit
            </span>
            <Seg
              labelId="m2d-toolkit"
              full
              value={toolkit}
              onChange={setToolkit}
              options={TOOLKIT_OPTIONS}
            />
            <p className="pc-note">Different toolkits may produce different 2D layouts.</p>
          </div>
          <p className="pc-tool-note">
            Generates 2D coordinates and returns them as an MDL molfile (V2000).
          </p>
        </div>
      </Cell>

      <Cell
        id="m2d-results"
        title={done ? "2D molfile" : "Results"}
        className="pc-t-wide"
        meta={
          done &&
          molblock && (
            <div className="pc-row">
              <CopyButton text={molblock} label="Copy molblock" />
              <Button variant="outline" onClick={downloadMolblock}>
                <Download />
                Download (.mol)
              </Button>
            </div>
          )
        }
      >
        {isLoading && <ToolSkeleton variant="conversion" />}
        {error && !isLoading && (
          <GlassErrorCard
            message={error}
            onRetry={() => {
              setError(null);
              document.getElementById("smiles-input")?.focus();
            }}
          />
        )}
        {!molblock && !isLoading && !error && (
          <p className="pc-empty">Enter a SMILES string to generate its 2D coordinates.</p>
        )}
        {done && molblock && (
          <div className="grid gap-2">
            <pre className="max-h-96 overflow-auto rounded-xl border border-border bg-muted p-3 font-mono text-xs whitespace-pre text-foreground">
              {molblock}
            </pre>
            <p className="pc-note">
              Generated with{" "}
              {TOOLKIT_OPTIONS.find((opt) => opt.value === usedToolkit)?.label ?? usedToolkit}.
            </p>
          </div>
        )}
      </Cell>
    </div>
  );
};

export default Mol2DView;
