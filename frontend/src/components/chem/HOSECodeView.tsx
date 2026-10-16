import { useState, type FormEvent } from "react";
import { FileSearch, Loader2 } from "lucide-react";
import SMILESInput from "../common/SMILESInput";
import { generateHOSECodes } from "../../services/chemService";
import { ToolSkeleton } from "@/components/feedback/ToolSkeleton";
import { GlassErrorCard } from "@/components/feedback/GlassErrorCard";
import { getErrorMessage } from "@/lib/error-messages";
import { Button } from "@/components/ui/button";
import { Cell, CopyButton, Seg, Select, Switch } from "@/components/common/porcelain-ui";
import { useRevealOnResult } from "@/hooks/useRevealOnResult";

type Toolkit = "rdkit" | "cdk";

/** HOSE codes: one spherical-environment code per atom of a structure. */
const HOSECodeView = () => {
  const [smiles, setSmiles] = useState("");
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [hoseCodes, setHoseCodes] = useState<string[]>([]);
  const [spheres, setSpheres] = useState(2);
  const [toolkit, setToolkit] = useState<Toolkit>("rdkit");
  const [ringsize, setRingsize] = useState(false);

  const handleSubmit = async (e: FormEvent) => {
    e.preventDefault();
    const trimmedSmiles = smiles.trim();
    if (!trimmedSmiles) {
      setError("Please enter a SMILES string.");
      setHoseCodes([]);
      return;
    }
    setLoading(true);
    setError(null);
    setHoseCodes([]);
    try {
      const result = await generateHOSECodes(trimmedSmiles, spheres, toolkit, ringsize);
      setHoseCodes(Array.isArray(result) ? result : []);
    } catch (err) {
      console.error("HOSE code generation error:", err);
      setError(getErrorMessage("chem", err));
      setHoseCodes([]);
    } finally {
      setLoading(false);
    }
  };

  const hasCodes = hoseCodes.length > 0 && !loading && !error;
  const title = hasCodes
    ? `${hoseCodes.length} HOSE code${hoseCodes.length === 1 ? "" : "s"}`
    : "Results";

  useRevealOnResult("hose-results", loading, hoseCodes.length > 0, !!error);

  return (
    <div className="pc-bento">
      <Cell id="hose-input" title="Structure input" className="pc-t-input">
        <form onSubmit={handleSubmit} className="pc-tool-form">
          <SMILESInput value={smiles} onChange={setSmiles} required />
          <div className="pc-field-row">
            <Select
              label="Spheres"
              value={String(spheres)}
              onChange={(v) => setSpheres(Number(v))}
              options={[
                ["1", "1"],
                ["2", "2"],
                ["3", "3"],
                ["4", "4"],
              ]}
            />
          </div>
          <Seg
            label="Toolkit"
            value={toolkit}
            onChange={setToolkit}
            options={[
              { value: "rdkit", label: "RDKit" },
              { value: "cdk", label: "CDK" },
            ]}
          />
          <Switch checked={ringsize} onChange={setRingsize} label="Include ring size" />
          <div className="pc-row pc-actions">
            <Button type="submit" disabled={!smiles.trim() || loading}>
              {loading ? <Loader2 className="animate-spin" /> : <FileSearch />}
              {loading ? "Generating…" : "Generate HOSE codes"}
            </Button>
          </div>
        </form>
      </Cell>

      <Cell id="hose-about" title="About HOSE codes" className="pc-t-side">
        <p className="pc-tool-note">
          HOSE (Hierarchically Ordered Spherical Environment) codes describe the chemical
          environment around each atom. They are used in NMR shift prediction, QSAR/QSPR studies and
          structure elucidation.
        </p>
        <dl className="pc-defs">
          <div>
            <dt>Spheres</dt>
            <dd>
              How many bonds out from the atom the code reaches. More spheres, more specific codes.
            </dd>
          </div>
        </dl>
      </Cell>

      <Cell
        id="hose-results"
        title={title}
        className="pc-t-wide"
        meta={hasCodes && <CopyButton text={hoseCodes.join("\n")} label="Copy all" />}
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
        {!loading && !error && !hoseCodes.length && (
          <p className="pc-empty">Enter a SMILES string to generate its HOSE codes.</p>
        )}
        {hasCodes && (
          <div className="max-h-96 overflow-y-auto rounded-xl border border-border">
            <table className="w-full border-collapse">
              <thead className="sticky top-0 bg-muted">
                <tr>
                  <th
                    scope="col"
                    className="px-4 py-2 text-left text-xs font-medium uppercase tracking-wider text-muted-foreground"
                  >
                    Atom index
                  </th>
                  <th
                    scope="col"
                    className="px-4 py-2 text-left text-xs font-medium uppercase tracking-wider text-muted-foreground"
                  >
                    HOSE code
                  </th>
                </tr>
              </thead>
              <tbody className="divide-y divide-border">
                {hoseCodes.map((code, index) => (
                  <tr key={index}>
                    <td className="px-4 py-2 text-sm whitespace-nowrap text-muted-foreground">
                      {index}
                    </td>
                    <td className="px-4 py-2 font-mono text-xs break-all text-foreground">
                      {code}
                    </td>
                  </tr>
                ))}
              </tbody>
            </table>
          </div>
        )}
      </Cell>
    </div>
  );
};

export default HOSECodeView;
