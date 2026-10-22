import { useState, type FormEvent } from "react";
import { Calculator, Loader2 } from "lucide-react";
import SMILESInput from "../common/SMILESInput";
import MoleculeCard from "../common/MoleculeCard";
import { useAppContext } from "../../context/AppContext";
import api from "../../services/api";
import { ToolSkeleton } from "@/components/feedback/ToolSkeleton";
import { GlassErrorCard } from "@/components/feedback/GlassErrorCard";
import { getErrorMessage } from "@/lib/error-messages";
import { Button } from "@/components/ui/button";
import { Cell, Seg, Select } from "@/components/common/porcelain-ui";
import { useRevealOnResult } from "@/hooks/useRevealOnResult";

type Toolkit = "rdkit" | "cdk";

const getSimilarityColor = (value: number) => {
  if (isNaN(value)) return "text-muted-foreground";
  if (value >= 0.7) return "text-primary";
  if (value >= 0.4) return "text-saffron-ink";
  return "text-madder-ink";
};

const getGaugeBgColor = (value: number) => {
  if (isNaN(value)) return "bg-accent";
  if (value >= 0.7) return "bg-primary";
  if (value >= 0.4) return "bg-saffron";
  return "bg-destructive";
};

const getSimilarityDescription = (value: number) => {
  if (isNaN(value)) return "Enter SMILES and calculate similarity.";
  if (value >= 0.85) return "High similarity - very close structural match";
  if (value >= 0.7) return "Good similarity - significant structural overlap";
  if (value >= 0.4) return "Moderate similarity - some shared features";
  return "Low similarity - structurally different compounds";
};

/** Fingerprints each toolkit offers. */
const getFingerprinterOptions = (tk: Toolkit): [string, string][] =>
  tk === "rdkit"
    ? [
        ["ECFP", "ECFP (Extended Connectivity)"],
        ["RDKit", "RDKit Topological"],
        ["Atompairs", "Atom Pairs"],
        ["MACCS", "MACCS Keys"],
        ["MAPC", "MAPC (MinHashed Atom-Pair Chiral)"],
      ]
    : [
        ["ECFP", "ECFP (Extended Connectivity)"],
        ["PubChem", "PubChem Fingerprints"],
      ];

/** Tanimoto similarity of two structures from a chosen fingerprint. */
const TanimotoView = () => {
  const [smiles1, setSmiles1] = useState("");
  const [smiles2, setSmiles2] = useState("");
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [result, setResult] = useState<number | null>(null);
  const [fingerprinter, setFingerprinter] = useState("ECFP");
  const [toolkit, setToolkit] = useState<Toolkit>("rdkit");
  const [nBits, setNBits] = useState(2048);
  const [radius, setRadius] = useState(2);
  const { addRecentMolecule } = useAppContext();

  const calculateTanimotoSimilarity = async () => {
    const trimmedSmiles1 = smiles1.trim();
    const trimmedSmiles2 = smiles2.trim();
    if (!trimmedSmiles1 || !trimmedSmiles2) {
      throw new Error("Please provide both SMILES strings");
    }
    const response = await api.get(`/chem/tanimoto`, {
      params: {
        smiles: `${trimmedSmiles1},${trimmedSmiles2}`,
        toolkit,
        fingerprinter,
        nBits,
        radius,
      },
    });
    if (typeof response.data !== "number") {
      throw new Error("Received invalid similarity score from API.");
    }
    return response.data;
  };

  const handleCalculate = async (e?: FormEvent) => {
    e?.preventDefault();
    const trimmedSmiles1 = smiles1.trim();
    const trimmedSmiles2 = smiles2.trim();
    if (!trimmedSmiles1 || !trimmedSmiles2) {
      setError("Please provide both SMILES strings");
      setResult(null);
      return;
    }
    setLoading(true);
    setError(null);
    setResult(null);
    try {
      const similarity = await calculateTanimotoSimilarity();
      setResult(similarity);
      addRecentMolecule({
        smiles: trimmedSmiles1,
        name: "Molecule 1 (Tanimoto)",
        timestamp: new Date().toISOString(),
      });
      addRecentMolecule({
        smiles: trimmedSmiles2,
        name: "Molecule 2 (Tanimoto)",
        timestamp: new Date().toISOString(),
      });
    } catch (err) {
      console.error("Tanimoto calculation error:", err);
      setError(getErrorMessage("chem", err));
      setResult(null);
    } finally {
      setLoading(false);
    }
  };

  // A new toolkit may not offer the current fingerprint: fall back to its first one.
  const handleToolkitChange = (newToolkit: Toolkit) => {
    setToolkit(newToolkit);
    const newOptions = getFingerprinterOptions(newToolkit);
    if (!newOptions.some(([value]) => value === fingerprinter)) {
      setFingerprinter(newOptions[0]?.[0] ?? "ECFP");
    }
  };

  const usesRadius = fingerprinter === "ECFP" || fingerprinter === "MAPC";
  const hasResult = typeof result === "number" && !isNaN(result) && !loading;

  useRevealOnResult("tc-results", loading, result !== null, !!error);

  return (
    <div className="pc-bento">
      <Cell id="tc-input" title="Structures to compare" className="pc-t-input">
        <form onSubmit={handleCalculate} className="pc-tool-form">
          <SMILESInput
            value={smiles1}
            onChange={setSmiles1}
            id="smiles-input-1"
            label="Molecule 1 SMILES"
            required={true}
          />
          <SMILESInput
            value={smiles2}
            onChange={setSmiles2}
            id="smiles-input-2"
            label="Molecule 2 SMILES"
            required={true}
          />
          <div className="pc-opt-group">
            <span className="pc-lbl" id="tc-toolkit">
              Toolkit
            </span>
            <Seg
              labelId="tc-toolkit"
              value={toolkit}
              onChange={handleToolkitChange}
              options={[
                { value: "rdkit", label: "RDKit" },
                { value: "cdk", label: "CDK" },
              ]}
            />
          </div>
          <div className="pc-field-row">
            <Select
              label="Fingerprinter"
              value={fingerprinter}
              onChange={setFingerprinter}
              options={getFingerprinterOptions(toolkit)}
            />
            <Select
              label="Bits"
              value={String(nBits)}
              onChange={(v) => setNBits(Number(v))}
              options={[
                ["1024", "1024"],
                ["2048", "2048"],
                ["4096", "4096"],
              ]}
            />
            <Select
              label="Radius"
              value={String(radius)}
              onChange={(v) => setRadius(Number(v))}
              disabled={!usesRadius}
              options={[
                ["2", "2 (ECFP4)"],
                ["3", "3 (ECFP6)"],
                ["4", "4 (ECFP8)"],
              ]}
            />
          </div>
          {!usesRadius && <p className="pc-note">Radius only applies to ECFP and MAPC.</p>}
          <div className="pc-row pc-actions">
            <Button type="submit" disabled={loading || !smiles1.trim() || !smiles2.trim()}>
              {loading ? <Loader2 className="animate-spin" /> : <Calculator />}
              {loading ? "Calculating…" : "Calculate similarity"}
            </Button>
          </div>
        </form>
      </Cell>

      <Cell id="tc-about" title="About Tanimoto similarity" className="pc-t-side">
        <p className="pc-tool-note">
          The Tanimoto coefficient measures structural similarity between two molecules from their
          fingerprints. It runs from 0 (no common features) to 1 (identical fingerprints).
        </p>
        <dl className="pc-defs">
          <div>
            <dt>0.85 and above</dt>
            <dd>High similarity</dd>
          </div>
          <div>
            <dt>0.7 to 0.85</dt>
            <dd>Good similarity</dd>
          </div>
          <div>
            <dt>0.4 to 0.7</dt>
            <dd>Moderate similarity</dd>
          </div>
          <div>
            <dt>Below 0.4</dt>
            <dd>Low similarity</dd>
          </div>
        </dl>
      </Cell>

      <Cell
        id="tc-results"
        title={hasResult ? "Tanimoto similarity" : "Results"}
        className="pc-t-wide"
        meta={
          hasResult && (
            <span className="pc-tag">
              {fingerprinter} · {toolkit.toUpperCase()}
              {usesRadius ? ` · radius ${radius}` : ""} · {nBits} bits
            </span>
          )
        }
      >
        {loading && <ToolSkeleton variant="general" />}
        {error && !loading && (
          <GlassErrorCard
            message={error}
            onRetry={() => {
              setError(null);
              document.getElementById("smiles-input-1")?.focus();
            }}
          />
        )}
        {result === null && !loading && !error && (
          <p className="pc-empty">Enter two SMILES strings to compare them.</p>
        )}
        {hasResult && result !== null && (
          <div className="pc-tool-form">
            <div className="text-center">
              <div
                className={`mb-2 text-5xl font-bold whitespace-nowrap ${getSimilarityColor(result)}`}
              >
                {result.toFixed(4)}
              </div>
              <p className="text-sm text-muted-foreground">{getSimilarityDescription(result)}</p>
            </div>
            <div>
              <div className="h-3 w-full overflow-hidden rounded-full bg-accent">
                <div
                  className={`h-full rounded-full transition-all duration-500 ease-out ${getGaugeBgColor(result)}`}
                  style={{ width: `${Math.max(1, result * 100)}%` }}
                />
              </div>
              <div className="mt-1 flex justify-between px-1 text-xs text-muted-foreground">
                <span>0.0</span>
                <span>0.5</span>
                <span>1.0</span>
              </div>
            </div>
            <div className="pc-result-grid">
              <MoleculeCard smiles={smiles1} title="Molecule 1" size="lg" />
              <MoleculeCard smiles={smiles2} title="Molecule 2" size="lg" />
            </div>
          </div>
        )}
      </Cell>
    </div>
  );
};

export default TanimotoView;
