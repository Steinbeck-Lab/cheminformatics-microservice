// Structure Explorer: PubChem lookup by name, CAS, formula, InChI or SMILES, then 2D and 3D views.
import { useState, type FormEvent } from "react";
import { Loader2, Search } from "lucide-react";
import SMILESDisplay from "../common/SMILESDisplay";
import { ToolSkeleton } from "@/components/feedback/ToolSkeleton";
import { GlassErrorCard } from "@/components/feedback/GlassErrorCard";
import { getErrorMessage } from "@/lib/error-messages";
import { useAppContext } from "../../context/AppContext";
import { lookupPubChem } from "../../services/chemService";
import type { PubChemLookupResult } from "@/types/api";
import MoleculeDepiction2D from "./MoleculeDepiction2D";
import MoleculeDepiction3D from "./MoleculeDepiction3D";
import { Button } from "@/components/ui/button";
import { Cell, Seg } from "@/components/common/porcelain-ui";
import { useRevealOnResult } from "@/hooks/useRevealOnResult";

type Toolkit2D = "rdkit" | "cdk";
type Toolkit3D = "openbabel" | "rdkit";
/** The lookup may also return the compound name. */
type LookupResult = PubChemLookupResult & { name?: string };

const TOOLKITS_2D: { value: Toolkit2D; label: string }[] = [
  { value: "rdkit", label: "RDKit" },
  { value: "cdk", label: "CDK" },
];
const TOOLKITS_3D: { value: Toolkit3D; label: string }[] = [
  { value: "openbabel", label: "Open Babel" },
  { value: "rdkit", label: "RDKit" },
];

const EXAMPLES = [
  { name: "Caffeine", value: "caffeine", description: "Common stimulant" },
  { name: "Aspirin", value: "aspirin", description: "Common pain reliever" },
  { name: "Ibuprofen", value: "ibuprofen", description: "NSAID pain reliever" },
  { name: "Paracetamol", value: "paracetamol", description: "Also known as acetaminophen" },
  { name: "CAS: 50-78-2", value: "50-78-2", description: "CAS for aspirin" },
  { name: "Formula: C8H10N4O2", value: "C8H10N4O2", description: "Caffeine formula" },
];

const StructureVisualizerView = () => {
  const [identifier, setIdentifier] = useState("");
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [result, setResult] = useState<LookupResult | null>(null);
  const { addRecentMolecule } = useAppContext();

  const [depictToolkit, setDepictToolkit] = useState<Toolkit2D>("rdkit");
  const [vis3DToolkit, setVis3DToolkit] = useState<Toolkit3D>("openbabel");

  const handleSubmit = async (e: FormEvent) => {
    e.preventDefault();
    const trimmedIdentifier = identifier.trim();
    if (!trimmedIdentifier) {
      setError("Please enter a chemical identifier");
      return;
    }

    setLoading(true);
    setError(null);
    setResult(null);

    try {
      const data: LookupResult = await lookupPubChem(trimmedIdentifier);
      setResult(data);

      // Add to recent molecules if the lookup was successful
      if (data.success && data.canonical_smiles) {
        addRecentMolecule({
          smiles: data.canonical_smiles,
          name: data.name || trimmedIdentifier,
          timestamp: new Date().toISOString(),
        });
      }
    } catch (err) {
      console.error("PubChem lookup error:", err);
      setError(getErrorMessage("depict", err));
    } finally {
      setLoading(false);
    }
  };

  const found = !!result?.success && !!result.canonical_smiles && !loading;
  const label = result?.name || result?.input;

  useRevealOnResult("sv-results", loading, !!result, !!error);

  return (
    <div className="pc-bento">
      <Cell id="sv-input" title="Structure input" className="pc-t-input">
        <form onSubmit={handleSubmit} className="pc-tool-form">
          <div className="grid gap-2">
            <label htmlFor="identifier-input" className="pc-lbl">
              Chemical identifier
            </label>
            <input
              id="identifier-input"
              type="text"
              className="pc-well"
              value={identifier}
              onChange={(e) => setIdentifier(e.target.value)}
              placeholder="Name, CAS number, formula, InChI or SMILES"
              spellCheck={false}
              autoComplete="off"
              required
            />
          </div>
          <div className="pc-smiles-chips">
            <span className="pc-lbl2" id="sv-examples">
              Examples
            </span>
            <div className="pc-row" role="group" aria-labelledby="sv-examples">
              {EXAMPLES.map((example) => (
                <button
                  key={example.value}
                  type="button"
                  className="clay clay-outline pc-chip"
                  aria-pressed={identifier === example.value}
                  title={`${example.name}: ${example.description}`}
                  onClick={() => setIdentifier(example.value)}
                >
                  {example.name}
                </button>
              ))}
            </div>
          </div>
          <div className="pc-row pc-actions">
            <Button type="submit" disabled={!identifier.trim() || loading}>
              {loading ? <Loader2 className="animate-spin" /> : <Search />}
              {loading ? "Searching…" : "Search and visualize"}
            </Button>
          </div>
        </form>
      </Cell>

      <Cell id="sv-opts" title="Options" className="pc-t-side">
        <div className="pc-tool-form">
          <div className="pc-opt-group">
            <span className="pc-lbl" id="sv-tk2d">
              2D toolkit
            </span>
            <Seg
              labelId="sv-tk2d"
              full
              value={depictToolkit}
              onChange={setDepictToolkit}
              options={TOOLKITS_2D}
            />
          </div>
          <div className="pc-opt-group">
            <span className="pc-lbl" id="sv-tk3d">
              3D toolkit
            </span>
            <Seg
              labelId="sv-tk3d"
              full
              value={vis3DToolkit}
              onChange={setVis3DToolkit}
              options={TOOLKITS_3D}
            />
          </div>
          <p className="pc-tool-note">
            Looks the structure up in PubChem by name, CAS number, formula, InChI or SMILES. Shows
            it as an exact 2D diagram and an interactive 3D model.
          </p>
        </div>
      </Cell>

      <Cell id="sv-results" title={found ? "Structure" : "Results"} className="pc-t-wide">
        {loading && <ToolSkeleton variant="molecule" />}
        {error && !loading && (
          <GlassErrorCard
            message={error}
            onRetry={() => {
              setError(null);
              document.getElementById("identifier-input")?.focus();
            }}
          />
        )}
        {!result && !loading && !error && (
          <p className="pc-empty">
            Search PubChem by name, formula, CAS number or SMILES to see the structure in 2D and 3D.
          </p>
        )}
        {result && !loading && (
          <div className="pc-tool-form">
            <div className="grid gap-1">
              <p
                className={`m-0 font-semibold ${result.success ? "text-primary" : "text-saffron-ink"}`}
              >
                {result.success ? `Found: ${label}` : "No structure found in PubChem"}
              </p>
              <p className="pc-note m-0">
                Input identifier: <span className="pc-mono">{result.input}</span>{" "}
                <span className="pc-tag">{result.input_type}</span>
              </p>
            </div>
            {found && (
              <>
                <SMILESDisplay
                  smiles={result.canonical_smiles}
                  label="Canonical SMILES"
                  showDownload={true}
                />
                <div className="grid gap-4 lg:grid-cols-2">
                  <MoleculeDepiction2D
                    smiles={result.canonical_smiles}
                    title={`2D Structure: ${label}`}
                    toolkit={depictToolkit}
                  />
                  <MoleculeDepiction3D
                    smiles={result.canonical_smiles}
                    title={`3D Structure: ${label}`}
                    toolkit={vis3DToolkit}
                  />
                </div>
              </>
            )}
          </div>
        )}
      </Cell>
    </div>
  );
};

export default StructureVisualizerView;
