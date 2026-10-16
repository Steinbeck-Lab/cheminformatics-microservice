import { useState, type FormEvent } from "react";
import { Check, ExternalLink, Info, Loader2, Search } from "lucide-react";
import MoleculeCard from "../common/MoleculeCard";
import SMILESDisplay from "../common/SMILESDisplay";
import { useAppContext } from "../../context/AppContext";
import { lookupPubChem } from "../../services/chemService";
import type { PubChemLookupResult } from "@/types/api";
import { ToolSkeleton } from "@/components/feedback/ToolSkeleton";
import { GlassErrorCard } from "@/components/feedback/GlassErrorCard";
import { getErrorMessage } from "@/lib/error-messages";
import { Button } from "@/components/ui/button";
import { Cell } from "@/components/common/porcelain-ui";
import { useRevealOnResult } from "@/hooks/useRevealOnResult";

/** The endpoint also sends the matching CIDs and their PubChem pages. */
type LookupResult = PubChemLookupResult & { cids?: number[]; pubchem_links?: string[] };

const EXAMPLES = [
  { name: "Chemical Name", value: "aspirin", description: "Common or IUPAC name" },
  { name: "CID", value: "2244", description: "PubChem Compound ID" },
  {
    name: "InChI",
    value: "InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)",
    description: "International Chemical Identifier",
  },
  { name: "InChIKey", value: "BSYNRYMUTXBXSQ-UHFFFAOYSA-N", description: "Hashed InChI" },
  { name: "CAS", value: "50-78-2", description: "Chemical Abstracts Service registry number" },
  { name: "Formula", value: "C9H8O4", description: "Molecular formula" },
  {
    name: "SMILES",
    value: "CC(=O)OC1=CC=CC=C1C(=O)O",
    description: "Simplified molecular-input line-entry system",
  },
];

/** PubChem lookup: a structure from a name, CID, InChI, InChIKey, CAS number, formula or SMILES. */
const PubChemLookupView = () => {
  const [identifier, setIdentifier] = useState("");
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [result, setResult] = useState<LookupResult | null>(null);
  const { addRecentMolecule } = useAppContext();

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
      if (data.success && data.canonical_smiles) {
        addRecentMolecule({
          smiles: data.canonical_smiles,
          name: trimmedIdentifier,
          timestamp: new Date().toISOString(),
        });
      }
    } catch (err) {
      console.error("PubChem lookup error:", err);
      setError(getErrorMessage("chem", err));
    } finally {
      setLoading(false);
    }
  };

  const done = !!result && !loading;
  const title = !done ? "Results" : result.success ? "PubChem match" : "No match";
  const links = result?.pubchem_links ?? [];

  useRevealOnResult("pc-lookup-results", loading, !!result, !!error);

  return (
    <div className="pc-bento">
      <Cell id="pc-lookup-input" title="Structure input" className="pc-t-input">
        <form onSubmit={handleSubmit} className="pc-tool-form">
          <div className="pc-smiles">
            <label htmlFor="identifier-input" className="pc-lbl">
              Chemical identifier <span className="text-madder-ink">*</span>
            </label>
            <div className="pc-smiles-field">
              <input
                id="identifier-input"
                type="text"
                className="pc-well"
                value={identifier}
                onChange={(e) => setIdentifier(e.target.value)}
                placeholder="Name, CID, InChI, InChIKey, CAS, formula or SMILES"
                spellCheck={false}
                autoComplete="off"
                autoCapitalize="off"
                required
                aria-required
              />
            </div>
            <div className="pc-smiles-chips">
              <span className="pc-lbl2" id="identifier-ex">
                Examples
              </span>
              <div className="pc-row" role="group" aria-labelledby="identifier-ex">
                {EXAMPLES.map((example) => (
                  <button
                    key={example.name}
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
          </div>
          <div className="pc-row pc-actions">
            <Button type="submit" disabled={!identifier.trim() || loading}>
              {loading ? <Loader2 className="animate-spin" /> : <Search />}
              {loading ? "Searching…" : "Find structure"}
            </Button>
          </div>
        </form>
      </Cell>

      <Cell id="pc-lookup-about" title="About PubChem lookup" className="pc-t-side">
        <p className="pc-tool-note">
          Finds a structure in PubChem, which holds over 100 million compounds. The input type is
          detected for you.
        </p>
        <dl className="pc-defs">
          <div>
            <dt>Accepted input</dt>
            <dd>
              Names (common, trade or IUPAC), PubChem CIDs, InChI, InChIKey, CAS numbers, molecular
              formulas and SMILES.
            </dd>
          </div>
        </dl>
      </Cell>

      <Cell id="pc-lookup-results" title={title} className="pc-t-wide">
        {loading && <ToolSkeleton variant="descriptors" />}
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
          <p className="pc-empty">Enter a name or identifier to look up its structure.</p>
        )}
        {done && (
          <div className="pc-tool-form">
            <div className={result.success ? "text-primary" : "text-saffron-ink"}>
              <p className="pc-row text-sm font-medium">
                {result.success ? (
                  <Check className="size-4 shrink-0" aria-hidden="true" />
                ) : (
                  <Info className="size-4 shrink-0" aria-hidden="true" />
                )}
                {result.success
                  ? "Successfully retrieved structure from PubChem"
                  : "No structure found in PubChem for this identifier"}
              </p>
            </div>
            <p className="text-sm text-muted-foreground">
              Input identifier: <span className="font-mono text-foreground">{result.input}</span>{" "}
              (detected as {result.input_type})
            </p>
            {result.success && result.cids && (
              <div className="pc-tool-form gap-2">
                <p className="text-sm font-medium text-foreground">
                  PubChem CID{result.cids.length > 1 ? "s" : ""}:{" "}
                  {result.cids.slice(0, 5).join(", ")}
                  {result.cids.length > 5 && ` and ${result.cids.length - 5} more`}
                </p>
                <div className="pc-row">
                  {links.slice(0, 3).map((link, index) => (
                    <a
                      key={link}
                      href={link}
                      target="_blank"
                      rel="noopener noreferrer"
                      className="pc-row gap-1 text-sm font-medium text-primary underline-offset-4 hover:underline"
                    >
                      View CID {result.cids?.[index]} in PubChem
                      <ExternalLink className="size-3.5" aria-hidden="true" />
                    </a>
                  ))}
                  {links.length > 3 && (
                    <span className="pc-tag">+{links.length - 3} more compounds</span>
                  )}
                </div>
              </div>
            )}
            {result.success && result.canonical_smiles && (
              <>
                <SMILESDisplay
                  smiles={result.canonical_smiles}
                  label="Canonical SMILES"
                  showDownload={true}
                />
                <div className="flex justify-center">
                  <MoleculeCard smiles={result.canonical_smiles} title={result.input} size="lg" />
                </div>
              </>
            )}
          </div>
        )}
      </Cell>
    </div>
  );
};

export default PubChemLookupView;
