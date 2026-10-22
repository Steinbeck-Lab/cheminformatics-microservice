import { useState, type FormEvent } from "react";
import { FlaskConical, Loader2 } from "lucide-react";
import SMILESInput from "../common/SMILESInput";
import { coconutPreprocessing } from "../../services/chemService";
import { ToolSkeleton } from "@/components/feedback/ToolSkeleton";
import { GlassErrorCard } from "@/components/feedback/GlassErrorCard";
import { getErrorMessage } from "@/lib/error-messages";
import { Button } from "@/components/ui/button";
import { Cell, Seg, Switch } from "@/components/common/porcelain-ui";
import { useRevealOnResult } from "@/hooks/useRevealOnResult";

/** One structure version (original, standardized, parent) from /chem/coconut/pre-processing. */
interface StructureVersion {
  representations?: Record<string, string>;
  has_stereo?: boolean;
  has_stereo_defined?: boolean;
  has_stereogenic_elements?: boolean;
  descriptors?: Record<string, unknown>;
  errors?: unknown;
}
type PreprocessData = Record<string, StructureVersion | undefined>;

/** 'some_key' -> 'Some Key' */
const formatKey = (key: string) => key.replace(/_/g, " ").replace(/\b\w/g, (l) => l.toUpperCase());

const formatValue = (value: unknown) => {
  if (typeof value === "boolean") return value ? "Yes" : "No";
  if (typeof value === "number") return parseFloat(value.toFixed(4)).toString();
  if (typeof value === "object" && value !== null) return JSON.stringify(value);
  return String(value);
};

const subhead = "text-sm font-medium text-muted-foreground";

/** Yes/No property row. */
function Flag({ label, value }: { label: string; value?: boolean }) {
  return (
    <div>
      <dt>{label}</dt>
      <dd className={value ? "text-primary" : "text-madder-ink"}>{value ? "Yes" : "No"}</dd>
    </div>
  );
}

/** COCONUT pre-processing: original, standardized and parent structure for COCONUT submission. */
const CoconutPreProcessingView = () => {
  const [smiles, setSmiles] = useState("");
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [preprocessData, setPreprocessData] = useState<PreprocessData | null>(null);
  const [generate3D, setGenerate3D] = useState(false);
  const [generateDescriptors, setGenerateDescriptors] = useState(false);
  const [activeTab, setActiveTab] = useState("original");

  const handleSubmit = async (e: FormEvent) => {
    e.preventDefault();
    const trimmedSmiles = smiles.trim();
    if (!trimmedSmiles) {
      setError("Please enter a SMILES string");
      return;
    }

    setLoading(true);
    setError(null);
    setPreprocessData(null);
    setActiveTab("original");

    try {
      const data = await coconutPreprocessing(trimmedSmiles, generate3D, generateDescriptors);
      setPreprocessData(data as PreprocessData);
    } catch (err) {
      console.error("Preprocessing error:", err);
      setError(getErrorMessage("chem", err));
    } finally {
      setLoading(false);
    }
  };

  const renderTabContent = (tabData?: StructureVersion) => {
    if (!tabData) return <p className="pc-empty">No data available for this structure.</p>;

    const reps = Object.entries(tabData.representations ?? {});
    const molblocks = reps.filter(([key]) => key.includes("MOL"));
    const otherReps = reps.filter(([key]) => !key.includes("MOL"));
    const descriptors =
      typeof tabData.descriptors === "object" && tabData.descriptors !== null
        ? Object.entries(tabData.descriptors)
        : [];

    return (
      <div className="pc-tool-form">
        <div className="grid gap-4 lg:grid-cols-2">
          <div className="grid content-start gap-2">
            <h3 className={subhead}>Representations</h3>
            {otherReps.length > 0 ? (
              <dl className="pc-ids">
                {otherReps.map(([key, value]) => (
                  <div className="pc-id" key={key}>
                    <dt>{formatKey(key)}</dt>
                    <dd>{value}</dd>
                  </div>
                ))}
              </dl>
            ) : (
              <p className="pc-note">No other representations available.</p>
            )}
          </div>

          <div className="grid content-start gap-2">
            <h3 className={subhead}>Properties</h3>
            <dl className="pc-defs">
              <Flag label="Has stereo" value={tabData.has_stereo} />
              {tabData.has_stereo_defined !== undefined && (
                <Flag label="Has stereo defined" value={tabData.has_stereo_defined} />
              )}
              {tabData.has_stereogenic_elements !== undefined && (
                <Flag label="Has stereogenic elements" value={tabData.has_stereogenic_elements} />
              )}
            </dl>
            {tabData.errors != null && (
              <>
                <h3 className={subhead}>Validation errors</h3>
                {Array.isArray(tabData.errors) && tabData.errors.length > 0 ? (
                  <ul className="list-disc space-y-1 pl-5 text-sm text-madder-ink">
                    {tabData.errors.map((err, idx) => (
                      <li key={idx}>{formatValue(err)}</li>
                    ))}
                  </ul>
                ) : (
                  <p className="text-sm text-primary">No errors found for this structure.</p>
                )}
              </>
            )}
          </div>
        </div>

        {molblocks.length > 0 && (
          <div className="grid gap-3">
            <h3 className={subhead}>Molecular structures</h3>
            {molblocks.map(([key, value]) => (
              <div key={key} className="grid gap-1">
                <span className="pc-lbl">{formatKey(key)}</span>
                <pre className="max-h-80 overflow-auto rounded-xl border border-border bg-muted p-3 font-mono text-xs whitespace-pre text-foreground">
                  {value}
                </pre>
              </div>
            ))}
          </div>
        )}

        {generateDescriptors && tabData.descriptors && (
          <div className="grid gap-2">
            <h3 className={subhead}>Descriptors</h3>
            {descriptors.length > 0 ? (
              <div className="overflow-x-auto rounded-xl border border-border">
                <table className="w-full min-w-[400px] border-collapse">
                  <thead className="bg-muted">
                    <tr>
                      <th
                        scope="col"
                        className="px-4 py-2 text-left text-xs font-medium tracking-wider text-muted-foreground uppercase"
                      >
                        Descriptor
                      </th>
                      <th
                        scope="col"
                        className="px-4 py-2 text-right text-xs font-medium tracking-wider text-muted-foreground uppercase"
                      >
                        Value
                      </th>
                    </tr>
                  </thead>
                  <tbody className="divide-y divide-border">
                    {descriptors.map(([key, value]) => (
                      <tr key={key}>
                        <td className="px-4 py-2 text-sm text-foreground">{formatKey(key)}</td>
                        <td className="px-4 py-2 text-right font-mono text-sm text-foreground tabular-nums">
                          {formatValue(value)}
                        </td>
                      </tr>
                    ))}
                  </tbody>
                </table>
              </div>
            ) : (
              <p className="pc-note">No descriptors were calculated or available.</p>
            )}
          </div>
        )}
      </div>
    );
  };

  const done = !!preprocessData && !loading;
  const tabs = preprocessData ? Object.keys(preprocessData) : [];

  useRevealOnResult("cp-results", loading, !!preprocessData, !!error);

  return (
    <div className="pc-bento">
      <Cell id="cp-input" title="Structure input" className="pc-t-input">
        <form onSubmit={handleSubmit} className="pc-tool-form">
          <SMILESInput value={smiles} onChange={setSmiles} required />
          <div className="pc-opt-list">
            <Switch checked={generate3D} onChange={setGenerate3D} label="Generate 3D coordinates" />
            <Switch
              checked={generateDescriptors}
              onChange={setGenerateDescriptors}
              label="Calculate descriptors"
            />
          </div>
          <div className="pc-row pc-actions">
            <Button type="submit" disabled={!smiles.trim() || loading}>
              {loading ? <Loader2 className="animate-spin" /> : <FlaskConical />}
              {loading ? "Processing…" : "Process for COCONUT"}
            </Button>
          </div>
        </form>
      </Cell>

      <Cell id="cp-about" title="About COCONUT pre-processing" className="pc-t-side">
        <p className="pc-tool-note">
          Prepares a molecule for submission to COCONUT, the COlleCtion of Open Natural prodUcTs.
          Returns the original, the ChEMBL-standardized and the parent (salt and solvent stripped)
          structure.
        </p>
        <dl className="pc-defs">
          <div>
            <dt>Representations</dt>
            <dd>SMILES, InChI and others for each structure.</dd>
          </div>
          <div>
            <dt>Optional</dt>
            <dd>3D coordinates and molecular descriptors.</dd>
          </div>
        </dl>
      </Cell>

      <Cell
        id="cp-results"
        title={done ? `${formatKey(activeTab)} structure` : "Results"}
        className="pc-t-wide"
        meta={
          done &&
          tabs.length > 0 && (
            <Seg
              label="Structure version"
              value={activeTab}
              onChange={setActiveTab}
              options={tabs.map((key) => ({ value: key, label: formatKey(key) }))}
            />
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
        {!preprocessData && !loading && !error && (
          <p className="pc-empty">Enter a SMILES string to pre-process it for COCONUT.</p>
        )}
        {done && preprocessData && renderTabContent(preprocessData[activeTab])}
      </Cell>
    </div>
  );
};

export default CoconutPreProcessingView;
