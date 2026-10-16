import { useState, useEffect, useCallback, type FormEvent } from "react";
import { Check, Loader2, RefreshCw, Search } from "lucide-react";
import SMILESInput from "../common/SMILESInput";
import MoleculeCard from "../common/MoleculeCard";
import { useAppContext } from "../../context/AppContext";
import { ToolSkeleton } from "@/components/feedback/ToolSkeleton";
import { GlassErrorCard } from "@/components/feedback/GlassErrorCard";
import { getErrorMessage } from "@/lib/error-messages";
import { Button } from "@/components/ui/button";
import { Cell, CopyButton, Seg } from "@/components/common/porcelain-ui";
import { useRevealOnResult } from "@/hooks/useRevealOnResult";

const POLLING_INTERVAL = 5000; // 5 seconds

const SERVICE_DOWN =
  "ClassyFire did not answer. The external service may be busy; try again in a few minutes.";

type Mode = "classify" | "jobid";

interface Level {
  name: string;
  description?: string;
}
interface Entity {
  smiles?: string;
  kingdom?: Level;
  superclass?: Level;
  class_?: Level;
  direct_parent?: Level;
  molecular_framework?: string;
  substituents?: string[];
  ancestors?: string[];
}
interface ClassyfireResult {
  classification_status?: string;
  entities?: Entity[];
}

/** An HTTP failure from the ClassyFire endpoints, with its status kept for the message. */
class ClassyfireHttpError extends Error {
  constructor(
    message: string,
    readonly status: number
  ) {
    super(message);
  }
}

/** A calm message: a missing job keeps its own text, bad input the parse text, the rest is "busy". */
function describeFailure(err: unknown): string {
  if (err instanceof ClassyfireHttpError) {
    if (err.status === 404) return err.message;
    if (err.status === 400 || err.status === 422) return getErrorMessage("chem", err);
    return SERVICE_DOWN;
  }
  // fetch rejects with a TypeError when the service cannot be reached at all
  if (err instanceof TypeError) return SERVICE_DOWN;
  return getErrorMessage("chem", err);
}

async function httpError(response: Response, fallback?: string) {
  let errorMsg = `Error ${response.status}: ${response.statusText}`;
  try {
    const errorData = await response.json();
    errorMsg = errorData.detail || errorMsg;
  } catch {
    // Ignore if response is not JSON
  }
  return new ClassyfireHttpError(fallback ?? errorMsg, response.status);
}

const LEVELS: [keyof Entity, string][] = [
  ["kingdom", "Kingdom"],
  ["superclass", "Superclass"],
  ["class_", "Class"],
  ["direct_parent", "Parent"],
];

/** ClassyFire: submit a structure to the external ClassyFire service and poll for its taxonomy. */
const ClassyfireView = () => {
  const [smiles, setSmiles] = useState("");
  const [jobId, setJobId] = useState<string | null>(null);
  const [jobStatus, setJobStatus] = useState<string | null>(null);
  const [classificationResults, setClassificationResults] = useState<ClassyfireResult | null>(null);
  const [loading, setLoading] = useState(false);
  const [polling, setPolling] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [activeTab, setActiveTab] = useState<Mode>("classify");
  const [manualJobId, setManualJobId] = useState("");
  const { addRecentMolecule, apiConfig } = useAppContext();

  // Submit classification request
  const handleSubmitClassification = async (e: FormEvent) => {
    e.preventDefault();
    const trimmedSmiles = smiles.trim();
    if (!trimmedSmiles) {
      setError("Please enter a SMILES string");
      return;
    }

    setLoading(true);
    setError(null);
    setJobId(null);
    setJobStatus(null);
    setClassificationResults(null);

    try {
      const response = await fetch(
        `${apiConfig.baseUrl}/chem/classyfire/classify?smiles=${encodeURIComponent(trimmedSmiles)}`
      );
      if (!response.ok) throw await httpError(response);

      const data = await response.json();

      if (data && data.id) {
        setJobId(String(data.id));
        setJobStatus("submitted");
        setPolling(true); // Start polling after successful submission

        addRecentMolecule({
          smiles: trimmedSmiles,
          name: "ClassyFire Molecule",
          timestamp: new Date().toISOString(),
        });
      } else {
        throw new Error("Invalid response from server: Missing job ID.");
      }
    } catch (err) {
      console.error("Classification submission error:", err);
      setError(describeFailure(err));
      setJobId(null);
      setJobStatus(null);
    } finally {
      setLoading(false);
    }
  };

  // Fetch job results
  const fetchClassificationResults = useCallback(
    async (id: string) => {
      try {
        const response = await fetch(`${apiConfig.baseUrl}/chem/classyfire/${id}/result`);

        if (!response.ok) {
          if (response.status === 404) {
            setJobId(null); // Clear invalid job ID
            setPolling(false); // Stop polling for invalid ID
            throw await httpError(
              response,
              `Job ID "${id}" not found. Please check the ID and try again.`
            );
          }
          throw await httpError(response);
        }

        const data = (await response.json()) as ClassyfireResult | null;

        if (data) {
          if (data.classification_status === "Done") {
            setClassificationResults(data);
            setJobStatus("completed");
            setPolling(false); // Stop polling only when successfully completed
          } else {
            setJobStatus(data.classification_status || "processing");
            setPolling(true);
          }
        } else {
          throw new Error("Received empty or invalid results data.");
        }
      } catch (err) {
        console.error(`Error fetching results for Job ID ${id}:`, err);
        setError(describeFailure(err));
        setPolling(false); // Stop polling on error
        setJobStatus("error");
      }
    },
    [apiConfig.baseUrl]
  );

  // Poll for results
  useEffect(() => {
    let intervalId: ReturnType<typeof setInterval> | null = null;

    if (polling && jobId) {
      const poll = () => {
        console.log(`Polling for Job ID: ${jobId}, Status: ${jobStatus}`);
        fetchClassificationResults(jobId);
      };

      // Fetch immediately when polling starts or jobId changes while polling is true
      poll();
      intervalId = setInterval(poll, POLLING_INTERVAL);
    }

    return () => {
      if (intervalId) {
        console.log(`Clearing interval for Job ID: ${jobId}`);
        clearInterval(intervalId);
      }
    };
  }, [polling, jobId, fetchClassificationResults, jobStatus]);

  // Handle manual job ID submission
  const handleSubmitJobId = (e: FormEvent) => {
    e.preventDefault();
    const trimmedJobId = manualJobId.trim();
    if (!trimmedJobId) {
      setError("Please enter a valid job ID");
      return;
    }
    setJobId(trimmedJobId);
    setPolling(true);
    setError(null);
    setClassificationResults(null);
    setJobStatus("fetching");
  };

  // Manually restart polling when it has stopped
  const handleManualRefresh = () => {
    if (jobId && !polling) {
      setError(null);
      setJobStatus("fetching");
      setPolling(true);
    }
  };

  const checkingThisJob = polling && jobId === manualJobId.trim();
  const entities = classificationResults?.entities ?? [];
  const busy = polling || ["fetching", "submitted", "processing"].includes(jobStatus ?? "");

  useRevealOnResult("cf-results", loading, !!jobId, !!classificationResults, !!error);

  return (
    <div className="pc-bento">
      <Cell id="cf-input" title="Structure input" className="pc-t-input">
        <div className="pc-tool-form">
          <div className="pc-opt-group">
            <span className="pc-lbl" id="cf-mode">
              Mode
            </span>
            <Seg
              labelId="cf-mode"
              full
              value={activeTab}
              onChange={setActiveTab}
              options={[
                { value: "classify", label: "Classify molecule" },
                { value: "jobid", label: "Check job status" },
              ]}
            />
          </div>

          {activeTab === "classify" ? (
            <form onSubmit={handleSubmitClassification} className="pc-tool-form">
              <SMILESInput value={smiles} onChange={setSmiles} required />
              <div className="pc-row pc-actions">
                <Button type="submit" disabled={!smiles.trim() || loading}>
                  {loading ? <Loader2 className="animate-spin" /> : <Search />}
                  {loading ? "Submitting…" : "Classify molecule"}
                </Button>
              </div>
            </form>
          ) : (
            <form onSubmit={handleSubmitJobId} className="pc-tool-form">
              <div className="grid gap-2">
                <label htmlFor="jobIdInput" className="pc-lbl">
                  ClassyFire job ID
                </label>
                <input
                  id="jobIdInput"
                  type="text"
                  className="pc-input"
                  value={manualJobId}
                  onChange={(e) => setManualJobId(e.target.value)}
                  placeholder="Enter job ID..."
                  spellCheck={false}
                  autoComplete="off"
                  required
                  aria-required="true"
                />
              </div>
              <div className="pc-row pc-actions">
                <Button type="submit" disabled={!manualJobId.trim() || loading || checkingThisJob}>
                  <Search />
                  {checkingThisJob ? "Checking…" : "Get classification results"}
                </Button>
              </div>
            </form>
          )}
        </div>
      </Cell>

      <Cell id="cf-about" title="About ClassyFire" className="pc-t-side">
        <p className="pc-tool-note">
          ClassyFire classifies chemical structures automatically with a computable chemical
          taxonomy from the Wishart Lab at the University of Alberta. Classification may take
          several minutes; keep the job ID to check the results later.
        </p>
        <dl className="pc-defs">
          <div>
            <dt>Levels</dt>
            <dd>
              Kingdom, superclass, class, subclass and further levels, down to the direct parent.
            </dd>
          </div>
        </dl>
      </Cell>

      <Cell
        id="cf-results"
        title={classificationResults ? "Classification" : "Results"}
        className="pc-t-wide"
        meta={
          classificationResults && (
            <span className="pc-tag pc-row text-primary">
              <Check className="size-3.5" aria-hidden="true" />
              {classificationResults.classification_status || "Completed"}
            </span>
          )
        }
      >
        {loading && !classificationResults && <ToolSkeleton variant="molecule" />}
        {error && !loading && (
          <GlassErrorCard
            message={error}
            onRetry={() => {
              setError(null);
              document
                .getElementById(activeTab === "classify" ? "smiles-input" : "jobIdInput")
                ?.focus();
            }}
          />
        )}

        {jobId && !loading && !classificationResults && (
          <div className="pc-tool-form" role="status">
            <p className="pc-row text-sm text-foreground">
              {busy ? (
                <Loader2 className="size-4 animate-spin text-primary" aria-hidden="true" />
              ) : (
                <RefreshCw className="size-4 text-primary" aria-hidden="true" />
              )}
              <span>
                Job <span className="pc-mono">{jobId}</span> · status{" "}
                <span className="font-medium capitalize">{jobStatus || "waiting"}</span>
              </span>
            </p>
            <p className="pc-note">
              {jobStatus === "error"
                ? "The job could not be checked; see the message above."
                : polling
                  ? "ClassyFire is working on it. The results appear here on their own."
                  : "Checking has stopped. Refresh to check the status again."}
            </p>
            <div className="pc-row">
              <CopyButton text={jobId} label="Copy job ID" />
              {!polling && jobStatus !== "completed" && jobStatus !== "error" && (
                <Button variant="outline" onClick={handleManualRefresh}>
                  <RefreshCw />
                  Refresh status
                </Button>
              )}
            </div>
          </div>
        )}

        {!classificationResults && !loading && !error && !jobId && (
          <p className="pc-empty">
            Enter a SMILES string to classify it, or check an earlier job by its ID.
          </p>
        )}

        {classificationResults &&
          (entities.length > 0 ? (
            <div className="pc-tool-form">
              {entities.map((entity, index) => {
                const levels = LEVELS.filter(([key]) => entity[key]);
                return (
                  <div
                    key={entity.smiles || index}
                    className="grid gap-4 border-b border-border pb-6 last:border-b-0 last:pb-0 lg:grid-cols-3"
                  >
                    <MoleculeCard
                      smiles={entity.smiles ?? ""}
                      title={`Input Structure ${index + 1}`}
                      size="sm"
                    />
                    <div className="grid content-start gap-4 lg:col-span-2">
                      <h3 className="text-sm font-medium text-muted-foreground">
                        Chemical ontology
                      </h3>
                      {levels.length > 0 ? (
                        <dl className="pc-defs">
                          {levels.map(([key, label]) => {
                            const level = entity[key] as Level;
                            return (
                              <div key={key}>
                                <dt>{label}</dt>
                                <dd>
                                  <span className="font-medium text-foreground">{level.name}</span>
                                  {level.description && (
                                    <span className="mt-1 block text-sm text-muted-foreground">
                                      {level.description}
                                    </span>
                                  )}
                                </dd>
                              </div>
                            );
                          })}
                        </dl>
                      ) : (
                        <p className="pc-note">No classification levels found.</p>
                      )}

                      {(entity.molecular_framework ||
                        (entity.substituents && entity.substituents.length > 0) ||
                        (entity.ancestors && entity.ancestors.length > 0)) && (
                        <>
                          <h3 className="text-sm font-medium text-muted-foreground">
                            Additional details
                          </h3>
                          <dl className="pc-defs">
                            {entity.molecular_framework && (
                              <div>
                                <dt>Molecular framework</dt>
                                <dd>{entity.molecular_framework}</dd>
                              </div>
                            )}
                            {entity.substituents && entity.substituents.length > 0 && (
                              <div>
                                <dt>Substituents</dt>
                                <dd className="pc-row mt-1 gap-2">
                                  {entity.substituents.map((substituent, idx) => (
                                    <span key={idx} className="pc-tag">
                                      {substituent}
                                    </span>
                                  ))}
                                </dd>
                              </div>
                            )}
                            {entity.ancestors && entity.ancestors.length > 0 && (
                              <div>
                                <dt>Classification lineage</dt>
                                <dd>
                                  <ol className="m-0 flex list-none flex-wrap items-center gap-x-1 p-0 text-sm">
                                    {entity.ancestors.map((ancestor, idx, all) => (
                                      <li key={idx} className="flex items-center">
                                        <span className="text-foreground">{ancestor}</span>
                                        {idx < all.length - 1 && (
                                          <span
                                            className="mx-1 text-muted-foreground"
                                            aria-hidden="true"
                                          >
                                            &gt;
                                          </span>
                                        )}
                                      </li>
                                    ))}
                                  </ol>
                                </dd>
                              </div>
                            )}
                          </dl>
                        </>
                      )}
                    </div>
                  </div>
                );
              })}
            </div>
          ) : (
            <p className="pc-empty">No classification entities found in the results.</p>
          ))}
      </Cell>
    </div>
  );
};

export default ClassyfireView;
