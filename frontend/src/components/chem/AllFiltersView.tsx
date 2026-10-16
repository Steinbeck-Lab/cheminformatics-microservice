import { useState, type FormEvent } from "react";
import { Check, FileDown, Filter, FlaskConical, Loader2, X } from "lucide-react";
import { useAppContext } from "../../context/AppContext";
import api from "../../services/api";
import { applyChemicalFiltersDetailed } from "../../services/toolsService";
import { ToolSkeleton } from "@/components/feedback/ToolSkeleton";
import { GlassErrorCard } from "@/components/feedback/GlassErrorCard";
import { getErrorMessage } from "@/lib/error-messages";
import { Button } from "@/components/ui/button";
import { Cell, CopyButton, Select, Seg, Switch } from "@/components/common/porcelain-ui";
import { useRevealOnResult } from "@/hooks/useRevealOnResult";

type BoolFilter = "pains" | "lipinski" | "veber" | "reos" | "ghose" | "ruleofthree";
type RangeFilter = "qedscore" | "sascore" | "nplikeness";
type FilterOptions = Record<BoolFilter, boolean> & Record<RangeFilter, string>;
type Operator = "AND" | "OR";

interface FilterDetail {
  passes: boolean;
  details?: string;
  violations?: string[];
  properties?: Record<string, unknown>;
}
interface DetailedResult {
  smiles: string;
  valid: boolean;
  error?: string;
  overall_pass?: boolean;
  filters: Record<string, FilterDetail>;
}
interface DetailedResults {
  total_molecules?: number;
  passing_molecules?: number;
  filter_operator?: string;
  results: DetailedResult[];
}

/** Table columns, in the order the backend returns them. */
const filterColumns: { key: BoolFilter | RangeFilter; label: string }[] = [
  { key: "pains", label: "PAINS" },
  { key: "lipinski", label: "Lipinski" },
  { key: "veber", label: "Veber" },
  { key: "reos", label: "REOS" },
  { key: "ghose", label: "Ghose" },
  { key: "ruleofthree", label: "Rule of 3" },
  { key: "qedscore", label: "QED" },
  { key: "sascore", label: "SA Score" },
  { key: "nplikeness", label: "NP-like" },
];

const RULES: { key: BoolFilter; label: string; hint: string }[] = [
  {
    key: "pains",
    label: "PAINS",
    hint: "Pan-assay interference compounds; a match is undesirable",
  },
  {
    key: "lipinski",
    label: "Lipinski",
    hint: "Rule of 5: weight, logP, H-bond donors and acceptors",
  },
  {
    key: "veber",
    label: "Veber",
    hint: "Oral bioavailability: rotatable bonds, polar surface area",
  },
  { key: "reos", label: "REOS", hint: "Rapid elimination of swill: lead-like properties" },
  { key: "ghose", label: "Ghose", hint: "logP, molecular weight and number of atoms" },
  { key: "ruleofthree", label: "Rule of 3", hint: "Criteria for fragment-based drug discovery" },
];

const RANGES: { key: RangeFilter; label: string; options: [string, string][] }[] = [
  {
    key: "qedscore",
    label: "QED drug-likeness (0-1)",
    options: [
      ["0-1", "0-1 (All)"],
      ["0-0.3", "0-0.3 (Low)"],
      ["0.3-0.6", "0.3-0.6 (Medium)"],
      ["0.6-1", "0.6-1 (High)"],
    ],
  },
  {
    key: "sascore",
    label: "Synthetic accessibility (1-10)",
    options: [
      ["1-10", "1-10 (All)"],
      ["1-4", "1-4 (Easy synthesis)"],
      ["4-7", "4-7 (Moderate synthesis)"],
      ["7-10", "7-10 (Difficult synthesis)"],
    ],
  },
  {
    key: "nplikeness",
    label: "NP-likeness score (-5 to 5)",
    options: [
      ["-5-5", "-5 to 5 (All)"],
      ["-5-0", "-5 to 0 (Synthetic-like)"],
      ["0-2", "0 to 2 (Moderate NP-like)"],
      ["2-5", "2 to 5 (High NP-like)"],
    ],
  },
];

const EXAMPLES = [
  "CC1=CC=C(C=C1)C(=O)NC=CC(=O)O", // Acetaminophen
  "CN1C=NC2=C1C(=O)N(C(=O)N2C)C", // Caffeine
  "CC(=O)OC1=CC=CC=C1C(=O)O", // Aspirin
  "COC1=CC=C(CCN(C)C)C=C1", // Metoclopramide
  "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O", // Ibuprofen
  "CCN(CC)CCOC(=O)C1=CC=CC=C1N", // Procaine
].join("\n");

const th =
  "px-3 py-2 text-xs font-medium tracking-wider whitespace-nowrap text-muted-foreground uppercase";

/** Pass/fail tag in the detailed view. */
const PassTag = ({ pass, children }: { pass?: boolean; children: string }) => (
  <span className={`pc-tag ${pass ? "text-primary" : "text-madder-ink"}`}>{children}</span>
);

/** Chemical filters: screen a list of molecules against medicinal-chemistry rules and scores. */
const AllFiltersView = () => {
  const [smilesInput, setSmilesInput] = useState("");
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [results, setResults] = useState<string[]>([]);
  const [detailedResults, setDetailedResults] = useState<DetailedResults | null>(null);
  const [showDetailedView, setShowDetailedView] = useState(false);
  const [filterOperator, setFilterOperator] = useState<Operator>("OR");
  const [ran, setRan] = useState(false);
  const { addRecentMolecule } = useAppContext();

  const [filterOptions, setFilterOptions] = useState<FilterOptions>({
    pains: true,
    lipinski: true,
    veber: true,
    reos: true,
    ghose: true,
    ruleofthree: true,
    qedscore: "0-1",
    sascore: "1-10",
    nplikeness: "-5-5",
  });

  const applyFilters = async (
    e?: FormEvent | null,
    overrideOperator: Operator | null = null,
    explicitOptions: FilterOptions | null = null
  ) => {
    e?.preventDefault();

    const trimmedInput = smilesInput.trim();
    if (!trimmedInput) {
      setError("Please enter at least one SMILES string.");
      return;
    }

    setLoading(true);
    setError(null);
    setResults([]);
    setDetailedResults(null);
    setRan(true);

    const optionsForApi = explicitOptions || filterOptions;

    try {
      const apiParams = {
        ...optionsForApi,
        filterOperator: overrideOperator || filterOperator,
      };

      console.log("Applying filters with parameters:", apiParams);
      console.log(`Filter operator: ${apiParams.filterOperator}`);
      console.log(`Detailed view: ${showDetailedView}`);

      let responseData: string[] = [];

      if (showDetailedView) {
        const detailedData = (await applyChemicalFiltersDetailed(
          trimmedInput,
          apiParams
        )) as DetailedResults;
        console.log("Received detailed filter results:", detailedData);
        setDetailedResults(detailedData);

        // Also set simple results for the table
        const simpleResults = detailedData.results.map((result) => {
          if (!result.valid) return `${result.smiles} : ERROR`;
          const filterResults = Object.values(result.filters).map((f) => (f.passes ? "T" : "F"));
          return `${result.smiles} : ${filterResults.join(", ")}`;
        });
        setResults(simpleResults);
        responseData = simpleResults;
      } else {
        const response = await api.post("/chem/all_filters", trimmedInput, {
          params: apiParams,
          headers: {
            "Content-Type": "text/plain",
            Accept: "application/json",
          },
        });

        if (Array.isArray(response.data)) {
          console.log("Received filter results:", response.data);
          setResults(response.data);
          responseData = response.data;
        } else {
          console.warn("Received non-array response:", response.data);
          setResults([]);
          setError("Received unexpected data format from the server.");
        }
      }

      // Add molecules to the recent list
      if (responseData.length > 0) {
        trimmedInput
          .split(/[\n\s,;]+/)
          .filter((s) => s.trim())
          .forEach((smiles) => {
            const smilesOnly = smiles.trim().split(/[\s:]+/)[0];
            if (smilesOnly) {
              addRecentMolecule({
                smiles: smilesOnly,
                name: `Filtered Molecule (${new Date().toLocaleTimeString()})`,
                timestamp: new Date().toISOString(),
              });
            }
          });
      }
    } catch (err) {
      console.error("Filter API error:", err);
      setError(getErrorMessage("chem", err));
      setResults([]);
    } finally {
      setLoading(false);
    }
  };

  // Change one filter option and re-apply when there is input
  const handleFilterChange = <K extends keyof FilterOptions>(name: K, value: FilterOptions[K]) => {
    const newFilterOptions = { ...filterOptions, [name]: value };
    setFilterOptions(newFilterOptions);
    if (smilesInput.trim() && !loading) {
      setTimeout(() => {
        applyFilters(null, null, newFilterOptions);
      }, 50);
    }
  };

  // Change AND/OR logic and re-apply when there is input
  const handleFilterOperatorChange = (newOperator: Operator) => {
    setFilterOperator(newOperator);
    if (smilesInput.trim()) {
      const updateAndApplyFilters = async () => {
        await new Promise((resolve) => setTimeout(resolve, 10));
        applyFilters(null, newOperator);
      };
      updateAndApplyFilters();
    }
  };

  const downloadResults = () => {
    if (results.length === 0) return;
    const text = results.join("\n");
    const blob = new Blob([text], { type: "text/plain;charset=utf-8" });
    const url = URL.createObjectURL(blob);
    const a = document.createElement("a");
    a.href = url;
    a.download = "filtered_molecules_results.txt";
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
    URL.revokeObjectURL(url);
  };

  // "CCO : T, F, T, T, F, T, 0.5, 3.2, 1.1" -> SMILES plus one value per column
  const parseFilterResult = (resultString: string) => {
    const resultArray: (boolean | string)[] = Array(9).fill("N/A");
    if (typeof resultString !== "string" || !resultString.includes(":")) {
      console.warn("Parsing unexpected result format:", resultString);
      return { smiles: resultString || "Invalid Data", results: resultArray };
    }

    const parts = resultString.split(" : ");
    const smiles = (parts[0] ?? "").trim();
    const filterValues = (parts[1] || "").split(",").map((item) => item.trim());

    // The six rule filters come first, only for the enabled ones; the three scores follow.
    let currentValueIndex = 0;
    RULES.forEach(({ key }, columnIndex) => {
      if (filterOptions[key]) {
        if (currentValueIndex < filterValues.length) {
          resultArray[columnIndex] = (filterValues[currentValueIndex] ?? "").toUpperCase() === "T";
          currentValueIndex++;
        } else {
          resultArray[columnIndex] = false;
        }
      } else {
        resultArray[columnIndex] = "disabled";
      }
    });
    [0, 1, 2].forEach((index) => {
      if (currentValueIndex < filterValues.length) {
        resultArray[index + 6] = filterValues[currentValueIndex] ?? "N/A";
        currentValueIndex++;
      }
    });

    console.log("Parsed filter result:", { smiles, values: resultArray });
    return { smiles, results: resultArray };
  };

  const hasResults = !loading && results.length > 0;

  useRevealOnResult("af-results", loading, ran, !!error);

  return (
    <div className="pc-bento">
      <div className="pc-col pc-t-input">
        <Cell id="af-input" title="Structure input" className="">
          <form onSubmit={applyFilters} className="pc-tool-form">
            <div className="grid gap-2">
              <label htmlFor="smilesInputArea" className="pc-lbl">
                Molecule list (SMILES)
              </label>
              <textarea
                id="smilesInputArea"
                className="pc-well resize-y"
                rows={6}
                value={smilesInput}
                onChange={(e) => setSmilesInput(e.target.value)}
                placeholder="Enter SMILES strings, one per line or separated by space/comma/semicolon..."
                spellCheck={false}
                autoComplete="off"
                required
                aria-required="true"
              />
              <p className="pc-note">
                One or more SMILES, separated by new lines, spaces, commas or semicolons.
              </p>
            </div>
            <div className="pc-row pc-actions">
              <Button type="submit" disabled={loading || !smilesInput.trim()}>
                {loading ? <Loader2 className="animate-spin" /> : <Filter />}
                {loading ? "Applying filters…" : "Apply filters"}
              </Button>
              <Button type="button" variant="outline" onClick={() => setSmilesInput(EXAMPLES)}>
                <FlaskConical />
                Load examples
              </Button>
            </div>
          </form>
        </Cell>

        <Cell id="af-about" title="About chemical filters" className="">
          <p className="pc-tool-note">
            Filters a list of molecules with medicinal-chemistry and drug-likeness rules. Detailed
            mode shows which properties failed and by how much, plus matched PAINS substructure
            families.
          </p>
          <dl className="pc-defs">
            <div>
              <dt>Match all (AND)</dt>
              <dd>A molecule must pass every selected filter.</dd>
            </div>
            <div>
              <dt>Match any (OR)</dt>
              <dd>A molecule must pass at least one selected filter.</dd>
            </div>
          </dl>
        </Cell>
      </div>

      <Cell id="af-opts" title="Options" className="pc-t-side">
        <div className="grid gap-6">
          <div className="pc-opt-group">
            <span className="pc-lbl">Rule-based filters</span>
            <div className="pc-opt-list">
              {RULES.map((rule) => (
                <Switch
                  key={rule.key}
                  checked={filterOptions[rule.key]}
                  onChange={(v) => handleFilterChange(rule.key, v)}
                  label={rule.label}
                  hint={rule.hint}
                />
              ))}
            </div>
          </div>

          <div className="pc-opt-group">
            <span className="pc-lbl">Score filters</span>
            <div className="pc-tool-form">
              {RANGES.map((range) => (
                <Select
                  key={range.key}
                  label={range.label}
                  value={filterOptions[range.key]}
                  options={range.options}
                  onChange={(v) => handleFilterChange(range.key, v)}
                />
              ))}
            </div>
          </div>

          <div className="pc-opt-group">
            <span className="pc-lbl" id="af-logic">
              Filter match logic
            </span>
            <Seg
              labelId="af-logic"
              full
              value={filterOperator}
              onChange={handleFilterOperatorChange}
              options={[
                { value: "AND", label: "Match all (AND)" },
                { value: "OR", label: "Match any (OR)" },
              ]}
            />
            <span className="pc-lbl">Result detail</span>
            <Switch
              checked={showDetailedView}
              onChange={setShowDetailedView}
              label="Show detailed violation information"
              hint="Property values and PAINS substructure details"
            />
          </div>
        </div>
      </Cell>

      <Cell
        id="af-results"
        title={hasResults ? "Filter results" : "Results"}
        className="pc-t-wide"
        meta={
          hasResults && (
            <div className="pc-row">
              <span className="pc-tag">
                {results.length} molecule{results.length === 1 ? "" : "s"}
              </span>
              <CopyButton text={results.join("\n")} label="Copy" />
              <Button variant="outline" onClick={downloadResults}>
                <FileDown />
                Download
              </Button>
            </div>
          )
        }
      >
        {loading && <ToolSkeleton variant="descriptors" />}
        {error && !loading && (
          <GlassErrorCard
            message={error}
            onRetry={() => {
              setError(null);
              document.getElementById("smilesInputArea")?.focus();
            }}
          />
        )}
        {!ran && !loading && !error && (
          <p className="pc-empty">Enter SMILES strings to screen them against the filters.</p>
        )}
        {ran && !loading && !error && results.length === 0 && (
          <p className="pc-empty">
            No molecules passed the selected filters; try other settings or switch AND/OR logic.
          </p>
        )}
        {hasResults && (
          <div className="overflow-x-auto rounded-xl border border-border">
            <table className="w-full min-w-[800px] border-collapse">
              <thead className="bg-muted">
                <tr>
                  <th
                    scope="col"
                    className={`${th} sticky left-0 z-10 border-r border-border bg-muted text-left`}
                  >
                    #
                  </th>
                  <th
                    scope="col"
                    className={`${th} sticky left-[calc(var(--spacing-12))] z-10 border-r border-border bg-muted text-left`}
                  >
                    SMILES
                  </th>
                  {filterColumns.map((col) => (
                    <th key={col.key} scope="col" className={`${th} text-center`}>
                      {col.label}
                    </th>
                  ))}
                </tr>
              </thead>
              <tbody className="divide-y divide-border">
                {results.map((resultString, index) => {
                  const { smiles, results: filterValues } = parseFilterResult(resultString);
                  return (
                    <tr key={index} className="transition-colors hover:bg-muted">
                      <td className="sticky left-0 z-10 border-r border-border bg-card px-3 py-2 font-mono text-sm text-muted-foreground tabular-nums">
                        {index + 1}
                      </td>
                      <td className="sticky left-[calc(var(--spacing-12))] z-10 border-r border-border bg-card px-4 py-2 font-mono text-sm whitespace-nowrap text-foreground">
                        {smiles}
                      </td>
                      {filterValues.map((value, filterIndex) => {
                        const isDisabled =
                          !!RULES[filterIndex] && !filterOptions[RULES[filterIndex]!.key];
                        return (
                          <td
                            key={filterIndex}
                            className={`px-3 py-2 text-center text-sm ${isDisabled ? "bg-muted" : ""}`}
                          >
                            {value === "disabled" || isDisabled ? (
                              <span className="text-muted-foreground">Disabled</span>
                            ) : typeof value === "boolean" ? (
                              value ? (
                                <Check
                                  className="inline-block size-5 text-primary"
                                  aria-label="Pass"
                                />
                              ) : (
                                <X
                                  className="inline-block size-5 text-madder-ink"
                                  aria-label="Fail"
                                />
                              )
                            ) : (
                              <span className="font-mono text-foreground tabular-nums">
                                {value === "N/A" ? "-" : value}
                              </span>
                            )}
                          </td>
                        );
                      })}
                    </tr>
                  );
                })}
              </tbody>
            </table>
          </div>
        )}
      </Cell>

      {!loading && showDetailedView && detailedResults && (
        <Cell
          id="af-detail"
          title="Detailed filter results"
          className="pc-t-wide"
          meta={
            <span className="pc-tag">
              {detailedResults.passing_molecules} of {detailedResults.total_molecules} passed (
              {detailedResults.filter_operator} logic)
            </span>
          }
        >
          <div className="pc-tool-form">
            {detailedResults.results.map((result, index) => (
              <div key={index} className="grid gap-3 rounded-xl border border-border p-4">
                <div className="pc-row">
                  <span className="font-mono text-sm break-all text-foreground">
                    {result.smiles}
                  </span>
                  <PassTag pass={result.overall_pass}>
                    {result.overall_pass ? "PASS" : "FAIL"}
                  </PassTag>
                </div>

                {!result.valid ? (
                  <p className="text-sm text-madder-ink">Error: {result.error}</p>
                ) : (
                  <div className="grid gap-3 md:grid-cols-2">
                    {Object.entries(result.filters).map(([filterName, filterResult]) => (
                      <div
                        key={filterName}
                        className="grid gap-2 rounded-xl border border-border p-3"
                      >
                        <div className="flex items-center justify-between gap-2">
                          <h3 className="text-sm font-medium text-foreground capitalize">
                            {filterName.replace(/_/g, " ")}
                          </h3>
                          <PassTag pass={filterResult.passes}>
                            {filterResult.passes ? "Pass" : "Fail"}
                          </PassTag>
                        </div>
                        <p className="pc-note">{filterResult.details}</p>

                        {filterResult.violations && filterResult.violations.length > 0 && (
                          <div className="text-xs text-madder-ink">
                            <p className="mb-1 font-medium">Violations</p>
                            <ul className="list-disc space-y-1 pl-4">
                              {filterResult.violations.map((violation, i) => (
                                <li key={i}>{violation}</li>
                              ))}
                            </ul>
                          </div>
                        )}

                        {filterResult.properties && (
                          <details className="text-xs text-muted-foreground">
                            <summary className="cursor-pointer hover:text-foreground">
                              View properties
                            </summary>
                            <dl className="mt-1 grid gap-0.5 border-l-2 border-border pl-2">
                              {Object.entries(filterResult.properties).map(([prop, value]) => (
                                <div key={prop} className="flex justify-between gap-2">
                                  <dt>{prop.replace(/_/g, " ")}</dt>
                                  <dd className="font-mono text-foreground tabular-nums">
                                    {String(value)}
                                  </dd>
                                </div>
                              ))}
                            </dl>
                          </details>
                        )}
                      </div>
                    ))}
                  </div>
                )}
              </div>
            ))}
          </div>
        </Cell>
      )}
    </div>
  );
};

export default AllFiltersView;
