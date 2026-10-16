import { useState, type FormEvent, type SyntheticEvent } from "react";
import { useRevealOnResult } from "@/hooks/useRevealOnResult";
import { FlaskConical, Loader2, RefreshCw, Search, Trash2 } from "lucide-react";
import SMILESInput from "../common/SMILESInput";
import MolFileUpload from "../common/MolFileUpload";
import MoleculeCard from "../common/MoleculeCard";
import SMILESDisplay from "../common/SMILESDisplay";
import { ToolSkeleton } from "@/components/feedback/ToolSkeleton";
import { GlassErrorCard } from "@/components/feedback/GlassErrorCard";
import { getErrorMessage } from "@/lib/error-messages";
import depictService from "../../services/depictService";
import {
  getSugarInfo,
  removeLinearSugars,
  removeCircularSugars,
  removeAllSugars,
  extractAglyconeAndSugars,
  getAglyconeAndSugarIndices,
} from "../../services/toolsService";
import { Button } from "@/components/ui/button";
import { Cell, Seg, Select, Switch } from "@/components/common/porcelain-ui";

// Default options configuration
const DEFAULT_OPTIONS = {
  // Basic options (apply to both sugar types)
  only_terminal: true,
  preservation_mode: 2,
  preservation_threshold: 5,
  mark_attach_points: false,

  // Circular sugar options
  gly_bond: false,
  oxygen_atoms: true,
  oxygen_atoms_threshold: 0.5,
  spiro_sugars: false,
  keto_sugars: false,

  // Linear sugar options
  linear_sugars_in_rings: false,
  linear_sugars_min_size: 4,
  linear_sugars_max_size: 7,
  linear_acidic_sugars: false,

  // Extract-specific options
  post_process_sugars: false,
  limit_post_process_by_size: false,
};

type SugarOptions = typeof DEFAULT_OPTIONS;
type OptionKey = keyof SugarOptions;
type Mode = "detect" | "remove" | "extract";
type SugarType = "circular" | "linear" | "both";
type SectionKey = "basic" | "circular" | "linear" | "extract" | "circular_detect" | "linear_detect";

type SugarResults =
  | {
      type: "detect";
      message: string;
      originalSmiles: string;
      hasCircular: boolean;
      hasLinear: boolean;
      circularIndices: number[][];
      linearIndices: number[][];
    }
  | { type: "remove"; removeType: SugarType; originalSmiles: string; resultSmiles: string }
  | {
      type: "extract";
      extractType: SugarType;
      originalSmiles: string;
      aglycone: string;
      sugars: string[];
    };

/** One option control: a switch, a number field or the preservation-mode select. */
interface OptionSpec {
  key: OptionKey;
  label: string;
  kind: "boolean" | "number" | "select";
  min?: number;
  max?: number;
  step?: number;
  hint: string;
}

const MODES: { value: Mode; label: string; hint: string }[] = [
  { value: "detect", label: "Detect", hint: "Identify sugar moieties." },
  { value: "remove", label: "Remove", hint: "Remove sugar moieties and keep the aglycone." },
  { value: "extract", label: "Extract", hint: "Extract the aglycone and each sugar moiety." },
];

const SUGAR_TYPES: { value: SugarType; label: string }[] = [
  { value: "circular", label: "Circular" },
  { value: "linear", label: "Linear" },
  { value: "both", label: "Both" },
];

const SUGAR_TYPE_LABEL: Record<SugarType, string> = {
  both: "Circular & linear sugars",
  circular: "Circular sugars only",
  linear: "Linear sugars only",
};

const PRESERVATION_MODES: [string, string][] = [
  ["1", "All: preserve all disconnected structures"],
  ["2", "Heavy atom count: remove structures below threshold"],
  ["3", "Molecular weight: remove structures below threshold"],
];

const circularOptions = (ratioHint: string): OptionSpec[] => [
  {
    key: "gly_bond",
    label: "Require O-glycosidic bonds",
    kind: "boolean",
    hint: "Consider only circular sugars with glycosidic bonds",
  },
  {
    key: "oxygen_atoms",
    label: "Require sufficient exocyclic oxygen atoms",
    kind: "boolean",
    hint: "Consider only circular sugars with enough exocyclic oxygen atoms",
  },
  {
    key: "oxygen_atoms_threshold",
    label: "Exocyclic oxygen atoms ratio",
    kind: "number",
    min: 0,
    max: 1,
    step: 0.1,
    hint: ratioHint,
  },
  {
    key: "spiro_sugars",
    label: "Include spiro sugars",
    kind: "boolean",
    hint: "Include rings that share one atom with another cycle",
  },
  {
    key: "keto_sugars",
    label: "Detect keto sugars",
    kind: "boolean",
    hint: "Detect circular sugars with keto groups",
  },
];

const LINEAR_OPTIONS: OptionSpec[] = [
  {
    key: "linear_sugars_in_rings",
    label: "Detect linear sugars in rings",
    kind: "boolean",
    hint: "Consider linear sugar patterns that are part of ring structures",
  },
  {
    key: "linear_sugars_min_size",
    label: "Minimum chain length",
    kind: "number",
    min: 0,
    max: 20,
    step: 1,
    hint: "Minimum carbon chain length to be considered a linear sugar",
  },
  {
    key: "linear_sugars_max_size",
    label: "Maximum chain length",
    kind: "number",
    min: 1,
    max: 20,
    step: 1,
    hint: "Maximum carbon chain length to be considered a linear sugar",
  },
  {
    key: "linear_acidic_sugars",
    label: "Detect acidic linear sugars",
    kind: "boolean",
    hint: "Include linear sugar moieties with acidic functional groups",
  },
];

const basicOptions = (mode: Mode): OptionSpec[] => [
  {
    key: "only_terminal",
    label: "Only terminal sugars",
    kind: "boolean",
    hint: "Only process sugars that can be removed without disconnecting the aglycone",
  },
  {
    key: "preservation_mode",
    label: "Preservation mode",
    kind: "select",
    hint: `Determines which disconnected structures to preserve after sugar ${mode === "remove" ? "removal" : "extraction"}`,
  },
  {
    key: "preservation_threshold",
    label: "Preservation threshold",
    kind: "number",
    min: 0,
    max: 100,
    step: 1,
    hint: "Threshold value for the selected preservation mode (e.g., minimum heavy atoms)",
  },
  {
    key: "mark_attach_points",
    label: "Mark attachment points",
    kind: "boolean",
    hint: `Mark the attachment points of ${mode === "remove" ? "removed" : "extracted"} sugars with a dummy atom`,
  },
];

const EXTRACT_OPTIONS: OptionSpec[] = [
  {
    key: "post_process_sugars",
    label: "Post-process extracted sugars",
    kind: "boolean",
    hint: "Split ether, ester, and peroxide bonds in extracted sugars",
  },
  {
    key: "limit_post_process_by_size",
    label: "Limit post-processing by size",
    kind: "boolean",
    hint: "Only split groups in post-processing that are above the preservation mode threshold",
  },
];

const BROKEN_IMAGE =
  "data:image/svg+xml;base64,PHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHZpZXdCb3g9IjAgMCAxMDAgMTAwIj48dGV4dCB4PSI1MCIgeT0iNTAiIGRvbWluYW50LWJhc2VsaW5lPSJtaWRkbGUiIHRleHQtYW5jaG9yPSJtaWRkbGUiIGZvbnQtc2l6ZT0iMTAiIGZpbGw9IiM4ODg4ODgiPkVycm9yIGxvYWRpbmcgc3RydWN0dXJlPC90ZXh0Pjwvc3ZnPg==";

/** Sugar Removal Utility: detect, remove or extract circular and linear sugar moieties. */
const SugarRemovalView = () => {
  const [smiles, setSmiles] = useState("");
  const [isLoading, setIsLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [uploadedFileName, setUploadedFileName] = useState<string | null>(null);
  const [operationMode, setOperationMode] = useState<Mode>("detect");
  // Which sugar types to remove / extract
  const [removeType, setRemoveType] = useState<SugarType>("circular");
  const [extractType, setExtractType] = useState<SugarType>("circular");
  const [options, setOptions] = useState<SugarOptions>(DEFAULT_OPTIONS);
  const [results, setResults] = useState<SugarResults | null>(null);
  useRevealOnResult("sr-results", isLoading, !!results, error);
  // Sugar highlighting toggle for detect mode
  const [highlightMode, setHighlightMode] = useState<SugarType>("circular");

  const resetAllOptions = () => setOptions(DEFAULT_OPTIONS);

  const resetOption = (key: OptionKey) => {
    setOptions((prev) => ({ ...prev, [key]: DEFAULT_OPTIONS[key] }));
  };

  const updateOption = (key: OptionKey, value: number | boolean) => {
    setOptions((prev) => ({ ...prev, [key]: value }));
  };

  const handleFileConversionSuccess = (
    convertedSmiles: string,
    _molblock: string,
    filename: string
  ) => {
    setSmiles(convertedSmiles);
    setUploadedFileName(filename);
    setError(null);
  };

  const handleFileConversionError = (errorMessage: string) => {
    setError(errorMessage);
    setUploadedFileName(null);
  };

  // Clear the uploaded-file indicator when the SMILES is edited by hand
  const handleSmilesChange = (newSmiles: string) => {
    setSmiles(newSmiles);
    if (uploadedFileName && newSmiles !== smiles) {
      setUploadedFileName(null);
    }
  };

  const handleExecute = async () => {
    const trimmedSmiles = smiles.trim();
    if (!trimmedSmiles) {
      setError("Please enter a SMILES string");
      setResults(null);
      return;
    }

    setIsLoading(true);
    setError(null);
    setResults(null);

    try {
      switch (operationMode) {
        case "detect": {
          const result = await getSugarInfo(trimmedSmiles, options);

          // Which sugars are present, from the message
          const hasCircular = result.toLowerCase().includes("circular");
          const hasLinear = result.toLowerCase().includes("linear");

          let circularIndices: number[][] = [];
          let linearIndices: number[][] = [];

          // Atom indices for highlighting
          if (hasCircular || hasLinear) {
            try {
              if (hasCircular) {
                const circularIndicesResult = await getAglyconeAndSugarIndices(trimmedSmiles, {
                  ...options,
                  extract_circular_sugars: true,
                  extract_linear_sugars: false,
                  only_terminal: false, // all sugars for highlighting
                });
                // Skip the first element (aglycone) and take the sugar indices
                if (Array.isArray(circularIndicesResult) && circularIndicesResult.length > 1) {
                  circularIndices = circularIndicesResult.slice(1);
                }
              }

              if (hasLinear) {
                const linearIndicesResult = await getAglyconeAndSugarIndices(trimmedSmiles, {
                  ...options,
                  extract_circular_sugars: false,
                  extract_linear_sugars: true,
                  only_terminal: false, // all sugars for highlighting
                });
                if (Array.isArray(linearIndicesResult) && linearIndicesResult.length > 1) {
                  linearIndices = linearIndicesResult.slice(1);
                }
              }
            } catch (extractError) {
              console.error("Could not get sugar indices for highlighting:", extractError);
            }
          }

          setResults({
            type: "detect",
            message: result,
            originalSmiles: trimmedSmiles,
            hasCircular,
            hasLinear,
            circularIndices,
            linearIndices,
          });
          break;
        }

        case "remove": {
          const result =
            removeType === "both"
              ? await removeAllSugars(trimmedSmiles, options)
              : removeType === "linear"
                ? await removeLinearSugars(trimmedSmiles, options)
                : await removeCircularSugars(trimmedSmiles, options);

          setResults({
            type: "remove",
            removeType,
            originalSmiles: trimmedSmiles,
            resultSmiles: result,
          });
          break;
        }

        case "extract": {
          const result = await extractAglyconeAndSugars(trimmedSmiles, {
            ...options,
            extract_circular_sugars: extractType === "circular" || extractType === "both",
            extract_linear_sugars: extractType === "linear" || extractType === "both",
          });
          if (Array.isArray(result) && result.length > 0) {
            setResults({
              type: "extract",
              extractType,
              originalSmiles: trimmedSmiles,
              aglycone: result[0] ?? "",
              sugars: result.slice(1),
            });
          } else {
            setError("No sugars were extracted from the molecule.");
          }
          break;
        }

        default:
          setError("Invalid operation mode selected.");
      }
    } catch (err) {
      console.error(`Error during ${operationMode} operation:`, err);
      setError(getErrorMessage("tools", err));
    } finally {
      setIsLoading(false);
    }
  };

  const submit = (e: FormEvent) => {
    e.preventDefault();
    handleExecute();
  };

  const renderOption = (spec: OptionSpec) => {
    const { key, label, kind, min, max, step, hint } = spec;
    const value = options[key];
    const isModified = value !== DEFAULT_OPTIONS[key];
    const fieldId = `sugar-option-${key}`;

    let control;
    if (kind === "boolean") {
      control = (
        <Switch
          checked={Boolean(value)}
          onChange={(v) => updateOption(key, v)}
          label={label}
          hint={hint}
        />
      );
    } else if (kind === "select") {
      control = (
        <div className="grid gap-1">
          <Select
            label={label}
            value={String(value)}
            options={PRESERVATION_MODES}
            onChange={(v) => updateOption(key, parseInt(v))}
          />
          <p className="pc-note">{hint}</p>
        </div>
      );
    } else {
      control = (
        <div className="grid gap-1">
          <label className="pc-field" htmlFor={fieldId}>
            <span>{label}</span>
            <input
              id={fieldId}
              type="number"
              className="pc-input"
              value={Number(value)}
              onChange={(e) => updateOption(key, parseFloat(e.target.value))}
              min={min}
              max={max}
              step={step || 1}
            />
          </label>
          <p className="pc-note">{hint}</p>
        </div>
      );
    }

    return (
      <div key={key} className="flex items-start gap-1">
        <div className="min-w-0 flex-1">{control}</div>
        <div className="flex items-center gap-1 pt-1">
          {isModified && <span className="text-xs text-primary">Modified</span>}
          <Button
            variant="ghost"
            size="icon-sm"
            type="button"
            onClick={() => resetOption(key)}
            disabled={!isModified}
            title="Reset to default"
            aria-label={`Reset ${label} to default`}
          >
            <RefreshCw />
          </Button>
        </div>
      </div>
    );
  };

  // every option stays visible: no collapsing sections
  const renderSection = (sectionKey: SectionKey, title: string, specs: OptionSpec[]) => {
    const headId = `sugar-section-${sectionKey}`;
    return (
      <div className="pc-opt-group" key={sectionKey} role="group" aria-labelledby={headId}>
        <span className="pc-lbl" id={headId}>
          {title}
        </span>
        <div className="pc-opt-list">{specs.map(renderOption)}</div>
      </div>
    );
  };

  const sugarType = operationMode === "remove" ? removeType : extractType;
  const showCircular = sugarType === "circular" || sugarType === "both";
  const showLinear = sugarType === "linear" || sugarType === "both";
  const ModeIcon =
    operationMode === "detect" ? Search : operationMode === "remove" ? Trash2 : FlaskConical;
  const actionLabel =
    operationMode === "detect"
      ? "Detect sugars"
      : operationMode === "remove"
        ? "Remove sugars"
        : "Extract aglycone and sugars";

  const done = !!results && !isLoading && !error;
  const resultTitle = !done
    ? "Results"
    : results.type === "detect"
      ? "Sugar detection"
      : results.type === "remove"
        ? "Sugar removal"
        : "Aglycone and sugars";
  const resultMeta =
    done && results.type === "remove" ? (
      <span className="pc-tag">{SUGAR_TYPE_LABEL[results.removeType]}</span>
    ) : done && results.type === "extract" ? (
      <span className="pc-tag">{SUGAR_TYPE_LABEL[results.extractType]}</span>
    ) : undefined;

  const renderDetect = (r: Extract<SugarResults, { type: "detect" }>) => {
    // Which atom indices to highlight, from the highlight mode
    const atomIndicesToHighlight =
      r.hasCircular && r.hasLinear
        ? highlightMode === "circular"
          ? r.circularIndices
          : highlightMode === "linear"
            ? r.linearIndices
            : [...r.circularIndices, ...r.linearIndices]
        : r.hasCircular
          ? r.circularIndices
          : r.hasLinear
            ? r.linearIndices
            : [];

    const depictionUrl = depictService.get2DDepictionUrl(r.originalSmiles, {
      toolkit: "cdk",
      width: 600,
      height: 400,
      atomIds: atomIndicesToHighlight,
      CIP: false,
      unicolor: false,
    });

    const description =
      r.hasCircular && r.hasLinear
        ? highlightMode === "circular"
          ? "Circular sugars highlighted"
          : highlightMode === "linear"
            ? "Linear sugars highlighted"
            : "Both circular and linear sugars highlighted"
        : r.hasCircular
          ? "Circular sugars highlighted"
          : r.hasLinear
            ? "Linear sugars highlighted"
            : r.message;

    return (
      <div className="pc-tool-form">
        <p className="pc-tool-note">{r.message}</p>
        {r.hasCircular && r.hasLinear && (
          <div className="pc-opt-group">
            <span className="pc-lbl" id="sugar-highlight">
              Highlight
            </span>
            <Seg
              labelId="sugar-highlight"
              value={highlightMode}
              onChange={setHighlightMode}
              options={[
                { value: "circular", label: "Circular only" },
                { value: "linear", label: "Linear only" },
                { value: "both", label: "Both" },
              ]}
            />
          </div>
        )}
        <SMILESDisplay smiles={r.originalSmiles} label="SMILES" />
        <figure className="m-0 grid gap-2">
          <figcaption className="pc-lbl">{description}</figcaption>
          <div className="flex min-h-[400px] items-center justify-center rounded-xl border border-border bg-card p-2">
            <img
              src={depictionUrl}
              alt="Molecule with highlighted sugars"
              className="max-h-full max-w-full object-contain"
              onError={(e: SyntheticEvent<HTMLImageElement>) => {
                console.error("Error loading depiction image");
                const img = e.currentTarget;
                img.onerror = null;
                img.src = BROKEN_IMAGE;
              }}
            />
          </div>
          <p className="pc-note">
            Sugar moieties are highlighted using precise atom indices from the Sugar Removal Utility
            analysis.
          </p>
        </figure>
      </div>
    );
  };

  return (
    <div className="pc-bento">
      <div className="pc-col pc-t-input">
        <Cell id="sr-input" title="Structure input" className="">
          <form onSubmit={submit} className="pc-tool-form">
            <SMILESInput
              value={smiles}
              onChange={handleSmilesChange}
              label="SMILES String"
              placeholder="Enter SMILES (e.g., C=CC1C(C[C@@H]2NCCC3=C2NC2=CC=CC=C32)C(C(=O)O)=CO[C@H]1O[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O)"
              required
            />

            <div className="grid gap-2">
              <MolFileUpload
                onConversionSuccess={handleFileConversionSuccess}
                onConversionError={handleFileConversionError}
                toolkit="cdk"
                allowMultipleMolecules={false}
              />
              {uploadedFileName && (
                <p className="pc-note">
                  Loaded from <span className="pc-tag">{uploadedFileName}</span>
                </p>
              )}
            </div>

            <div className="pc-opt-group">
              <span className="pc-lbl" id="sr-mode">
                Operation
              </span>
              <Seg
                labelId="sr-mode"
                full
                value={operationMode}
                onChange={setOperationMode}
                options={MODES}
              />
              <p className="pc-note mt-2">{MODES.find((m) => m.value === operationMode)?.hint}</p>
            </div>

            {operationMode !== "detect" && (
              <div className="pc-opt-group">
                <span className="pc-lbl" id="sr-types">
                  {operationMode === "remove" ? "Sugar types to remove" : "Sugar types to extract"}
                </span>
                <Seg
                  labelId="sr-types"
                  full
                  value={sugarType}
                  onChange={operationMode === "remove" ? setRemoveType : setExtractType}
                  options={SUGAR_TYPES}
                />
              </div>
            )}

            <div className="pc-row pc-actions">
              <Button type="submit" disabled={!smiles.trim() || isLoading}>
                {isLoading ? <Loader2 className="animate-spin" /> : <ModeIcon />}
                {isLoading ? "Processing…" : actionLabel}
              </Button>
            </div>
          </form>
        </Cell>

        <Cell id="sr-results" title={resultTitle} meta={resultMeta} className="">
          {isLoading && <ToolSkeleton variant="molecule" />}
          {error && !isLoading && (
            <GlassErrorCard
              message={error}
              onRetry={() => {
                setError(null);
                document.getElementById("smiles-input")?.focus();
              }}
            />
          )}
          {!results && !isLoading && !error && (
            <p className="pc-empty">
              Enter a SMILES string, choose an operation and run it to see the sugars.
            </p>
          )}

          {done && results.type === "detect" && renderDetect(results)}

          {done && results.type === "remove" && (
            <div className="pc-result-grid">
              <div className="grid content-start gap-3">
                <MoleculeCard
                  smiles={results.originalSmiles}
                  title="Before Removal"
                  description="Original molecule with sugars"
                  showActions={false}
                />
                <SMILESDisplay smiles={results.originalSmiles} label="Original SMILES" />
              </div>
              <div className="grid content-start gap-3">
                <MoleculeCard
                  smiles={results.resultSmiles}
                  title="After Removal"
                  description={`${
                    results.removeType === "both"
                      ? "Circular & linear sugars"
                      : results.removeType === "circular"
                        ? "Circular sugars"
                        : "Linear sugars"
                  } removed`}
                  showActions={false}
                />
                <SMILESDisplay smiles={results.resultSmiles} label="Aglycone SMILES" />
              </div>
            </div>
          )}

          {done && results.type === "extract" && (
            <div className="pc-tool-form">
              <div className="pc-result-grid">
                <div className="grid content-start gap-3">
                  <MoleculeCard
                    smiles={results.aglycone}
                    title="Aglycone"
                    description="Core structure after sugar removal"
                    showActions
                  />
                  <SMILESDisplay smiles={results.aglycone} label="Aglycone SMILES" />
                </div>
                {results.sugars.map((sugarSmiles, index) => (
                  <div key={index} className="grid content-start gap-3">
                    <MoleculeCard
                      smiles={sugarSmiles}
                      title={`Sugar ${index + 1}`}
                      description="Extracted sugar moiety"
                      showActions
                      size="md"
                    />
                    <SMILESDisplay smiles={sugarSmiles} label={`Sugar ${index + 1} SMILES`} />
                  </div>
                ))}
              </div>
              {results.sugars.length === 0 && (
                <p className="pc-empty">No sugar moieties were extracted from the molecule.</p>
              )}
            </div>
          )}
        </Cell>
      </div>

      <div className="pc-col pc-t-side">
        <Cell
          id="sr-opts"
          title="Options"
          className=""
          meta={
            <Button variant="ghost" size="sm" type="button" onClick={resetAllOptions}>
              <RefreshCw />
              Reset all to default
            </Button>
          }
        >
          <div className="grid items-start gap-4">
            {operationMode === "detect" ? (
              <>
                {renderSection(
                  "circular_detect",
                  "Circular sugar detection options",
                  circularOptions(
                    "Minimum ratio of exocyclic oxygen atoms to ring atoms (0.5 = 50%)"
                  )
                )}
                {renderSection("linear_detect", "Linear sugar detection options", LINEAR_OPTIONS)}
              </>
            ) : (
              <>
                {renderSection("basic", "Basic options", basicOptions(operationMode))}
                {showCircular &&
                  renderSection(
                    "circular",
                    "Circular sugar detection options",
                    circularOptions("0.5 = a pyranose needs at least 3 hydroxy groups")
                  )}
                {showLinear &&
                  renderSection("linear", "Linear sugar detection options", LINEAR_OPTIONS)}
                {operationMode === "extract" &&
                  renderSection("extract", "Advanced extraction options", EXTRACT_OPTIONS)}
              </>
            )}
          </div>
        </Cell>

        <Cell id="sr-about" title="About sugar removal" className="">
          <p className="pc-tool-note">
            Detects, removes or extracts circular and linear sugar moieties with the Sugar Removal
            Utility. Removal keeps the aglycone; extraction returns the aglycone and each sugar
            separately.
          </p>
          <dl className="pc-defs">
            <div>
              <dt>Reference</dt>
              <dd>
                Schaub, J., Zielesny, A., Steinbeck, C., Sorokina, M. Too sweet: cheminformatics for
                deglycosylation in natural products. J Cheminform 12, 67 (2020).{" "}
                <a
                  href="https://doi.org/10.1186/s13321-020-00467-y"
                  target="_blank"
                  rel="noopener noreferrer"
                  className="text-link font-semibold"
                >
                  View publication
                </a>
              </dd>
            </div>
            <div>
              <dt>Batch processing</dt>
              <dd>
                For many molecules, use{" "}
                <a
                  href="https://felixbaensch.github.io/MORTAR"
                  target="_blank"
                  rel="noopener noreferrer"
                  className="text-link font-semibold"
                >
                  MORTAR
                </a>{" "}
                (MOlecule fragmenTation fRamework) or the{" "}
                <a
                  href="https://github.com/JonasSchaub/SugarRemoval"
                  target="_blank"
                  rel="noopener noreferrer"
                  className="text-link font-semibold"
                >
                  SugarRemoval CMD app
                </a>
                .
              </dd>
            </div>
          </dl>
        </Cell>
      </div>
    </div>
  );
};

export default SugarRemovalView;
