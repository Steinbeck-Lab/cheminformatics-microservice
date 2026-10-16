import {
  useState,
  useEffect,
  useRef,
  useCallback,
  type ChangeEvent,
  type FormEvent,
  type ReactNode,
} from "react";
import { useKetcherEditor } from "@/hooks/useKetcherEditor";
import {
  DEFAULT_VERSION,
  VERSION_OPTIONS,
  loadInchiModule,
  convertMolfileToInchi,
  generateInchiKey,
  convertInchiToMolfile,
  convertAuxinfoToMolfile,
} from "@/utils/inchiUtils";
import {
  Check,
  ChevronDown,
  ChevronUp,
  ClipboardCopy,
  FileText,
  Loader2,
  RefreshCw,
  Upload,
  X,
} from "lucide-react";
import { GlassErrorCard } from "@/components/feedback/GlassErrorCard";
import { Button } from "@/components/ui/button";
import { Dialog, DialogContent, DialogTitle } from "@/components/ui/dialog";
import { Input } from "@/components/ui/input";
import { Cell, Seg, Select, Switch } from "@/components/common/porcelain-ui";
import { useRevealOnResult } from "@/hooks/useRevealOnResult";
import { errorText } from "@/lib/error-messages";

type InputType = "structure" | "smiles" | "molfile" | "inchi";

const ENHANCED_STEREO = "Dev-EnhancedStereo";
const POLYMER = "Polymer";
// Builds whose options include Molecular Inorganics (on by default), per the IUPAC demo.
const MOIN = new Set(["Latest-MoIn", "Dev-NoMetalH", "Dev-ExplicitZeroValence"]);

/** An empty Ketcher canvas still returns a molfile: header plus a counts line with 0 atoms. */
function hasAtoms(molfile: string | undefined): boolean {
  if (!molfile?.trim()) return false;
  const v2000 = molfile.split(/\r?\n/).find((l) => /V2000\s*$/.test(l));
  if (v2000) return parseInt(v2000.slice(0, 3), 10) > 0;
  const v3000 = molfile.match(/M {2}V30 COUNTS (\d+)/);
  return v3000 ? parseInt(v3000[1]!, 10) > 0 : true;
}

const STEREO_MODES: [string, string][] = [
  ["SAbs", "Absolute"],
  ["SRel", "Relative"],
  ["SRac", "Racemic"],
  ["SUCF", "From chiral flag"],
];

const INPUT_TYPES: { value: InputType; label: string }[] = [
  { value: "structure", label: "Draw" },
  { value: "smiles", label: "SMILES" },
  { value: "molfile", label: "Mol/Aux" },
  { value: "inchi", label: "InChI" },
];

/** InChI options: builds the "-Option" string the InChI library expects. */
const InChIOptions = ({
  onChange,
  inchiVersion,
  setInchiVersion,
  children,
}: {
  onChange: (options: string) => void;
  inchiVersion: string;
  setInchiVersion: (version: string) => void;
  children?: ReactNode;
}) => {
  const [mobileH, setMobileH] = useState(true);
  const [includeStereo, setIncludeStereo] = useState(true);
  const [stereoMode, setStereoMode] = useState("SAbs");
  const [recMet, setRecMet] = useState(false);
  const [ket, setKet] = useState(false);
  const [t15, setT15] = useState(false);
  const [polymers, setPolymers] = useState(false);
  const [NPZz, setNPZz] = useState(false);
  const [molecularInorganics, setMolecularInorganics] = useState(MOIN.has(inchiVersion));
  const [enhancedStereo, setEnhancedStereo] = useState(false);
  const [showTautomerism, setShowTautomerism] = useState(true);

  // Additional stereo options
  const [SUU, setSUU] = useState(false);
  const [SLUUD, setSLUUD] = useState(false);
  const [NEWPSOFF, setNEWPSOFF] = useState(false);

  // Additional tautomer options
  const [PT_06_00, setPT_06_00] = useState(false); // 1,3 heteroatom H-shift
  const [PT_13_00, setPT_13_00] = useState(false); // keten-inol exchange
  const [PT_16_00, setPT_16_00] = useState(false); // nitroso/oxime
  const [PT_18_00, setPT_18_00] = useState(false); // cyanic/iso-cyanic acids
  const [PT_22_00, setPT_22_00] = useState(false); // imine/imine
  const [PT_39_00, setPT_39_00] = useState(false); // nitrone/azoxy or Behrend rearrangement

  // Polymer detail options
  const [NoEdits, setNoEdits] = useState(false);
  const [FoldCRU, setFoldCRU] = useState(false);
  const [NoFrameShift, setNoFrameShift] = useState(false);

  // Generate and propagate the options string when settings change
  useEffect(() => {
    const options = [];

    if (!mobileH) options.push("FixedH");
    if (!includeStereo) options.push("SNon");

    if (includeStereo) {
      if (stereoMode === "SRel") options.push("SRel");
      if (stereoMode === "SRac") options.push("SRac");
      if (stereoMode === "SUCF") options.push("SUCF");

      if (SUU) options.push("SUU");
      if (SLUUD) options.push("SLUUD");
      if (NEWPSOFF) options.push("NEWPSOFF");
    }

    if (recMet) options.push("RecMet");
    if (ket) options.push("KET");
    if (t15) options.push("15T");

    if (PT_06_00) options.push("PT_06_00");
    if (PT_13_00) options.push("PT_13_00");
    if (PT_16_00) options.push("PT_16_00");
    if (PT_18_00) options.push("PT_18_00");
    if (PT_22_00) options.push("PT_22_00");
    if (PT_39_00) options.push("PT_39_00");

    if (polymers) {
      options.push("Polymers");
      if (NoEdits) options.push("NoEdits");
      if (FoldCRU) options.push("FoldCRU");
      if (NoFrameShift) options.push("NoFrameShift");
    }

    if (NPZz) options.push("NPZz");

    // MolecularInorganics exists in the MoIn-based builds only
    if (molecularInorganics && MOIN.has(inchiVersion)) options.push("MolecularInorganics");

    // EnhancedStereochemistry only exists in the enhanced-stereo build
    if (enhancedStereo && inchiVersion === ENHANCED_STEREO) options.push("EnhancedStereochemistry");

    // The options string format the API expects
    const optionsString = options.map((opt) => `-${opt}`).join(" ");
    onChange(optionsString);
  }, [
    mobileH,
    includeStereo,
    stereoMode,
    recMet,
    ket,
    t15,
    polymers,
    NPZz,
    molecularInorganics,
    enhancedStereo,
    SUU,
    SLUUD,
    NEWPSOFF,
    PT_06_00,
    PT_13_00,
    PT_16_00,
    PT_18_00,
    PT_22_00,
    PT_39_00,
    NoEdits,
    FoldCRU,
    NoFrameShift,
    inchiVersion,
    onChange,
  ]);

  // Version defaults from the IUPAC demo: the polymer build starts with polymer handling and
  // NPZz on; the MoIn build starts with NPZz on.
  useEffect(() => {
    if (inchiVersion === POLYMER) {
      setPolymers(true);
      setNPZz(true);
    } else if (inchiVersion === "Latest-MoIn") {
      setNPZz(true);
    }
    if (inchiVersion !== ENHANCED_STEREO) setEnhancedStereo(false);
  }, [inchiVersion]);

  // MolecularInorganics follows the InChI version
  useEffect(() => {
    if (MOIN.has(inchiVersion)) {
      setMolecularInorganics(true);
    } else if (molecularInorganics) {
      setMolecularInorganics(false);
    }
  }, [inchiVersion, molecularInorganics]);

  // Polymer details reset when polymers are turned off
  useEffect(() => {
    if (!polymers) {
      setNoEdits(false);
      setFoldCRU(false);
      setNoFrameShift(false);
    }
  }, [polymers]);

  const resetOptions = () => {
    setMobileH(true);
    setIncludeStereo(true);
    setStereoMode("SAbs");
    setRecMet(false);
    setKet(false);
    setT15(false);
    setPolymers(inchiVersion === POLYMER);
    setNPZz(inchiVersion === POLYMER || inchiVersion === "Latest-MoIn");
    setMolecularInorganics(MOIN.has(inchiVersion));
    setEnhancedStereo(false);
    setSUU(false);
    setSLUUD(false);
    setNEWPSOFF(false);
    setPT_06_00(false);
    setPT_13_00(false);
    setPT_16_00(false);
    setPT_18_00(false);
    setPT_22_00(false);
    setPT_39_00(false);
    setNoEdits(false);
    setFoldCRU(false);
    setNoFrameShift(false);
  };

  const isLatest = inchiVersion !== "1.06";

  return (
    <Cell
      id="inchi-opts"
      title="Options"
      className="pc-t-side"
      meta={
        <Button variant="ghost" size="sm" onClick={resetOptions}>
          <RefreshCw />
          Reset
        </Button>
      }
    >
      <div className="pc-tool-form">
        {children}
        <Select
          label="InChI version"
          value={inchiVersion}
          options={VERSION_OPTIONS}
          onChange={setInchiVersion}
        />

        <div className="pc-opt-list">
          <Switch
            checked={mobileH}
            onChange={setMobileH}
            label="Mobile H perception"
            hint="Off uses the FixedH option, which perceives fixed H positions"
          />
          <Switch
            checked={includeStereo}
            onChange={setIncludeStereo}
            label="Include stereo"
            hint="Off uses the SNon option, which ignores stereochemistry"
          />
          {includeStereo && (
            <div className="grid gap-1 pl-4">
              <Select
                label="Stereo"
                value={stereoMode}
                options={STEREO_MODES}
                onChange={setStereoMode}
              />
              <Switch
                checked={SUU}
                onChange={setSUU}
                label="Always include omitted/undefined stereo (SUU)"
              />
              <Switch
                checked={SLUUD}
                onChange={setSLUUD}
                label="Different marks for unknown/undefined stereo (SLUUD)"
              />
              <Switch
                checked={NEWPSOFF}
                onChange={setNEWPSOFF}
                label="Both ends of wedge point to stereocenters (NEWPSOFF)"
                hint="Old stereocenter recognition: the sharp end points to a stereocenter"
              />
            </div>
          )}
          <Switch
            checked={recMet}
            onChange={setRecMet}
            label="Include bonds to metal"
            hint="Include metal bonds, which are otherwise disconnected"
          />
        </div>

        <div className="pc-opt-group">
          <div className="flex items-center justify-between gap-2">
            <span className="pc-lbl" id="inchi-taut">
              Tautomerism
            </span>
            <Button
              variant="ghost"
              size="sm"
              onClick={() => setShowTautomerism(!showTautomerism)}
              aria-expanded={showTautomerism}
              aria-controls="inchi-taut-list"
            >
              {showTautomerism ? <ChevronUp /> : <ChevronDown />}
              {showTautomerism ? "Hide" : "Show"}
            </Button>
          </div>
          {showTautomerism && (
            <div id="inchi-taut-list" className="pc-opt-list" aria-labelledby="inchi-taut">
              <Switch checked={ket} onChange={setKet} label="Keto-enol tautomerism (KET)" />
              <Switch checked={t15} onChange={setT15} label="1,5-tautomerism (15T)" />
              {/* Further tautomer options exist from InChI 1.07 */}
              {isLatest && (
                <>
                  <Switch
                    checked={PT_06_00}
                    onChange={setPT_06_00}
                    label="1,3 heteroatom H-shift (PT_06_00)"
                  />
                  <Switch
                    checked={PT_13_00}
                    onChange={setPT_13_00}
                    label="Keten-inol exchange (PT_13_00)"
                  />
                  <Switch
                    checked={PT_16_00}
                    onChange={setPT_16_00}
                    label="Nitroso/oxime (PT_16_00)"
                  />
                  <Switch
                    checked={PT_18_00}
                    onChange={setPT_18_00}
                    label="Cyanic/iso-cyanic acids (PT_18_00)"
                  />
                  <Switch
                    checked={PT_22_00}
                    onChange={setPT_22_00}
                    label="Imine/imine (PT_22_00)"
                  />
                  <Switch
                    checked={PT_39_00}
                    onChange={setPT_39_00}
                    label="Nitrone/azoxy or Behrend rearrangement (PT_39_00)"
                  />
                </>
              )}
            </div>
          )}
        </div>

        <div className="pc-opt-list">
          <Switch
            checked={polymers}
            onChange={setPolymers}
            label="Treat polymers"
            hint="Process polymers with the InChI polymer extension"
          />
          {polymers && (
            <div className="grid gap-1 pl-4">
              <Switch
                checked={NoEdits}
                onChange={setNoEdits}
                label="No pre-edits of original polymer structure"
              />
              <Switch
                checked={FoldCRU}
                onChange={setFoldCRU}
                label="Enable CRU folding"
                hint="Constitutional repeating unit folding"
              />
              <Switch
                checked={NoFrameShift}
                onChange={setNoFrameShift}
                label="Disable CRU frame shift"
              />
            </div>
          )}
          <Switch
            checked={NPZz}
            onChange={setNPZz}
            label="Allow non-polymer Zz pseudoatoms"
            hint="Allow Zz or * (star) pseudoatoms outside the polymer context"
          />
          {inchiVersion === ENHANCED_STEREO && (
            <Switch
              checked={enhancedStereo}
              onChange={setEnhancedStereo}
              label="Enhanced stereochemistry"
              hint="Experimental. Reads the AND/OR stereo groups of a V3000 molfile"
            />
          )}
          {MOIN.has(inchiVersion) && (
            <Switch
              checked={molecularInorganics}
              onChange={setMolecularInorganics}
              label="Molecular inorganics"
              hint="Experimental option of the dev builds"
            />
          )}
        </div>
      </div>
    </Cell>
  );
};

/** InChI: draw or load a structure in Ketcher and convert it with the IUPAC InChI library. */
const InChIView = () => {
  const [inchi, setInchi] = useState("");
  const [auxInfo, setAuxInfo] = useState("");
  const [inchiKey, setInchiKey] = useState("");
  const [inputInchi, setInputInchi] = useState("");
  const [inputSmiles, setInputSmiles] = useState("");
  const [molfileContent, setMolfileContent] = useState("");
  const [logMessage, setLogMessage] = useState("");
  const [options, setOptions] = useState("");
  const [inchiVersion, setInchiVersion] = useState(DEFAULT_VERSION);
  const [activeInputType, setActiveInputType] = useState<InputType>("structure");

  const [isLoading, setIsLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [copySuccess, setCopySuccess] = useState(false);
  const [showCopyModal, setShowCopyModal] = useState(false);
  const [copyModalText, setCopyModalText] = useState("");
  const [inchiModuleLoaded, setInchiModuleLoaded] = useState(false);
  const [autoGenerate, setAutoGenerate] = useState(false);

  const copyTextRef = useRef<HTMLInputElement>(null);
  const fileInputRef = useRef<HTMLInputElement>(null);
  const previousOptionsRef = useRef("");

  const { ketcherFrame, isEditorReady, executeCommand, resetEditor } = useKetcherEditor();

  // Enhanced stereochemistry (AND/OR groups) only survives in V3000; other builds read V2000.
  const getMolfile = useCallback(
    (version: string) =>
      version === ENHANCED_STEREO
        ? executeCommand("getMolfile", ["v3000"])
        : executeCommand("getMolfile"),
    [executeCommand]
  );

  // Clear a stale error when the input, version or input method changes.
  // (`error` used to be a dependency too, which cleared every error the moment it was set.)
  useEffect(() => {
    setError(null);
  }, [inputInchi, inputSmiles, inchiVersion, activeInputType]);

  // Auto-select text in copy modal when it appears
  useEffect(() => {
    if (showCopyModal && copyTextRef.current) {
      copyTextRef.current.select();
    }
  }, [showCopyModal]);

  // Generate InChIKey from InChI
  const generateInChIKeyFromInChI = useCallback(
    async (inchiString: string) => {
      if (!inchiString) {
        console.warn("Cannot generate InChIKey: No InChI provided");
        return;
      }

      if (!inchiString.startsWith("InChI=")) {
        console.warn("Invalid InChI format for key generation:", inchiString);
        return;
      }

      try {
        const result = await generateInchiKey(inchiString, inchiVersion);

        if (!result) {
          console.error("No result returned from InChIKey generation");
          return;
        }

        if (result.return_code === -1) {
          console.error("Failed to generate InChIKey:", result.message);
          return;
        }

        if (!result.inchikey) {
          console.error("InChIKey not found in result");
          return;
        }

        setInchiKey(result.inchikey);

        const message = result.message;
        if (message) {
          setLogMessage((prevLog) => (prevLog ? `${prevLog}\n${message}` : message));
        }

        return result.inchikey;
      } catch (err) {
        console.error("Error generating InChIKey:", err);
        // A secondary operation: no error message
        return null;
      }
    },
    [inchiVersion]
  );

  // Generate InChI from a molfile
  const generateInChIFromMolfile = useCallback(
    async (molfile: string) => {
      try {
        if (!molfile || molfile.trim() === "") {
          throw new Error("No molecule data provided");
        }

        if (!inchiModuleLoaded) {
          throw new Error("InChI module not loaded. Please wait or refresh the page.");
        }

        const result = await convertMolfileToInchi(molfile, options, inchiVersion);

        if (!result) {
          throw new Error("No result returned from conversion function");
        }

        if (result.return_code === -1) {
          throw new Error(result.message || "Failed to generate InChI");
        }

        setInchi(result.inchi || "");
        setAuxInfo(result.auxinfo || "");

        setLogMessage(
          `InChI generated with options: ${options || "none"}\n${result.message || ""}`
        );

        if (result.inchi) {
          await generateInChIKeyFromInChI(result.inchi);
        }

        return result;
      } catch (err) {
        console.error("Failed to generate InChI from molfile:", err);
        setLogMessage(`Failed to generate InChI: ${errorText(err)}`);
        setError(`Failed to generate InChI: ${errorText(err)}`);
        throw err;
      } finally {
        setIsLoading(false);
      }
    },
    [inchiModuleLoaded, options, inchiVersion, generateInChIKeyFromInChI]
  );

  // Generate InChI from the structure in Ketcher
  const generateInChI = useCallback(async () => {
    if (!isEditorReady) {
      setError("Editor not ready. Please try again in a moment.");
      return;
    }

    if (!inchiModuleLoaded) {
      setError("InChI module not loaded. Please wait or refresh the page.");
      return;
    }

    setIsLoading(true);
    setError(null);
    setLogMessage("Generating InChI...");
    setInchi("");
    setAuxInfo("");
    setInchiKey("");

    try {
      const molfile = (await getMolfile(inchiVersion)) as string;

      if (!hasAtoms(molfile)) {
        setError("No structure drawn. Draw a molecule, or load a SMILES first.");
        setIsLoading(false);
        return;
      }

      // Reactions are not supported by InChI
      if (molfile.includes("$RXN")) {
        setError("Reactions are not supported. Please draw a single molecule.");
        setIsLoading(false);
        return;
      }

      setMolfileContent(molfile);
      await generateInChIFromMolfile(molfile);
    } catch (err) {
      console.error("Failed to generate InChI:", err);
      setLogMessage(`Error: ${errorText(err)}`);
      setError(`Failed to generate InChI: ${errorText(err)}`);
      setIsLoading(false);
    }
  }, [isEditorReady, inchiModuleLoaded, getMolfile, inchiVersion, generateInChIFromMolfile]);

  // Log option changes and auto-generate if enabled
  useEffect(() => {
    if (options !== previousOptionsRef.current) {
      if (previousOptionsRef.current && options) {
        const optionChangeLog = `InChI options changed: ${options || "none"}`;
        setLogMessage((prevLog) => (prevLog ? `${optionChangeLog}\n${prevLog}` : optionChangeLog));

        // Auto-generate if enabled and there is a structure
        if (autoGenerate && isEditorReady && molfileContent) {
          // Slight delay to let the UI update first
          setTimeout(() => {
            generateInChI();
          }, 100);
        }
      }

      previousOptionsRef.current = options;
    }
  }, [options, autoGenerate, isEditorReady, molfileContent, generateInChI]);

  // Copy with fallbacks: Clipboard API, execCommand, then a dialog to copy by hand
  const copyToClipboard = async (text: string | null = null, type = "inchi") => {
    let textToCopy;

    switch (type) {
      case "inchi":
        textToCopy = text || inchi;
        break;
      case "key":
        textToCopy = text || inchiKey;
        break;
      case "auxinfo":
        textToCopy = text || auxInfo;
        break;
      default:
        textToCopy = text || inchi;
    }

    if (!textToCopy) {
      setError(`No ${type.toUpperCase()} to copy. Generate ${type.toUpperCase()} first.`);
      return;
    }

    try {
      if (navigator.clipboard && navigator.clipboard.writeText) {
        await navigator.clipboard.writeText(textToCopy);
        setCopySuccess(true);
        setTimeout(() => setCopySuccess(false), 2000);
        return;
      }

      const textArea = document.createElement("textarea");
      textArea.value = textToCopy;
      textArea.style.position = "fixed";
      textArea.style.left = "-999999px";
      textArea.style.top = "-999999px";
      document.body.appendChild(textArea);
      textArea.focus();
      textArea.select();

      const successful = document.execCommand("copy");
      document.body.removeChild(textArea);

      if (successful) {
        setCopySuccess(true);
        setTimeout(() => setCopySuccess(false), 2000);
        return;
      } else {
        throw new Error("execCommand copy failed");
      }
    } catch (err) {
      console.error("Failed to copy text:", err);
      setCopyModalText(textToCopy);
      setShowCopyModal(true);
    }
  };

  // Load the InChI module for the selected version
  useEffect(() => {
    const checkModuleStatus = async () => {
      try {
        setInchiModuleLoaded(false);
        await loadInchiModule(inchiVersion);
        setInchiModuleLoaded(true);
      } catch (err) {
        console.error(`Failed to load InChI module ${inchiVersion}:`, err);
        setError(
          `Failed to load InChI module ${inchiVersion}. Please check your network connection.`
        );
      }
    };
    checkModuleStatus();
  }, [inchiVersion]);

  // Load an InChI into Ketcher
  const loadInChI = useCallback(async () => {
    if (!inputInchi.trim()) {
      setError("Please enter an InChI string");
      return;
    }

    if (!isEditorReady) {
      setError("Editor not ready. Please try again in a moment.");
      return;
    }

    if (!inputInchi.startsWith("InChI=")) {
      setError('Invalid InChI string. InChI should start with "InChI="');
      return;
    }

    setIsLoading(true);
    setError(null);
    setLogMessage(`Converting InChI to structure...`);

    try {
      const result = await convertInchiToMolfile(inputInchi, "", inchiVersion);

      if (result.return_code === -1) {
        throw new Error(result.message || "Failed to convert InChI to structure");
      }

      if (result.molfile) {
        await executeCommand("setMolecule", [result.molfile]);
        setMolfileContent(result.molfile);
        setInchi(inputInchi);
        setLogMessage(`Structure loaded from InChI successfully.\n${result.message || ""}`);
        generateInChIKeyFromInChI(inputInchi);
      } else {
        throw new Error("No structure data generated from InChI");
      }
    } catch (err) {
      console.error("Failed to load InChI:", err);
      setError(`Failed to convert InChI to structure: ${errorText(err)}`);
    } finally {
      setIsLoading(false);
    }
  }, [inputInchi, isEditorReady, inchiVersion, executeCommand, generateInChIKeyFromInChI]);

  // Load a SMILES string into Ketcher
  const loadSmiles = useCallback(async () => {
    if (!inputSmiles.trim()) {
      setError("Please enter a SMILES string");
      return;
    }

    if (!isEditorReady) {
      setError("Editor not ready. Please try again in a moment.");
      return;
    }

    setIsLoading(true);
    setError(null);
    setLogMessage("Loading SMILES structure...");
    setInchi("");
    setAuxInfo("");
    setInchiKey("");

    try {
      // Ketcher parses the SMILES, not the InChI module
      await executeCommand("setMolecule", [inputSmiles]);
      const molfile = (await getMolfile(inchiVersion)) as string;
      setMolfileContent(molfile);

      setLogMessage(`SMILES loaded successfully. Use "Generate InChI" to create InChI.`);

      if (autoGenerate && molfile) {
        await generateInChIFromMolfile(molfile);
      }
    } catch (err) {
      console.error("Failed to load SMILES:", err);
      setError(`Failed to load structure from SMILES: ${errorText(err)}`);
    } finally {
      setIsLoading(false);
    }
  }, [
    inputSmiles,
    isEditorReady,
    executeCommand,
    getMolfile,
    inchiVersion,
    autoGenerate,
    generateInChIFromMolfile,
  ]);

  // Load a molfile or AuxInfo into Ketcher
  const loadMolfile = useCallback(async () => {
    if (!molfileContent.trim()) {
      setError("Please enter or upload a molfile or AuxInfo");
      return;
    }

    if (!isEditorReady) {
      setError("Editor not ready. Please try again in a moment.");
      return;
    }

    setIsLoading(true);
    setError(null);
    setLogMessage("Processing input...");
    setInchi("");
    setAuxInfo("");
    setInchiKey("");

    try {
      if (molfileContent.startsWith("AuxInfo=")) {
        setLogMessage("Converting AuxInfo to structure...");

        try {
          const result = await convertAuxinfoToMolfile(molfileContent, 0, 0, inchiVersion);

          if (result.molfile) {
            await executeCommand("setMolecule", [result.molfile]);
            setMolfileContent(result.molfile);
            setLogMessage(`AuxInfo converted to structure successfully.\n${result.message || ""}`);

            if (autoGenerate && result.molfile) {
              await generateInChIFromMolfile(result.molfile);
            }
          } else {
            throw new Error("Failed to convert AuxInfo to structure");
          }
        } catch (err) {
          setError(`AuxInfo conversion error: ${errorText(err)}`);
        }
      } else if (molfileContent.includes("V2000") || molfileContent.includes("V3000")) {
        await executeCommand("setMolecule", [molfileContent]);
        setLogMessage("Molfile loaded successfully. Use 'Generate InChI' to create InChI.");

        if (autoGenerate) {
          await generateInChIFromMolfile(molfileContent);
        }
      } else {
        setError("Invalid format. Please enter a valid Molfile or AuxInfo string");
      }
    } catch (err) {
      console.error("Failed to load input:", err);
      setError(`Failed to load structure: ${errorText(err)}`);
    } finally {
      setIsLoading(false);
    }
  }, [
    molfileContent,
    isEditorReady,
    inchiVersion,
    executeCommand,
    autoGenerate,
    generateInChIFromMolfile,
  ]);

  // Read an uploaded molfile into the text field
  const handleMolfileUpload = useCallback((event: ChangeEvent<HTMLInputElement>) => {
    const file = event.target.files?.[0];
    if (!file) return;

    const reader = new FileReader();
    reader.onload = (e) => {
      setMolfileContent(String(e.target?.result ?? ""));
    };
    reader.readAsText(file);
  }, []);

  const handleRetryInit = useCallback(() => {
    resetEditor();
    setError(null);
  }, [resetEditor]);

  const clearEditor = useCallback(async () => {
    if (!isEditorReady) {
      console.warn("Editor not ready for clearing");
      return;
    }

    try {
      await executeCommand("setMolecule", [""]);
      setInchi("");
      setAuxInfo("");
      setInchiKey("");
      setMolfileContent("");
      setLogMessage("Editor cleared");
    } catch (err) {
      console.error("Failed to clear editor:", err);
      setError("Failed to clear the editor");
    }
  }, [isEditorReady, executeCommand]);

  const toggleAutoGenerate = useCallback(() => {
    setAutoGenerate((prev) => !prev);
    if (!autoGenerate) {
      setLogMessage((prev) =>
        prev
          ? `Auto-generate enabled. InChI will update automatically when options change.\n${prev}`
          : "Auto-generate enabled. InChI will update automatically when options change."
      );
    } else {
      setLogMessage((prev) =>
        prev
          ? `Auto-generate disabled. Use "Generate InChI" button to apply changes.\n${prev}`
          : 'Auto-generate disabled. Use "Generate InChI" button to apply changes.'
      );
    }
  }, [autoGenerate]);

  // Common molecules with their InChI and SMILES
  const examples = [
    {
      name: "Ethanol",
      inchi: "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3",
      smiles: "CCO",
      description: "Alcohol",
    },
    {
      name: "Aspirin",
      inchi: "InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)",
      smiles: "CC(=O)OC1=CC=CC=C1C(=O)O",
      description: "Pain reliever",
    },
    {
      name: "Caffeine",
      inchi: "InChI=1S/C8H10N4O2/c1-10-4-9-6-5(10)7(13)12(3)8(14)11(6)2/h4H,1-3H3",
      smiles: "CN1C=NC2=C1C(=O)N(C)C(=O)N2C",
      description: "Stimulant",
    },
    {
      name: "Ibuprofen",
      inchi:
        "InChI=1S/C13H18O2/c1-9(2)8-11-4-6-12(7-5-11)10(3)13(14)15/h4-7,9-10H,8H2,1-3H3,(H,14,15)",
      smiles: "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O",
      description: "Anti-inflammatory",
    },
  ];

  const exampleChips = (field: "smiles" | "inchi") => (
    <div className="pc-smiles-chips">
      <span className="pc-lbl2" id={`inchi-ex-${field}`}>
        Examples
      </span>
      <div className="pc-row" role="group" aria-labelledby={`inchi-ex-${field}`}>
        {examples.map((example) => {
          const value = field === "smiles" ? example.smiles : example.inchi;
          const current = field === "smiles" ? inputSmiles : inputInchi;
          return (
            <button
              key={example.name}
              type="button"
              className="clay clay-outline pc-chip"
              aria-pressed={current === value}
              title={`${example.name}: ${example.description}`}
              onClick={() => (field === "smiles" ? setInputSmiles : setInputInchi)(value)}
            >
              {example.name}
            </button>
          );
        })}
      </div>
    </div>
  );

  // Form submit runs the method's load action, or generation when drawing
  const submit = (e: FormEvent) => {
    e.preventDefault();
    if (activeInputType === "smiles") loadSmiles();
    else if (activeInputType === "molfile") loadMolfile();
    else if (activeInputType === "inchi") loadInChI();
    else generateInChI();
  };

  const ready = isEditorReady && inchiModuleLoaded;
  const hasResult = !!(inchi || inchiKey || auxInfo);
  const resultRows: { label: string; value: string; type: string }[] = [
    { label: "InChI", value: inchi, type: "inchi" },
    { label: "InChIKey", value: inchiKey, type: "key" },
    { label: "AuxInfo", value: auxInfo, type: "auxinfo" },
  ];

  useRevealOnResult(
    "inchi-results",
    isLoading && activeInputType === "structure",
    hasResult,
    !!error
  );

  return (
    <div className="pc-bento">
      {/* Fallback dialog for copying by hand */}
      <Dialog open={showCopyModal} onOpenChange={setShowCopyModal}>
        <DialogContent
          showCloseButton={false}
          className="block w-full max-w-lg gap-0 rounded-xl border-0 bg-card p-6"
        >
          <div className="mb-4 flex items-center justify-between">
            <DialogTitle className="text-lg font-bold text-foreground">Copy Text</DialogTitle>
            <Button
              variant="ghost"
              size="icon-sm"
              onClick={() => setShowCopyModal(false)}
              aria-label="Close"
            >
              <X />
            </Button>
          </div>
          <p className="mb-4 text-foreground">
            Automatic copying failed. Please select and copy this text manually:
          </p>
          <div className="mb-4">
            <Input
              type="text"
              ref={copyTextRef}
              value={copyModalText}
              readOnly
              className="w-full font-mono text-sm"
              onClick={(e) => (e.target as HTMLInputElement).select()}
            />
          </div>
          <div className="flex justify-end gap-3">
            <Button variant="outline" onClick={() => setShowCopyModal(false)}>
              Close
            </Button>
          </div>
        </DialogContent>
      </Dialog>

      <Cell id="inchi-input" title="Structure input" className="pc-t-input">
        <form onSubmit={submit} className="pc-tool-form">
          <div className="pc-opt-group">
            <span className="pc-lbl" id="inchi-method">
              Input method
            </span>
            <Seg
              labelId="inchi-method"
              full
              value={activeInputType}
              onChange={setActiveInputType}
              options={INPUT_TYPES}
            />
          </div>

          {activeInputType === "structure" && (
            <p className="pc-note">
              Draw the structure in the editor below, then generate its InChI. Clear the editor at
              any time to start over.
            </p>
          )}

          {activeInputType === "smiles" && (
            <div className="grid gap-2">
              <label htmlFor="smiles-input" className="pc-lbl">
                SMILES to edit
              </label>
              <input
                id="smiles-input"
                type="text"
                className="pc-well"
                value={inputSmiles}
                onChange={(e) => setInputSmiles(e.target.value)}
                placeholder="e.g., CCO for ethanol"
                spellCheck={false}
                autoComplete="off"
              />
              {exampleChips("smiles")}
            </div>
          )}

          {activeInputType === "molfile" && (
            <div className="grid gap-2">
              <label htmlFor="molfile-input" className="pc-lbl">
                Molfile or AuxInfo
              </label>
              <textarea
                id="molfile-input"
                className="pc-well resize-y"
                value={molfileContent}
                onChange={(e) => setMolfileContent(e.target.value)}
                placeholder="Paste Molfile content here... or AuxInfo=1/0/N:1,2,3/E:(1,2)/rA:3nCCO/rB:s1;s2;/rC:;;;"
                rows={6}
                spellCheck={false}
              />
              <input
                type="file"
                id="molfile-upload"
                accept=".mol,.sdf,.mdl"
                ref={fileInputRef}
                onChange={handleMolfileUpload}
                className="sr-only"
                tabIndex={-1}
                aria-hidden="true"
              />
            </div>
          )}

          {activeInputType === "inchi" && (
            <div className="grid gap-2">
              <label htmlFor="inchi-input-structure" className="pc-lbl">
                InChI to convert
              </label>
              <textarea
                id="inchi-input-structure"
                className="pc-well resize-y"
                value={inputInchi}
                onChange={(e) => setInputInchi(e.target.value)}
                placeholder="e.g., InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3"
                rows={4}
                spellCheck={false}
              />
              {exampleChips("inchi")}
            </div>
          )}

          <div className="pc-row pc-actions">
            <Button
              type={activeInputType === "structure" ? "submit" : "button"}
              onClick={activeInputType === "structure" ? undefined : generateInChI}
              disabled={isLoading || !isEditorReady || !inchiModuleLoaded}
            >
              {isLoading ? <Loader2 className="animate-spin" /> : <FileText />}
              {isLoading ? "Generating…" : "Generate InChI"}
            </Button>
            {activeInputType === "smiles" && (
              <Button
                type="submit"
                variant="outline"
                disabled={isLoading || !isEditorReady || !inputSmiles.trim()}
              >
                {isLoading ? "Loading…" : !isEditorReady ? "Initializing…" : "Load structure"}
              </Button>
            )}
            {activeInputType === "molfile" && (
              <>
                <Button
                  type="submit"
                  variant="outline"
                  disabled={isLoading || !isEditorReady || !molfileContent.trim()}
                >
                  {isLoading ? "Loading…" : "Load structure"}
                </Button>
                <Button
                  type="button"
                  variant="outline"
                  onClick={() => fileInputRef.current?.click()}
                >
                  <Upload />
                  Upload file
                </Button>
              </>
            )}
            {activeInputType === "inchi" && (
              <Button
                type="submit"
                variant="outline"
                disabled={isLoading || !isEditorReady || !inputInchi.trim() || !inchiModuleLoaded}
              >
                {isLoading
                  ? "Loading…"
                  : !isEditorReady
                    ? "Initializing…"
                    : !inchiModuleLoaded
                      ? "Loading InChI module…"
                      : "Convert to structure"}
              </Button>
            )}
          </div>
          <div className="grid gap-2" aria-label="Structure editor" role="region">
            <div className="pc-row">
              <span className={`pc-tag ${isEditorReady ? "text-primary" : "text-saffron-ink"}`}>
                {isEditorReady ? "Editor ready" : "Initializing editor…"}
              </span>
              <span className={`pc-tag ${inchiModuleLoaded ? "text-primary" : "text-saffron-ink"}`}>
                {inchiModuleLoaded ? "InChI module ready" : "Loading InChI module…"}
              </span>
              {!ready && (
                <Button type="button" size="sm" onClick={handleRetryInit}>
                  <RefreshCw />
                  Retry
                </Button>
              )}
              <Button
                type="button"
                variant="outline"
                size="sm"
                onClick={clearEditor}
                disabled={isLoading || !isEditorReady}
              >
                <RefreshCw />
                Clear editor
              </Button>
            </div>
            <div className="h-[75dvh] min-h-[480px] overflow-hidden rounded-xl border border-border bg-card md:h-[900px]">
              <iframe
                ref={ketcherFrame}
                src="/standalone/index.html"
                title="Ketcher Editor"
                className="h-full w-full"
                frameBorder="0"
              />
            </div>
          </div>
        </form>
      </Cell>

      <InChIOptions
        onChange={setOptions}
        inchiVersion={inchiVersion}
        setInchiVersion={setInchiVersion}
      >
        <p className="pc-tool-note">
          The IUPAC International Chemical Identifier (InChI) is a canonical, layered text
          identifier for chemical substances.
        </p>
        <Switch
          checked={autoGenerate}
          onChange={toggleAutoGenerate}
          label="Auto-generate"
          hint="Regenerate the InChI when the options change"
        />
      </InChIOptions>

      <Cell id="inchi-results" title={inchi ? "InChI" : "Results"} className="pc-t-wide">
        <div className="pc-tool-form">
          {error && <GlassErrorCard message={error} onRetry={() => setError(null)} />}
          {!hasResult && !error && (
            <p className="pc-empty">Draw or load a structure, then generate its InChI.</p>
          )}
          {hasResult && (
            <dl className="pc-ids">
              {resultRows.map(({ label, value, type }) => (
                <div className="pc-id" key={label}>
                  <dt>{label}</dt>
                  <dd>
                    {value || <span className="text-muted-foreground">Not generated</span>}
                    <Button
                      variant="outline"
                      size="icon-sm"
                      aria-label={`Copy ${label}`}
                      disabled={!value}
                      onClick={() => copyToClipboard(value, type)}
                    >
                      {copySuccess ? <Check /> : <ClipboardCopy />}
                    </Button>
                  </dd>
                </div>
              ))}
            </dl>
          )}
          <div className="grid gap-1">
            <span className="pc-lbl">Log</span>
            <pre className="max-h-48 overflow-auto rounded-xl border border-border bg-muted p-3 font-mono text-xs whitespace-pre-wrap text-muted-foreground">
              {logMessage || "No log messages yet."}
            </pre>
          </div>
        </div>
      </Cell>
    </div>
  );
};

export default InChIView;
