import { useState, useEffect, useRef, type ChangeEvent, type FormEvent } from "react";
import { useKetcherEditor } from "@/hooks/useKetcherEditor";
import {
  DEFAULT_VERSION,
  VERSION_OPTIONS,
  loadRinchiModule,
  convertRxnfileToRinchi,
  generateRinchiKey,
  convertRinchiToFileText,
} from "@/utils/rinchiUtils";
import {
  ArrowLeftRight,
  Check,
  ChevronDown,
  ChevronRight,
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

type InputType = "reaction" | "rxnfile" | "rinchi" | "rinchi-to-file";
type FileFormat = "RXN" | "RD";
type KeyType = "Long" | "Short" | "Web";

const INPUT_TYPES: { value: InputType; label: string }[] = [
  { value: "reaction", label: "Draw" },
  { value: "rxnfile", label: "RXN/RD file" },
  { value: "rinchi", label: "RInChI" },
  { value: "rinchi-to-file", label: "RInChI → file" },
];

/** One copyable identifier; the RInChIKeys fold away. */
const ResultRow = ({
  title,
  value,
  onCopy,
  copyState,
  collapsible = false,
}: {
  title: string;
  value: string;
  onCopy: () => void;
  copyState: boolean;
  collapsible?: boolean;
}) => {
  const [isCollapsed, setIsCollapsed] = useState(false);
  const open = !isCollapsed || !collapsible;

  return (
    <div className="pc-id">
      <dt>
        {collapsible ? (
          <button
            type="button"
            className="inline-flex cursor-pointer items-center gap-1"
            aria-expanded={!isCollapsed}
            onClick={() => setIsCollapsed(!isCollapsed)}
          >
            {isCollapsed ? <ChevronRight className="size-4" /> : <ChevronDown className="size-4" />}
            {title}
          </button>
        ) : (
          title
        )}
      </dt>
      <dd>
        {open && (value || <span className="text-muted-foreground">Not generated</span>)}
        <Button
          variant="outline"
          size="icon-sm"
          aria-label={`Copy ${title}`}
          disabled={!value}
          onClick={onCopy}
        >
          {copyState ? <Check /> : <ClipboardCopy />}
        </Button>
      </dd>
    </div>
  );
};

/** RInChI: draw or load a reaction in Ketcher and convert it with the IUPAC RInChI library. */
const RInChIView = () => {
  const [rinchi, setRinchi] = useState("");
  const [rauxInfo, setRauxInfo] = useState("");
  const [longRinchiKey, setLongRinchiKey] = useState("");
  const [shortRinchiKey, setShortRinchiKey] = useState("");
  const [webRinchiKey, setWebRinchiKey] = useState("");
  const [inputRinchi, setInputRinchi] = useState("");
  const [inputRauxInfo, setInputRauxInfo] = useState("");
  const [rxnfileContent, setRxnfileContent] = useState("");
  const [outputRxnFile, setOutputRxnFile] = useState("");
  const [logMessage, setLogMessage] = useState("");
  const [forceEquilibrium, setForceEquilibrium] = useState(false);
  const [outputFormat, setOutputFormat] = useState<FileFormat>("RXN");
  const [rinchiVersion, setRinchiVersion] = useState(DEFAULT_VERSION);
  const [activeInputType, setActiveInputType] = useState<InputType>("reaction");

  const [isLoading, setIsLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [copySuccess, setCopySuccess] = useState(false);
  const [showCopyModal, setShowCopyModal] = useState(false);
  const [copyModalText, setCopyModalText] = useState("");
  const [rinchiModuleLoaded, setRinchiModuleLoaded] = useState(false);

  const copyTextRef = useRef<HTMLInputElement>(null);
  const fileInputRef = useRef<HTMLInputElement>(null);

  const { ketcherFrame, isEditorReady, executeCommand, resetEditor } = useKetcherEditor();

  // Clear a stale error when the input, version or input method changes.
  // (`error` used to be a dependency too, which cleared every error the moment it was set.)
  useEffect(() => {
    setError(null);
  }, [inputRinchi, rinchiVersion, activeInputType]);

  // Auto-select text in copy modal when it appears
  useEffect(() => {
    if (showCopyModal && copyTextRef.current) {
      copyTextRef.current.select();
    }
  }, [showCopyModal]);

  // Copy with fallbacks: Clipboard API, execCommand, then a dialog to copy by hand
  const copyToClipboard = async (text: string | null = null, type = "rinchi") => {
    let textToCopy;

    switch (type) {
      case "rinchi":
        textToCopy = text || rinchi;
        break;
      case "longkey":
        textToCopy = text || longRinchiKey;
        break;
      case "shortkey":
        textToCopy = text || shortRinchiKey;
        break;
      case "webkey":
        textToCopy = text || webRinchiKey;
        break;
      case "rauxinfo":
        textToCopy = text || rauxInfo;
        break;
      case "rxnfile":
        textToCopy = text || outputRxnFile;
        break;
      default:
        textToCopy = text || rinchi;
    }

    if (!textToCopy) {
      setError(
        `No ${type.replace("key", " RInChIKey").toUpperCase()} to copy. Generate data first.`
      );
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

  // Load the RInChI module for the selected version
  useEffect(() => {
    const checkModuleStatus = async () => {
      try {
        setRinchiModuleLoaded(false);
        await loadRinchiModule(rinchiVersion);
        setRinchiModuleLoaded(true);
      } catch (err) {
        console.error(`Failed to load RInChI module ${rinchiVersion}:`, err);
        setError(
          `Failed to load RInChI module ${rinchiVersion}. Please check your network connection.`
        );
      }
    };
    checkModuleStatus();
  }, [rinchiVersion]);

  // Does Ketcher hold a reaction?
  const checkForReaction = async () => {
    try {
      // The containsReaction method first, if available
      try {
        const ketcher = (
          ketcherFrame.current?.contentWindow as
            (Window & { ketcher?: { containsReaction?: () => boolean } }) | null | undefined
        )?.ketcher;
        if (ketcher?.containsReaction) {
          return ketcher.containsReaction();
        }
      } catch (e) {
        console.debug("Direct Ketcher access failed in checkForReaction:", e);
      }

      // Otherwise: get the RXN file and look for reaction components
      const rxnFile = (await executeCommand("getRxn")) as string;
      return !!rxnFile && (rxnFile.includes("$RXN") || rxnFile.includes("$RDFILE"));
    } catch (err) {
      console.error("Error checking for reaction:", err);
      return false;
    }
  };

  // Does the reaction have an equilibrium arrow? (a simple text check)
  const hasEquilibriumArrow = async () => {
    try {
      const rxnFile = (await executeCommand("getRxn")) as string;
      return !!rxnFile && rxnFile.includes("EQUILIBRIUM");
    } catch (err) {
      console.error("Error checking for equilibrium arrow:", err);
      return false;
    }
  };

  // Generate RInChIKey of one type and store it
  const generateRinchiKeyAndUpdateState = async (rinchiString: string, keyType: KeyType) => {
    if (!rinchiString) {
      console.warn(`Cannot generate ${keyType}-RInChIKey: No RInChI provided`);
      return;
    }

    if (!rinchiString.startsWith("RInChI=")) {
      console.warn(`Invalid RInChI format for ${keyType}-RInChIKey generation:`, rinchiString);
      return;
    }

    try {
      const result = await generateRinchiKey(rinchiString, keyType, rinchiVersion);

      if (!result) {
        console.error(`No result returned from ${keyType}-RInChIKey generation`);
        return;
      }

      if (result.return_code !== 0) {
        console.error(`Failed to generate ${keyType}-RInChIKey:`, result.error);
        return;
      }

      const key = result.rinchikey ?? "";
      switch (keyType) {
        case "Long":
          setLongRinchiKey(key);
          break;
        case "Short":
          setShortRinchiKey(key);
          break;
        case "Web":
          setWebRinchiKey(key);
          break;
        default:
          console.warn(`Unknown RInChIKey type: ${keyType}`);
      }

      return result.rinchikey;
    } catch (err) {
      console.error(`Error generating ${keyType}-RInChIKey:`, err);
      return null;
    }
  };

  // Convert a reaction file to RInChI and store the results
  const convertReactionToRinchiAndWriteResults = async (rxnFile: string, isEquilibrium = false) => {
    try {
      if (!rxnFile || rxnFile.trim() === "") {
        throw new Error("No reaction data provided");
      }

      if (!rinchiModuleLoaded) {
        throw new Error("RInChI module not loaded. Please wait or refresh the page.");
      }

      const result = await convertRxnfileToRinchi(rxnFile, isEquilibrium, rinchiVersion);

      if (!result) {
        throw new Error("No result returned from conversion function");
      }

      if (result.return_code !== 0) {
        throw new Error(result.error || "Failed to generate RInChI");
      }

      setRinchi(result.rinchi || "");
      setRauxInfo(result.rauxinfo || "");

      setLogMessage(
        `RInChI generated successfully${
          isEquilibrium ? " (equilibrium forced)" : ""
        }\n${result.error || ""}`
      );

      if (result.rinchi) {
        await Promise.all([
          generateRinchiKeyAndUpdateState(result.rinchi, "Long"),
          generateRinchiKeyAndUpdateState(result.rinchi, "Short"),
          generateRinchiKeyAndUpdateState(result.rinchi, "Web"),
        ]);
      }

      return result;
    } catch (err) {
      console.error("Failed to generate RInChI from reaction:", err);
      const errorMsg = `Failed to generate RInChI: ${errorText(err)}`;
      setError(errorMsg);
      setLogMessage(errorMsg);
      throw err;
    } finally {
      setIsLoading(false);
    }
  };

  // Generate RInChI from the reaction in Ketcher
  const generateRInChI = async () => {
    if (!isEditorReady) {
      const errorMsg = "Editor not ready. Please try again in a moment.";
      setError(errorMsg);
      setLogMessage(errorMsg);
      return;
    }

    if (!rinchiModuleLoaded) {
      const errorMsg = "RInChI module not loaded. Please wait or refresh the page.";
      setError(errorMsg);
      setLogMessage(errorMsg);
      return;
    }

    setIsLoading(true);
    setError(null);
    setLogMessage("Generating RInChI...");
    setRinchi("");
    setRauxInfo("");
    setLongRinchiKey("");
    setShortRinchiKey("");
    setWebRinchiKey("");

    try {
      const hasReaction = await checkForReaction();
      if (!hasReaction) {
        const errorMsg = "No reaction found. Please draw a reaction with reactants and products.";
        setError(errorMsg);
        setLogMessage(errorMsg);
        setIsLoading(false);
        return;
      }

      const rxnFile = (await executeCommand("getRxn")) as string;

      if (!rxnFile || rxnFile.trim() === "") {
        const errorMsg =
          "Failed to get reaction data. Please ensure you've drawn a valid reaction.";
        setError(errorMsg);
        setLogMessage(errorMsg);
        setIsLoading(false);
        return;
      }

      setRxnfileContent(rxnFile);

      // Equilibrium: forced, or read from the arrow
      let isEquilibrium = forceEquilibrium;
      if (!isEquilibrium) {
        isEquilibrium = await hasEquilibriumArrow();
      }

      await convertReactionToRinchiAndWriteResults(rxnFile, isEquilibrium);
    } catch (err) {
      console.error("Failed to generate RInChI:", err);
      setLogMessage(`Failed to generate RInChI: ${errorText(err)}`);
      setError(`Failed to generate RInChI: ${errorText(err)}`);
      setIsLoading(false);
    }
  };

  // Read an uploaded RXN/RD file into the text field
  const handleRxnFileUpload = (event: ChangeEvent<HTMLInputElement>) => {
    const file = event.target.files?.[0];
    if (!file) return;

    const reader = new FileReader();
    reader.onload = (e) => {
      setRxnfileContent(String(e.target?.result ?? ""));
    };
    reader.readAsText(file);
  };

  // Convert a pasted or uploaded RXN/RD file
  const processRxnFile = async () => {
    if (!rxnfileContent.trim()) {
      const errorMsg = "Please enter or upload a RXN/RD file";
      setError(errorMsg);
      setLogMessage(errorMsg);
      return;
    }

    if (!rinchiModuleLoaded) {
      const errorMsg = "RInChI module not loaded. Please wait or refresh the page.";
      setError(errorMsg);
      setLogMessage(errorMsg);
      return;
    }

    setIsLoading(true);
    setError(null);
    setLogMessage("Processing RXN/RD file...");
    setRinchi("");
    setRauxInfo("");
    setLongRinchiKey("");
    setShortRinchiKey("");
    setWebRinchiKey("");

    try {
      if (!rxnfileContent.includes("$RXN") && !rxnfileContent.includes("$RDFILE")) {
        const errorMsg = "Invalid RXN/RD file format. Please check your input.";
        setError(errorMsg);
        setLogMessage(errorMsg);
        setIsLoading(false);
        return;
      }

      await convertReactionToRinchiAndWriteResults(rxnfileContent, forceEquilibrium);

      // Show the reaction in Ketcher too, when the editor is ready
      if (isEditorReady) {
        try {
          await executeCommand("setMolecule", [rxnfileContent]);
          setLogMessage((prevLog) => `Reaction loaded in editor.\n${prevLog}`);
        } catch (ketcherErr) {
          // Only the drawing failed: keep the result
          console.error("Failed to load reaction in Ketcher:", ketcherErr);
        }
      }
    } catch (err) {
      console.error("Failed to process RXN/RD file:", err);
      setLogMessage(`Failed to process RXN/RD file: ${errorText(err)}`);
      setError(`Failed to process RXN/RD file: ${errorText(err)}`);
    } finally {
      setIsLoading(false);
    }
  };

  // Load a RInChI into Ketcher
  const loadRInChI = async () => {
    if (!inputRinchi.trim()) {
      setError("Please enter a RInChI string");
      return;
    }

    if (!isEditorReady) {
      setError("Editor not ready. Please try again in a moment.");
      return;
    }

    if (!rinchiModuleLoaded) {
      setError("RInChI module not loaded. Please wait or refresh the page.");
      return;
    }

    if (!inputRinchi.startsWith("RInChI=")) {
      setError('Invalid RInChI string. RInChI should start with "RInChI="');
      return;
    }

    setIsLoading(true);
    setError(null);
    setLogMessage("Converting RInChI to reaction...");

    try {
      const result = await convertRinchiToFileText(
        inputRinchi,
        inputRauxInfo || "",
        "RXN",
        rinchiVersion
      );

      if (result.return_code !== 0) {
        throw new Error(result.error || "Failed to convert RInChI to reaction");
      }

      if (!result.fileText) {
        throw new Error("No reaction data generated from RInChI");
      }

      await executeCommand("setMolecule", [result.fileText]);

      setRinchi(inputRinchi);
      setRauxInfo(inputRauxInfo || "");
      setRxnfileContent(result.fileText);
      setLogMessage(
        `Reaction loaded from RInChI successfully${result.error ? `\n${result.error}` : ""}`
      );

      if (inputRinchi) {
        await Promise.all([
          generateRinchiKeyAndUpdateState(inputRinchi, "Long"),
          generateRinchiKeyAndUpdateState(inputRinchi, "Short"),
          generateRinchiKeyAndUpdateState(inputRinchi, "Web"),
        ]);
      }
    } catch (err) {
      console.error("Failed to load RInChI:", err);
      setLogMessage(`Failed to convert RInChI to reaction: ${errorText(err)}`);
      setError(`Failed to convert RInChI to reaction: ${errorText(err)}`);
    } finally {
      setIsLoading(false);
    }
  };

  // Convert a RInChI to a RXN or RD file
  const convertRInChIToFile = async () => {
    if (!inputRinchi.trim()) {
      setError("Please enter a RInChI string");
      return;
    }

    if (!rinchiModuleLoaded) {
      setError("RInChI module not loaded. Please wait or refresh the page.");
      return;
    }

    if (!inputRinchi.startsWith("RInChI=")) {
      setError('Invalid RInChI string. RInChI should start with "RInChI="');
      return;
    }

    setIsLoading(true);
    setError(null);
    setLogMessage(`Converting RInChI to ${outputFormat} file...`);
    setOutputRxnFile("");

    try {
      const result = await convertRinchiToFileText(
        inputRinchi,
        inputRauxInfo || "",
        outputFormat,
        rinchiVersion
      );

      if (result.return_code !== 0) {
        throw new Error(result.error || `Failed to convert RInChI to ${outputFormat} file`);
      }

      if (!result.fileText) {
        throw new Error(`No ${outputFormat} file data generated from RInChI`);
      }

      setOutputRxnFile(result.fileText);
      setLogMessage(
        `${outputFormat} file generated successfully${result.error ? `\n${result.error}` : ""}`
      );
    } catch (err) {
      console.error("Failed to convert RInChI to file:", err);
      setLogMessage(`Failed to convert RInChI to ${outputFormat} file: ${errorText(err)}`);
      setError(`Failed to convert RInChI to ${outputFormat} file: ${errorText(err)}`);
    } finally {
      setIsLoading(false);
    }
  };

  const clearEditor = async () => {
    if (!isEditorReady) {
      console.warn("Editor not ready for clearing");
      return;
    }

    try {
      await executeCommand("setMolecule", [""]);
      setRinchi("");
      setRauxInfo("");
      setLongRinchiKey("");
      setShortRinchiKey("");
      setWebRinchiKey("");
      setRxnfileContent("");
      setLogMessage("Editor cleared");
    } catch (err) {
      console.error("Failed to clear editor:", err);
      setError("Failed to clear the editor");
    }
  };

  const handleRetryInit = () => {
    resetEditor();
    setError(null);
  };

  // Common reactions with their RInChI
  const examples = [
    {
      name: "Esterification",
      rinchi:
        "RInChI=1.00.1S/C2H4O2/c1-2(3)4/h1H3,(H,3,4)!C2H6O/c1-2-3/h3H,2H2,1H3<>C4H8O2/c1-3-6-4(2)5/h3H2,1-2H3!H2O/h1H2/d+",
      description: "Acetic acid + ethanol → ethyl acetate + water",
    },
    {
      name: "Diels-Alder",
      rinchi:
        "RInChI=1.00.1S/C4H6/c1-3-4-2/h3-4H,1-2H2!C2H2/c1-2/h1-2H<>C6H8/c1-2-4-6-5-3-1/h1-6H/d-",
      description: "1,3-butadiene + ethyne → cyclohexadiene",
    },
    {
      name: "Equilibrium",
      rinchi:
        "RInChI=1.00.1S/C6H12O6/c7-1-2-3(8)4(9)5(10)6(11)12-2/h2-11H,1H2/t2-,3-,4+,5-,6?/m1/s1<>C6H12O6/c7-1-3-4(9)5(10)6(11)2-8/h3-11H,1-2H2/t3-,4+,5-,6-/m1/s1/d=",
      description: "Glucose ⇌ fructose",
    },
  ];

  // Form submit runs the method's action
  const submit = (e: FormEvent) => {
    e.preventDefault();
    if (activeInputType === "rxnfile") processRxnFile();
    else if (activeInputType === "rinchi") loadRInChI();
    else if (activeInputType === "rinchi-to-file") convertRInChIToFile();
    else generateRInChI();
  };

  const toFile = activeInputType === "rinchi-to-file";
  const ready = isEditorReady && rinchiModuleLoaded;
  const hasResult = !!(rinchi || rauxInfo || longRinchiKey || shortRinchiKey || webRinchiKey);
  const hasFile = toFile && !!outputRxnFile;
  const resultTitle = hasFile ? `${outputFormat} file` : rinchi ? "RInChI" : "Results";
  const fieldSuffix = toFile ? "file" : "structure";

  const rinchiFields = (
    <>
      <div className="grid gap-2">
        <label htmlFor={`rinchi-input-${fieldSuffix}`} className="pc-lbl">
          RInChI to convert
        </label>
        <textarea
          id={`rinchi-input-${fieldSuffix}`}
          className="pc-well resize-y"
          value={inputRinchi}
          onChange={(e) => setInputRinchi(e.target.value)}
          placeholder="e.g., RInChI=1.00.1S/C2H4O2/c1-2(3)4/h1H3,(H,3,4)!C2H6O/c1-2-3/h3H,2H2,1H3<>C4H8O2/c1-3-6-4(2)5/h3H2,1-2H3!H2O/h1H2/d+"
          rows={4}
          spellCheck={false}
        />
      </div>
      <div className="grid gap-2">
        <label htmlFor={`rauxinfo-input-${fieldSuffix}`} className="pc-lbl">
          RAuxInfo (optional)
        </label>
        <textarea
          id={`rauxinfo-input-${fieldSuffix}`}
          className="pc-well resize-y"
          value={inputRauxInfo}
          onChange={(e) => setInputRauxInfo(e.target.value)}
          placeholder="e.g., RAuxInfo=1.00.1/1/N:1,2,3,4/E:(3,4)/rA:4nCCOO/rB:s1;s2;d2;/rC:-3.8549,-.5552,0;..."
          rows={3}
          spellCheck={false}
        />
      </div>
      <div className="pc-smiles-chips">
        <span className="pc-lbl2" id="rinchi-ex">
          Examples
        </span>
        <div className="pc-row" role="group" aria-labelledby="rinchi-ex">
          {examples.map((example) => (
            <button
              key={example.name}
              type="button"
              className="clay clay-outline pc-chip"
              aria-pressed={inputRinchi === example.rinchi}
              title={`${example.name}: ${example.description}`}
              onClick={() => setInputRinchi(example.rinchi)}
            >
              {example.name}
            </button>
          ))}
        </div>
      </div>
    </>
  );

  useRevealOnResult("rinchi-results", isLoading, hasResult, hasFile, !!error);

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

      <Cell id="rinchi-input" title="Reaction input" className="pc-t-input">
        <form onSubmit={submit} className="pc-tool-form">
          <div className="pc-opt-group">
            <span className="pc-lbl" id="rinchi-method">
              Input method
            </span>
            <Seg
              labelId="rinchi-method"
              full
              value={activeInputType}
              onChange={setActiveInputType}
              options={INPUT_TYPES}
            />
          </div>

          {activeInputType === "reaction" && (
            <p className="pc-note">
              Draw reactants on the left, add a reaction arrow from the toolbar and draw products on
              the right, then generate the RInChI. Turn on Force equilibrium for equilibrium
              reactions.
            </p>
          )}

          {activeInputType === "rxnfile" && (
            <div className="grid gap-2">
              <label htmlFor="rxnfile-input" className="pc-lbl">
                RXN/RD file content
              </label>
              <textarea
                id="rxnfile-input"
                className="pc-well resize-y"
                value={rxnfileContent}
                onChange={(e) => setRxnfileContent(e.target.value)}
                placeholder="Paste RXN/RD file content here..."
                rows={6}
                spellCheck={false}
              />
              <input
                type="file"
                id="rxnfile-upload"
                accept=".rxn,.rd"
                ref={fileInputRef}
                onChange={handleRxnFileUpload}
                className="sr-only"
                tabIndex={-1}
                aria-hidden="true"
              />
            </div>
          )}

          {(activeInputType === "rinchi" || toFile) && rinchiFields}

          <div className="pc-row pc-actions">
            {activeInputType === "reaction" && (
              <Button type="submit" disabled={isLoading || !isEditorReady || !rinchiModuleLoaded}>
                {isLoading ? <Loader2 className="animate-spin" /> : <FileText />}
                {isLoading ? "Generating…" : "Generate RInChI"}
              </Button>
            )}
            {activeInputType === "rxnfile" && (
              <>
                <Button
                  type="submit"
                  disabled={isLoading || !rxnfileContent.trim() || !rinchiModuleLoaded}
                >
                  {isLoading ? <Loader2 className="animate-spin" /> : <FileText />}
                  {isLoading ? "Processing…" : "Process RXN/RD file"}
                </Button>
                <Button
                  type="button"
                  variant="outline"
                  onClick={() => fileInputRef.current?.click()}
                >
                  <Upload />
                  Upload RXN/RD file
                </Button>
              </>
            )}
            {activeInputType === "rinchi" && (
              <Button
                type="submit"
                disabled={isLoading || !isEditorReady || !inputRinchi.trim() || !rinchiModuleLoaded}
              >
                {isLoading
                  ? "Converting…"
                  : !isEditorReady
                    ? "Initializing…"
                    : !rinchiModuleLoaded
                      ? "Loading RInChI module…"
                      : "Convert to reaction"}
              </Button>
            )}
            {toFile && (
              <Button
                type="submit"
                disabled={isLoading || !inputRinchi.trim() || !rinchiModuleLoaded}
              >
                {isLoading ? <Loader2 className="animate-spin" /> : <ArrowLeftRight />}
                {isLoading
                  ? "Converting…"
                  : !rinchiModuleLoaded
                    ? "Loading RInChI module…"
                    : `Convert to ${outputFormat} file`}
              </Button>
            )}
            {activeInputType !== "reaction" && !toFile && (
              <Button
                type="button"
                variant="outline"
                onClick={generateRInChI}
                disabled={isLoading || !isEditorReady || !rinchiModuleLoaded}
              >
                <FileText />
                Generate RInChI
              </Button>
            )}
          </div>
          {/* The editor is unmounted for RInChI → file, as before. */}
          {!toFile && (
            <div className="grid gap-2" role="region" aria-label="Reaction editor">
              <div className="pc-row">
                <span className={`pc-tag ${isEditorReady ? "text-primary" : "text-saffron-ink"}`}>
                  {isEditorReady ? "Editor ready" : "Initializing editor…"}
                </span>
                <span
                  className={`pc-tag ${rinchiModuleLoaded ? "text-primary" : "text-saffron-ink"}`}
                >
                  {rinchiModuleLoaded ? "RInChI module ready" : "Loading RInChI module…"}
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
              <div className="h-[75dvh] min-h-[480px] overflow-hidden rounded-xl border border-border bg-card md:h-[800px]">
                <iframe
                  ref={ketcherFrame}
                  src="/standalone/index.html"
                  title="Ketcher Editor"
                  className="h-full w-full"
                  frameBorder="0"
                />
              </div>
            </div>
          )}
        </form>
      </Cell>

      <div className="pc-col pc-t-side">
        <Cell id="rinchi-opts" title="Options" className="">
          <div className="pc-tool-form">
            <p className="pc-tool-note">
              The IUPAC RInChI (Reaction InChI) is a text identifier for chemical reactions, built
              from the InChIs of reactants and products and the direction of the reaction.
            </p>
            <Select
              label="RInChI version"
              value={rinchiVersion}
              options={VERSION_OPTIONS}
              onChange={setRinchiVersion}
            />
            <Switch
              checked={forceEquilibrium}
              onChange={setForceEquilibrium}
              label="Force equilibrium"
              hint="Treat the reaction as an equilibrium, whatever the arrow type"
            />
            <div className="pc-opt-group">
              <span className="pc-lbl" id="rinchi-format">
                File format
              </span>
              <Seg
                labelId="rinchi-format"
                full
                value={outputFormat}
                onChange={setOutputFormat}
                options={[
                  { value: "RXN", label: "RXN" },
                  { value: "RD", label: "RD" },
                ]}
              />
              <p className="pc-note mt-2">Output format for RInChI to file conversion.</p>
            </div>
          </div>
        </Cell>

        <Cell id="rinchi-results" title={resultTitle} className="">
          <div className="pc-tool-form">
            {error && <GlassErrorCard message={error} onRetry={() => setError(null)} />}
            {!hasResult && !hasFile && !error && (
              <p className="pc-empty">
                {toFile
                  ? "Enter a RInChI to convert it to a reaction file."
                  : "Draw or load a reaction, then generate its RInChI."}
              </p>
            )}

            {hasFile && (
              <div className="grid gap-2">
                <span className="pc-lbl">{outputFormat} file output</span>
                {/* Fixed-width columns: keep every space */}
                <pre
                  data-testid="rinchi-file-output"
                  className="max-h-96 overflow-auto rounded-xl border border-border bg-muted p-3 font-mono text-xs whitespace-pre text-foreground"
                >
                  {outputRxnFile}
                </pre>
                <div className="pc-row">
                  <Button
                    variant="outline"
                    onClick={() => copyToClipboard(outputRxnFile, "rxnfile")}
                  >
                    {copySuccess ? <Check /> : <ClipboardCopy />}
                    {`Copy ${outputFormat} file`}
                  </Button>
                </div>
              </div>
            )}

            {hasResult && (
              <>
                <dl className="pc-ids">
                  <ResultRow
                    title="RInChI"
                    value={rinchi}
                    onCopy={() => copyToClipboard(rinchi, "rinchi")}
                    copyState={copySuccess}
                  />
                  <ResultRow
                    title="RAuxInfo"
                    value={rauxInfo}
                    onCopy={() => copyToClipboard(rauxInfo, "rauxinfo")}
                    copyState={copySuccess}
                  />
                </dl>
                <div className="grid gap-2">
                  <span className="pc-lbl">RInChIKeys</span>
                  <dl className="pc-ids">
                    <ResultRow
                      title="Long-RInChIKey"
                      value={longRinchiKey}
                      onCopy={() => copyToClipboard(longRinchiKey, "longkey")}
                      copyState={copySuccess}
                      collapsible
                    />
                    <ResultRow
                      title="Short-RInChIKey"
                      value={shortRinchiKey}
                      onCopy={() => copyToClipboard(shortRinchiKey, "shortkey")}
                      copyState={copySuccess}
                      collapsible
                    />
                    <ResultRow
                      title="Web-RInChIKey"
                      value={webRinchiKey}
                      onCopy={() => copyToClipboard(webRinchiKey, "webkey")}
                      copyState={copySuccess}
                      collapsible
                    />
                  </dl>
                </div>
              </>
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
    </div>
  );
};

export default RInChIView;
