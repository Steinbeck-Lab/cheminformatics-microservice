import { useState, useEffect, useRef, type FormEvent } from "react";
import { Check, ClipboardCopy, Loader2, PenLine, RefreshCw, X } from "lucide-react";
import { Button } from "@/components/ui/button";
import { Dialog, DialogContent, DialogTitle } from "@/components/ui/dialog";
import { Input } from "@/components/ui/input";
import { GlassErrorCard } from "@/components/feedback/GlassErrorCard";
import { Cell } from "@/components/common/porcelain-ui";
import { useRevealOnResult } from "@/hooks/useRevealOnResult";
import { useKetcherEditor } from "@/hooks/useKetcherEditor";

// Examples of common molecules
const EXAMPLES = [
  { name: "Ethanol", value: "CCO", description: "Alcohol" },
  { name: "Aspirin", value: "CC(=O)OC1=CC=CC=C1C(=O)O", description: "Pain reliever" },
  { name: "Caffeine", value: "CN1C=NC2=C1C(=O)N(C)C(=O)N2C", description: "Stimulant" },
  { name: "Ibuprofen", value: "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O", description: "Anti-inflammatory" },
];

/** Draw a structure in Ketcher and read it back as SMILES. */
const StructureDrawView = () => {
  const [smiles, setSmiles] = useState("");
  const [inputSmiles, setInputSmiles] = useState("");
  const [copySuccess, setCopySuccess] = useState(false);
  const [showCopyModal, setShowCopyModal] = useState(false);
  const [copyModalText, setCopyModalText] = useState("");
  const [error, setError] = useState<string | null>(null);
  const [isLoading, setIsLoading] = useState(false);
  const copyTextRef = useRef<HTMLInputElement>(null);

  const { ketcherFrame, isEditorReady, executeCommand, resetEditor } = useKetcherEditor();

  // Clear a stale error when the input changes.
  // (`error` used to be a dependency too, which cleared every error the moment it was set.)
  useEffect(() => {
    setError(null);
  }, [inputSmiles]);

  // Auto-select text in the copy dialog when it appears
  useEffect(() => {
    if (showCopyModal && copyTextRef.current) {
      copyTextRef.current.select();
    }
  }, [showCopyModal]);

  // Copy with the Clipboard API, then execCommand, then a dialog to copy by hand
  const copyToClipboard = async (text: string | null = null) => {
    const textToCopy = text || smiles;

    if (!textToCopy) {
      setError("No SMILES to copy. Generate SMILES first.");
      return;
    }

    try {
      // Method 1: the Clipboard API (modern browsers)
      if (navigator.clipboard && navigator.clipboard.writeText) {
        await navigator.clipboard.writeText(textToCopy);
        setCopySuccess(true);
        setTimeout(() => setCopySuccess(false), 2000);
        return;
      }

      // Method 2: execCommand (older browsers)
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
      // Method 3: a dialog with the text to copy by hand
      setCopyModalText(textToCopy);
      setShowCopyModal(true);
    }
  };

  // Load SMILES into Ketcher
  const loadSmiles = async (e?: FormEvent) => {
    e?.preventDefault();
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

    try {
      await executeCommand("setMolecule", [inputSmiles]);
      setSmiles(inputSmiles);
    } catch (err) {
      console.error("Failed to load SMILES:", err);
      setError("Invalid SMILES string or error loading structure. Please check your input.");
    } finally {
      setIsLoading(false);
    }
  };

  // Get SMILES from Ketcher
  const getSmiles = async () => {
    if (!isEditorReady) {
      setError("Editor not ready. Please try again in a moment.");
      return "";
    }

    setIsLoading(true);
    setError(null);

    try {
      const newSmiles = (await executeCommand("getSmiles")) as string;

      if (!newSmiles || newSmiles === "") {
        setError("No structure drawn. Please draw a molecule first.");
        setIsLoading(false);
        return "";
      }

      setSmiles(newSmiles);
      return newSmiles;
    } catch (err) {
      console.error("Failed to generate SMILES:", err);
      setError("Could not generate SMILES from the current structure");
    } finally {
      setIsLoading(false);
    }
    return "";
  };

  // Clear the editor
  const clearEditor = async () => {
    if (!isEditorReady) {
      console.warn("Editor not ready for clearing");
      return;
    }

    try {
      await executeCommand("setMolecule", [""]);
      setSmiles("");
    } catch (err) {
      console.error("Failed to clear editor:", err);
      setError("Failed to clear the editor");
    }
  };

  // Retry initialization
  const handleRetryInit = () => {
    resetEditor();
    setError(null);
  };

  useRevealOnResult("draw-results", !!smiles, !!error);

  return (
    <div className="pc-bento">
      {/* Fallback dialog for copying by hand */}
      <Dialog open={showCopyModal} onOpenChange={setShowCopyModal}>
        <DialogContent
          showCloseButton={false}
          className="block w-full max-w-lg gap-0 rounded-xl border-0 bg-card p-6"
        >
          <div className="mb-4 flex items-center justify-between">
            <DialogTitle className="text-lg font-bold text-foreground">Copy SMILES</DialogTitle>
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

      <div className="pc-col pc-t-stack">
        <Cell id="draw-input" title="Structure input" className="">
          <form onSubmit={loadSmiles} className="pc-tool-form">
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
            </div>
            <div className="pc-smiles-chips">
              <span className="pc-lbl2" id="draw-examples">
                Examples
              </span>
              <div className="pc-row" role="group" aria-labelledby="draw-examples">
                {EXAMPLES.map((example) => (
                  <button
                    key={example.value}
                    type="button"
                    className="clay clay-outline pc-chip"
                    aria-pressed={inputSmiles === example.value}
                    title={`${example.name}: ${example.description}`}
                    onClick={() => setInputSmiles(example.value)}
                  >
                    {example.name}
                  </button>
                ))}
              </div>
            </div>
            <div className="pc-row pc-actions">
              <Button
                type="submit"
                variant="outline"
                disabled={isLoading || !isEditorReady || !inputSmiles.trim()}
              >
                {isLoading ? "Loading…" : !isEditorReady ? "Initializing…" : "Load structure"}
              </Button>
            </div>
          </form>
        </Cell>

        <Cell
          id="draw-results"
          title={smiles ? "SMILES" : "Results"}
          className=""
          meta={
            <>
              <Button onClick={getSmiles} disabled={isLoading || !isEditorReady}>
                {isLoading ? <Loader2 className="animate-spin" /> : <PenLine />}
                Generate SMILES
              </Button>
              <Button variant="outline" onClick={() => copyToClipboard()} disabled={!smiles}>
                {copySuccess ? <Check /> : <ClipboardCopy />}
                {copySuccess ? "Copied" : "Copy SMILES"}
              </Button>
            </>
          }
        >
          <div className="pc-tool-form">
            {error && <GlassErrorCard message={error} onRetry={() => setError(null)} />}
            {smiles ? (
              <p className="pc-well pc-mono m-0 max-h-32 overflow-y-auto break-all">{smiles}</p>
            ) : (
              !error && (
                <p className="pc-empty">
                  Draw a structure in the editor, then generate its SMILES.
                </p>
              )
            )}
          </div>
        </Cell>

        <Cell id="draw-about" title="About the editor" className="">
          <p className="pc-tool-note">
            Draw a structure in the Ketcher editor and generate its SMILES, or load a SMILES string
            to edit it.
          </p>
          <p className="pc-note">
            Generate SMILES reads the drawing; Copy puts the SMILES on your clipboard.
          </p>
        </Cell>
      </div>

      <Cell
        id="draw-editor"
        title="Structure editor"
        className="pc-t-editor"
        meta={
          <>
            <span className={`pc-tag ${isEditorReady ? "text-primary" : "text-saffron-ink"}`}>
              {isEditorReady ? "Editor ready" : "Initializing editor…"}
            </span>
            {!isEditorReady && (
              <Button size="sm" onClick={handleRetryInit}>
                <RefreshCw />
                Retry
              </Button>
            )}
            <Button
              variant="outline"
              size="sm"
              onClick={clearEditor}
              disabled={isLoading || !isEditorReady}
            >
              <RefreshCw />
              Clear editor
            </Button>
          </>
        }
      >
        <div className="h-[75dvh] min-h-[480px] overflow-hidden rounded-xl border border-border bg-card md:h-[700px]">
          <iframe
            ref={ketcherFrame}
            src="/standalone/index.html"
            title="Ketcher Editor"
            className="h-full w-full"
            frameBorder="0"
          />
        </div>
      </Cell>
    </div>
  );
};

export default StructureDrawView;
