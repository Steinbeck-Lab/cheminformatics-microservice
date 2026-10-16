// OCSR: an image (upload, drag and drop, or URL) to SMILES with DECIMER.
import { useState, useCallback, useEffect, useMemo } from "react";
import { useDropzone, type FileRejection } from "react-dropzone";
import MoleculeCard from "../common/MoleculeCard";
import ocsrService from "../../services/ocsrService";
import { CheckCircle, Image, Loader2, RefreshCw, ScanLine, Upload, X } from "lucide-react";
import { ToolSkeleton } from "@/components/feedback/ToolSkeleton";
import { GlassErrorCard } from "@/components/feedback/GlassErrorCard";
import { Button } from "@/components/ui/button";
import { Cell, Seg, Switch } from "@/components/common/porcelain-ui";
import { useRevealOnResult } from "@/hooks/useRevealOnResult";

type InputMethod = "upload" | "url";

const INPUT_METHODS: { value: InputMethod; label: string }[] = [
  { value: "upload", label: "Upload image" },
  { value: "url", label: "Image URL" },
];

const OCRView = () => {
  const [files, setFiles] = useState<File[]>([]);
  const [imageUrl, setImageUrl] = useState("");
  const [results, setResults] = useState<string[]>([]);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [reference, setReference] = useState("");
  const [inputMethod, setInputMethod] = useState<InputMethod>("upload");
  const [handDrawn, setHandDrawn] = useState(false);
  const [processingStage, setProcessingStage] = useState("");

  const onDrop = useCallback((acceptedFiles: File[], rejectedFiles: FileRejection[]) => {
    setFiles([]);
    setImageUrl("");
    setError(null);
    setResults([]);

    const rejected = rejectedFiles[0];
    if (acceptedFiles.length > 0) {
      setFiles(acceptedFiles);
    } else if (rejected) {
      const errorMessages = rejected.errors.map((e) => e.message).join(", ");
      setError(
        `File rejected: ${errorMessages}. Please upload a valid image file (PNG, JPG, or GIF).`
      );
    }
  }, []);

  const { getRootProps, getInputProps, isDragActive, isDragReject } = useDropzone({
    onDrop,
    accept: {
      "image/png": [".png"],
      "image/jpeg": [".jpeg", ".jpg"],
      "image/gif": [".gif"],
    },
    maxFiles: 1,
    multiple: false,
    disabled: loading,
  });

  const changeInputMethod = (method: InputMethod) => {
    setInputMethod(method);
    if (method === "upload") setImageUrl("");
    else setFiles([]);
    setError(null);
  };

  const file = files[0];

  const handleProcessImage = async () => {
    const hasFile = !!file;
    const hasUrl = imageUrl.trim() !== "";

    if (!hasFile && !hasUrl) {
      setError("Please upload an image or provide an image URL.");
      return;
    }

    setLoading(true);
    setError(null);
    setResults([]);
    setProcessingStage("Uploading image...");

    try {
      let response;

      if (file) {
        setProcessingStage("Analyzing structure...");
        if (!ocsrService || typeof ocsrService.processImage !== "function") {
          throw new Error("Image processing service is not available.");
        }
        response = await ocsrService.processImage(file, reference || undefined, handDrawn);
      } else {
        setProcessingStage("Fetching image...");
        if (!ocsrService || typeof ocsrService.processImageUrl !== "function") {
          throw new Error("URL processing service is not available.");
        }
        response = await ocsrService.processImageUrl(
          imageUrl.trim(),
          reference || undefined,
          handDrawn
        );
      }

      setProcessingStage("Extracting structures...");

      if (response && response.smiles) {
        const smilesArray: unknown[] = Array.isArray(response.smiles)
          ? response.smiles
          : [response.smiles];
        const validSmiles = smilesArray.filter(
          (s): s is string => typeof s === "string" && s.trim() !== ""
        );

        if (validSmiles.length > 0) {
          setResults(validSmiles);
          setProcessingStage("");
        } else {
          setError(
            "No chemical structures were detected in the image. Please try a different image or adjust the settings."
          );
        }
      } else {
        setError("Processing completed, but no chemical structures were detected.");
      }
    } catch (err) {
      console.error("Error processing image:", err);
      const e = err as { response?: { data?: { detail?: string } }; message?: string };
      setError(`Error: ${e.response?.data?.detail || e.message || "Unable to process image"}`);
      setResults([]);
    } finally {
      setLoading(false);
      setProcessingStage("");
    }
  };

  const clearAll = () => {
    setFiles([]);
    setImageUrl("");
    setResults([]);
    setError(null);
    setReference("");
    setProcessingStage("");
  };

  // One preview URL per selected file (it used to be recreated on every render)
  const filePreviewUrl = useMemo(() => (file ? URL.createObjectURL(file) : null), [file]);
  const displayImageUrl = filePreviewUrl || (imageUrl.trim() ? imageUrl : null);

  useEffect(() => {
    return () => {
      if (filePreviewUrl) {
        URL.revokeObjectURL(filePreviewUrl);
      }
    };
  }, [filePreviewUrl]);

  const dropState = isDragReject
    ? "border-destructive bg-destructive/8"
    : isDragActive || files.length > 0
      ? "border-primary bg-primary/8"
      : "border-border bg-muted hover:border-primary";
  const done = results.length > 0 && !loading;

  useRevealOnResult("ocr-results", loading, results.length > 0, !!error);

  return (
    <div className="pc-bento">
      <Cell id="ocr-input" title="Image input" className="pc-t-input">
        <div className="pc-tool-form">
          <div className="pc-opt-group">
            <span className="pc-lbl" id="ocr-method">
              Input method
            </span>
            <Seg
              labelId="ocr-method"
              full
              value={inputMethod}
              onChange={changeInputMethod}
              options={INPUT_METHODS}
              disabled={loading}
            />
          </div>

          {inputMethod === "upload" ? (
            <div
              {...getRootProps()}
              className={`grid cursor-pointer justify-items-center gap-2 rounded-2xl border border-dashed p-8 text-center transition-colors ${dropState} ${
                loading ? "cursor-not-allowed opacity-50" : ""
              }`}
            >
              <input {...getInputProps({ "aria-label": "Upload image" })} />
              {file ? (
                <>
                  <CheckCircle className="size-10 text-primary" aria-hidden="true" />
                  <p className="m-0 text-sm font-medium text-primary">File selected: {file.name}</p>
                  <p className="pc-note m-0">{(file.size / 1024).toFixed(2)} KB</p>
                </>
              ) : isDragActive ? (
                <>
                  <Upload className="size-10 text-primary" aria-hidden="true" />
                  <p className="m-0 text-sm font-medium text-primary">Drop your image here</p>
                </>
              ) : (
                <>
                  <Image className="size-10 text-muted-foreground" aria-hidden="true" />
                  <p className="m-0 text-sm font-medium text-foreground">
                    Click to upload or drag and drop
                  </p>
                  <p className="pc-note m-0">PNG, JPG, or GIF (max 10MB)</p>
                </>
              )}
            </div>
          ) : (
            <div className="grid gap-2">
              <label htmlFor="image-url" className="pc-lbl">
                Image URL
              </label>
              <input
                id="image-url"
                type="url"
                className="pc-input"
                value={imageUrl}
                onChange={(e) => {
                  setImageUrl(e.target.value);
                  setFiles([]);
                  setError(null);
                  setResults([]);
                }}
                placeholder="https://example.com/structure.png"
                disabled={loading}
                spellCheck={false}
                autoComplete="off"
              />
            </div>
          )}

          {displayImageUrl && (
            <div className="relative overflow-hidden rounded-xl border border-border bg-muted">
              <img
                src={displayImageUrl}
                alt="Structure preview"
                className="h-auto max-h-96 w-full object-contain"
                onError={() => setError("Failed to load image. Please check the URL.")}
              />
              <Button
                variant="destructive"
                size="icon-sm"
                className="absolute top-2 right-2"
                aria-label="Remove image"
                title="Remove image"
                onClick={(e) => {
                  e.stopPropagation();
                  clearAll();
                }}
              >
                <X />
              </Button>
            </div>
          )}

          <div className="pc-row pc-actions">
            <Button
              type="button"
              onClick={handleProcessImage}
              disabled={loading || (!files.length && !imageUrl.trim())}
            >
              {loading ? <Loader2 className="animate-spin" /> : <ScanLine />}
              {loading ? "Processing…" : "Extract structures"}
            </Button>
            {(files.length > 0 || imageUrl.trim() || results.length > 0 || error) && (
              <Button variant="outline" type="button" onClick={clearAll} disabled={loading}>
                <RefreshCw />
                Reset
              </Button>
            )}
          </div>
        </div>
      </Cell>

      <Cell id="ocr-opts" title="Options" className="pc-t-side">
        <div className="pc-tool-form">
          <div className="grid gap-2">
            <label htmlFor="reference" className="pc-lbl">
              Reference name (optional)
            </label>
            <input
              id="reference"
              type="text"
              className="pc-input"
              value={reference}
              onChange={(e) => setReference(e.target.value)}
              placeholder="e.g., Figure 1A, Compound X"
              disabled={loading}
            />
          </div>
          <Switch
            checked={handDrawn}
            onChange={setHandDrawn}
            label="Hand-drawn model"
            hint="For hand-drawn structures"
            disabled={loading}
          />
          <p className="pc-tool-note">
            DECIMER reads chemical structures in an image and returns them as SMILES. Clear,
            high-contrast images of standard 2D depictions work best.
          </p>
        </div>
      </Cell>

      <Cell
        id="ocr-results"
        title={done ? "Detected structures" : "Results"}
        className="pc-t-wide"
        meta={
          done && (
            <span className="pc-tag">
              {results.length} {results.length === 1 ? "structure" : "structures"}
            </span>
          )
        }
      >
        {loading && (
          <div className="pc-tool-form">
            {processingStage && (
              <p className="pc-note m-0" role="status">
                {processingStage}
              </p>
            )}
            <ToolSkeleton variant="molecule" />
          </div>
        )}
        {error && !loading && <GlassErrorCard message={error} onRetry={() => setError(null)} />}
        {!results.length && !loading && !error && (
          <p className="pc-empty">
            Upload an image or give an image URL to extract its structures.
          </p>
        )}
        {done && (
          <div className="pc-result-grid">
            {results.map((smilesResult, index) => (
              <MoleculeCard
                key={`${index}-${smilesResult}`}
                smiles={smilesResult}
                title={reference || `Structure ${index + 1}`}
                showActions={true}
                size="md"
              />
            ))}
          </div>
        )}
      </Cell>
    </div>
  );
};

export default OCRView;
