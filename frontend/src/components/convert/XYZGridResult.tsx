import { useState } from "react";
import { AlertTriangle, Check, Clipboard, Download, Info } from "lucide-react";
import MoleculeDepiction2D from "../depict/MoleculeDepiction2D";
import { Button } from "@/components/ui/button";
import type { XYZBatchConversionResult, XYZStructureResult } from "@/types/api";

interface XYZGridResultProps {
  result: XYZBatchConversionResult;
}

const downloadBlob = (content: string, filename: string, mime: string) => {
  const blob = new Blob([content], { type: mime });
  const url = URL.createObjectURL(blob);
  const a = document.createElement("a");
  a.href = url;
  a.download = filename;
  document.body.appendChild(a);
  a.click();
  document.body.removeChild(a);
  URL.revokeObjectURL(url);
};

const escapeCsvField = (value: string): string => {
  // RFC 4180: wrap in quotes if contains comma, quote, or newline; double-up internal quotes.
  if (/[",\n\r]/.test(value)) {
    return `"${value.replace(/"/g, '""')}"`;
  }
  return value;
};

const XYZGridResult = ({ result }: XYZGridResultProps) => {
  const { structures, sdf, summary } = result;

  const handleDownloadSdf = () => {
    if (!sdf) return;
    downloadBlob(sdf, "structures.sdf", "chemical/x-mdl-sdfile");
  };

  const handleDownloadCsv = () => {
    const header = "index,title,success,method,canonicalsmiles,inchi,inchikey\n";
    const rows = structures
      .map((s) =>
        [
          String(s.index),
          escapeCsvField(s.title),
          String(s.success),
          escapeCsvField(s.method),
          escapeCsvField(s.canonicalsmiles),
          escapeCsvField(s.inchi),
          escapeCsvField(s.inchikey),
        ].join(",")
      )
      .join("\n");
    downloadBlob(header + rows + "\n", "structures.csv", "text/csv");
  };

  return (
    <div className="space-y-4">
      <div className="flex flex-wrap items-center justify-between gap-3 p-3 bg-primary/8 border border-primary/30 rounded-lg">
        <div data-testid="xyz-summary" className="text-sm text-foreground ">
          <span className="font-semibold">{summary.successful}</span> of {summary.total} structures
          succeeded
          {summary.connectivity_only_count > 0 && (
            <>
              {" "}
              &middot;{" "}
              <span className="text-saffron-ink ">
                {summary.connectivity_only_count} connectivity-only
              </span>
            </>
          )}
          {summary.failed > 0 && (
            <>
              {" "}
              &middot; <span className="text-madder-ink ">{summary.failed} failed</span>
            </>
          )}
        </div>
        <div className="flex gap-2">
          <Button
            variant="outline"
            type="button"
            onClick={handleDownloadSdf}
            disabled={summary.successful === 0}
          >
            <Download className="mr-2 h-4 w-4" />
            Download SDF (all)
          </Button>
          <Button variant="outline" type="button" onClick={handleDownloadCsv}>
            <Download className="mr-2 h-4 w-4" />
            Download CSV
          </Button>
        </div>
      </div>

      <div className="grid grid-cols-1 sm:grid-cols-2 lg:grid-cols-3 gap-4">
        {structures.map((s) => (
          <XYZGridCard key={s.index} structure={s} />
        ))}
      </div>
    </div>
  );
};

interface XYZGridCardProps {
  structure: XYZStructureResult;
}

const XYZGridCard = ({ structure }: XYZGridCardProps) => {
  const [copied, setCopied] = useState(false);

  const handleCopy = async (text: string) => {
    if (!text || !navigator.clipboard) return;
    try {
      await navigator.clipboard.writeText(text);
      setCopied(true);
      setTimeout(() => setCopied(false), 1500);
    } catch {
      // Clipboard write blocked — silently ignore.
    }
  };

  if (!structure.success) {
    return (
      <div className="flex flex-col gap-2 p-3 bg-destructive/8 border border-destructive/30 rounded-lg">
        <div className="flex items-center gap-2 text-sm font-semibold text-madder-ink ">
          <AlertTriangle className="h-4 w-4" />
          {structure.title || `Frame ${structure.index}`}
        </div>
        <div className="text-xs text-madder-ink break-words">{structure.error}</div>
      </div>
    );
  }

  return (
    <div className="flex flex-col gap-2 p-3 bg-card border border-border rounded-lg">
      <div className="flex items-center justify-between gap-2">
        <span className="text-sm font-semibold text-foreground truncate">
          {structure.title || `Frame ${structure.index}`}
        </span>
        {structure.method === "connectivity_only" && (
          <span
            title="Bond orders could not be perceived; only connectivity. Common for transition metals."
            className="inline-flex items-center gap-1 px-2 py-0.5 rounded-full text-[10px] font-medium bg-saffron/8 text-saffron-ink "
          >
            <Info className="h-3 w-3" />
            connectivity-only
          </span>
        )}
      </div>

      <div className="border border-border rounded-md overflow-hidden">
        <MoleculeDepiction2D
          smiles={structure.canonicalsmiles}
          title={structure.title}
          toolkit="cdk"
          showCIP={false}
        />
      </div>

      <div className="flex items-center gap-2">
        <code className="flex-1 px-2 py-1 bg-muted rounded text-xs font-mono break-all text-foreground ">
          {structure.canonicalsmiles}
        </code>
        <Button
          variant="ghost"
          type="button"
          onClick={() => handleCopy(structure.canonicalsmiles)}
          aria-label="Copy SMILES"
        >
          {copied ? <Check className="h-4 w-4" /> : <Clipboard className="h-4 w-4" />}
        </Button>
      </div>

      <details className="text-xs text-muted-foreground ">
        <summary className="cursor-pointer select-none">More</summary>
        <div className="mt-2 space-y-1 break-all">
          <div>
            <span className="font-semibold">InChI:</span> {structure.inchi}
          </div>
          <div>
            <span className="font-semibold">InChIKey:</span> {structure.inchikey}
          </div>
          {structure.warnings.length > 0 && (
            <ul className="list-disc list-inside text-saffron-ink ">
              {structure.warnings.map((w, i) => (
                <li key={i}>{w}</li>
              ))}
            </ul>
          )}
        </div>
      </details>
    </div>
  );
};

export default XYZGridResult;
