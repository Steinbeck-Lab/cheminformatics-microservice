// Description: A React component to display SMILES strings with options to copy to clipboard and download as a MOL file.
import { useState } from "react";
import { generate2DCoordinates } from "../../services/convertService";
import { Download } from "lucide-react";
import { Button } from "@/components/ui/button";
import { CopyButton } from "@/components/common/porcelain-ui";

interface SMILESDisplayProps {
  smiles?: string | null;
  label?: string;
  showDownload?: boolean;
}

const SMILESDisplay = ({ smiles, label = "SMILES", showDownload = true }: SMILESDisplayProps) => {
  const [downloading, setDownloading] = useState(false);

  // Handle download as MOL
  const handleDownload = async () => {
    if (!smiles) return;

    setDownloading(true);
    try {
      const molData = await generate2DCoordinates(smiles, "rdkit");

      const blob = new Blob([molData], { type: "chemical/x-mdl-molfile" });
      const url = URL.createObjectURL(blob);
      const a = document.createElement("a");
      a.href = url;
      a.download = "molecule.mol";
      document.body.appendChild(a);
      a.click();
      document.body.removeChild(a);
      URL.revokeObjectURL(url);
    } catch (err) {
      console.error("Failed to convert and download MOL:", err);
    } finally {
      setDownloading(false);
    }
  };

  return (
    <div className="grid w-full gap-1.5">
      <div className="flex items-center justify-between gap-2">
        <span className="text-[13px] font-semibold text-muted-foreground">{label}</span>
        <div className="flex gap-1">
          <CopyButton text={smiles || ""} label="Copy SMILES" icon />
          {showDownload && (
            <Button
              variant="ghost"
              size="icon-sm"
              onClick={handleDownload}
              disabled={downloading || !smiles}
              title="Download as MOL file"
              aria-label="Download as MOL file"
            >
              <Download aria-hidden="true" />
            </Button>
          )}
        </div>
      </div>
      <p className="pc-well pc-mono m-0 min-h-0 py-2.5 text-sm break-all">
        {smiles || "No SMILES data available"}
      </p>
    </div>
  );
};
export default SMILESDisplay;
