// Description: A component for 2D molecule depiction
import { useState, useEffect, useCallback, useId } from "react";

import depictService from "../../services/depictService";
import { errorText } from "@/lib/error-messages";
import { AlertCircle, Download, RefreshCw } from "lucide-react";
import { Button } from "@/components/ui/button";
import { RotationSlider, Select, Switch } from "@/components/common/porcelain-ui";

const H_DISPLAY_OPTIONS: [string, string][] = [
  ["Smart", "Smart"],
  ["Provided", "Provided"],
  ["Minimal", "Minimal"],
  ["Explicit", "Explicit"],
  ["Stereo", "Stereo"],
];

interface MoleculeDepiction2DProps {
  smiles: string;
  title?: string;
  toolkit?: string;
  showCIP?: boolean;
}

const MoleculeDepiction2D = ({
  smiles,
  title,
  toolkit = "rdkit",
  showCIP: initialShowCIP = false,
}: MoleculeDepiction2DProps) => {
  const uid = useId();
  const rotationId = `${uid}-rotation`;
  const [error, setError] = useState<string | null>(null);
  const [imageUrl, setImageUrl] = useState("");
  const [rotation, setRotation] = useState(0);
  // The slider label follows `rotation` live; the image only re-renders once dragging pauses
  const [renderedRotation, setRenderedRotation] = useState(0);
  useEffect(() => {
    const timer = setTimeout(() => setRenderedRotation(rotation), 200);
    return () => clearTimeout(timer);
  }, [rotation]);
  const [useUnicolor, setUseUnicolor] = useState(false);
  const [showCIP, setShowCIP] = useState(initialShowCIP);
  const [hydrogenDisplay, setHydrogenDisplay] = useState("Smart");

  // Fixed width and height - not using state to avoid ESLint warnings
  const imageWidth = 400;
  const imageHeight = 350;

  // Generate the depiction URL - using useCallback for proper memoization
  const generateDepiction = useCallback(() => {
    if (!smiles || !depictService || typeof depictService.get2DDepictionUrl !== "function") {
      setError("Cannot generate 2D depiction: missing SMILES or service unavailable");
      return;
    }

    try {
      const options = {
        toolkit,
        width: imageWidth,
        height: imageHeight,
        rotate: renderedRotation,
        CIP: showCIP,
        unicolor: useUnicolor,
        hydrogen_display: toolkit === "cdk" ? hydrogenDisplay : undefined,
      };

      const url = depictService.get2DDepictionUrl(smiles, options);
      setImageUrl(url);
      setError(null);
    } catch (err) {
      console.error("Error generating 2D depiction:", err);
      setError(`Failed to generate 2D depiction: ${errorText(err)}`);
    }
  }, [smiles, toolkit, renderedRotation, useUnicolor, showCIP, hydrogenDisplay]);

  // Generate depiction when props or settings change
  useEffect(() => {
    if (smiles) {
      generateDepiction();
    }
  }, [smiles, generateDepiction]);

  // Reset every option to the value the component started with
  const handleRefresh = () => {
    setRotation(0);
    setUseUnicolor(false);
    setShowCIP(initialShowCIP);
    setHydrogenDisplay("Smart");
  };

  const handleDownload = () => {
    if (!imageUrl) return;

    const extension = imageUrl.includes("format=png") ? "png" : "svg";
    const filenameSafe = (title || "molecule").replace(/[^a-z0-9]/gi, "_").toLowerCase();

    const a = document.createElement("a");
    a.href = imageUrl;
    a.download = `${filenameSafe}_2d.${extension}`;
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
  };

  return (
    <div className="grid h-full content-start gap-3 rounded-2xl border border-[var(--line)] bg-[var(--inset)] p-3">
      <div className="flex items-center justify-between gap-2">
        <h3 className="m-0 truncate text-sm font-semibold text-[var(--head)]">
          {title || "2D Structure"}
        </h3>
        <Button
          variant="ghost"
          size="icon-sm"
          onClick={handleDownload}
          disabled={!imageUrl}
          title="Download 2D structure"
          aria-label="Download 2D structure"
        >
          <Download aria-hidden="true" />
        </Button>
      </div>

      <div className="pc-stage h-[300px] sm:h-[350px]">
        {imageUrl ? (
          <img
            src={imageUrl}
            alt={title || "2D Molecular Structure"}
            width={imageWidth}
            height={imageHeight}
            onError={() => setError("Failed to load 2D depiction")}
          />
        ) : (
          <p className="pc-stage-empty">{error || "Loading depiction..."}</p>
        )}
      </div>
      {imageUrl && error && (
        <p className="pc-err on m-0" role="alert">
          <AlertCircle className="pc-i size-4 shrink-0" aria-hidden="true" />
          <span>{error}</span>
        </p>
      )}

      <RotationSlider
        id={rotationId}
        label="Rotation"
        value={rotation}
        onChange={setRotation}
        className="mt-0 [&_.pc-i]:size-4"
      />

      <div className="grid gap-x-6 sm:grid-cols-2">
        <Switch label="Black & White" checked={useUnicolor} onChange={setUseUnicolor} />
        <Switch label="Show R/S, E/Z" checked={showCIP} onChange={setShowCIP} />
      </div>

      <div className="pc-row items-end justify-between">
        {toolkit === "cdk" ? (
          <div className="min-w-[140px]">
            <Select
              label="H display"
              value={hydrogenDisplay}
              options={H_DISPLAY_OPTIONS}
              onChange={setHydrogenDisplay}
            />
          </div>
        ) : (
          <span />
        )}
        <Button
          variant="ghost"
          size="icon-sm"
          onClick={handleRefresh}
          title="Reset to original view"
          aria-label="Reset to original view"
        >
          <RefreshCw aria-hidden="true" />
        </Button>
      </div>
    </div>
  );
};
export default MoleculeDepiction2D;
