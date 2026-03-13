import { useState, useEffect } from "react";
import { useDebounce } from "@/hooks/useDebounce";
import { get2DDepictionUrl } from "@/services/depictService";
import { GlassSkeleton } from "@/components/feedback/GlassSkeleton";

interface SMILESPreviewProps {
  smiles: string;
}

export function SMILESPreview({ smiles }: SMILESPreviewProps) {
  const debouncedSmiles = useDebounce(smiles, 500);
  const [loading, setLoading] = useState(true);
  const [hasError, setHasError] = useState(false);

  // Reset loading and error state when debouncedSmiles changes
  useEffect(() => {
    if (debouncedSmiles.trim()) {
      setLoading(true);
      setHasError(false);
    }
  }, [debouncedSmiles]);

  // Render nothing when empty or whitespace-only
  if (!debouncedSmiles.trim()) {
    return null;
  }

  // Render nothing on error (silent hide)
  if (hasError) {
    return null;
  }

  const previewUrl = get2DDepictionUrl(debouncedSmiles.trim(), {
    width: 240,
    height: 240,
    toolkit: "rdkit",
  });

  return (
    <div
      data-testid="smiles-preview"
      className="glass-bold rounded-xl p-2 inline-block animate-in mt-2"
    >
      <div className="relative" style={{ width: 120, height: 120 }}>
        {loading && <GlassSkeleton className="absolute inset-0 h-[120px] w-[120px] rounded-lg" />}
        <img
          src={previewUrl}
          width={120}
          height={120}
          alt="Structure preview"
          className={`rounded-lg transition-opacity duration-200 ${
            loading ? "opacity-0" : "opacity-100"
          }`}
          onLoad={() => setLoading(false)}
          onError={() => {
            setHasError(true);
            setLoading(false);
          }}
        />
      </div>
    </div>
  );
}
