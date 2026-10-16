// 3D depiction: SMILES to a 3D molblock, shown in an interactive $3Dmol.js viewer.
import { useState, useRef, useEffect, useCallback, type FormEvent } from "react";
import SMILESInput from "../common/SMILESInput";
import convertService from "../../services/convertService";
import { Box, Camera, Check, Clipboard, Download, Loader2, LocateFixed } from "lucide-react";
import { ToolSkeleton } from "@/components/feedback/ToolSkeleton";
import { GlassErrorCard } from "@/components/feedback/GlassErrorCard";
import { Button } from "@/components/ui/button";
import { Cell, Seg, Select, Switch } from "@/components/common/porcelain-ui";
import { useRevealOnResult } from "@/hooks/useRevealOnResult";
import { errorText } from "@/lib/error-messages";

type Style3D = "stick" | "line" | "sphere";
type Toolkit3D = "openbabel" | "rdkit";

// Visualization styles
const VISUALIZATION_STYLES: { value: Style3D; label: string }[] = [
  { value: "stick", label: "Stick" },
  { value: "line", label: "Line" },
  { value: "sphere", label: "Sphere" },
];

// Background colours (hex values passed to $3Dmol). The dark presets are the Porcelain
// dark-mode --card (#242125) and --background (#171517) values in tailwind.css.
const BACKGROUND_COLORS = [
  { id: "black", label: "Black", value: "#000000" },
  { id: "white", label: "White", value: "#ffffff" },
  { id: "gray", label: "Dark Gray", value: "#242125" },
  { id: "blue", label: "Charcoal", value: "#171517" },
];

// Colour schemes (keys match $3Dmol options)
const COLOR_SCHEMES = [
  { id: "default", label: "Element (Default)" }, // Maps to undefined in $3Dmol style
  { id: "cyanCarbon", label: "Cyan Carbon" },
  { id: "greenCarbon", label: "Green Carbon" },
  { id: "magentaCarbon", label: "Magenta Carbon" },
  { id: "yellowCarbon", label: "Yellow Carbon" },
  { id: "whiteCarbon", label: "White Carbon" },
  { id: "spectrum", label: "Rainbow Spectrum" },
];

// Toolkits for 3D generation
const TOOLKIT_OPTIONS_3D: { value: Toolkit3D; label: string }[] = [
  { value: "openbabel", label: "Open Babel" },
  { value: "rdkit", label: "RDKit" },
];

// Example molecule SMILES
const EXAMPLE_MOLECULE = "CC1(C)OC2COC3(COS(N)(=O)=O)OC(C)(C)OC3C2O1"; // Example: (+)-Camphorsultam

const Depict3DView = ({ isActive = true }: { isActive?: boolean }) => {
  // isActive controls whether the viewer mounts
  const [smiles, setSmiles] = useState("");
  const [loading, setLoading] = useState(false); // For API calls
  const [error, setError] = useState<string | null>(null);
  const [toolkit, setToolkit] = useState<Toolkit3D>("openbabel");
  const [style, setStyle] = useState<Style3D>("stick");
  const [backgroundColor, setBackgroundColor] = useState("white"); // Background colour id
  const [colorScheme, setColorScheme] = useState("default");
  const [showLabels, setShowLabels] = useState(false);
  const [spin, setSpin] = useState(false);
  const [molData, setMolData] = useState<string | null>(null); // The 3D molblock
  const [viewerInitialized, setViewerInitialized] = useState(false); // $3Dmol viewer is ready
  const [viewerMounted, setViewerMounted] = useState(false); // Container is mounted
  const [copied, setCopied] = useState(false);

  const viewerRef = useRef<Viewer3D | null>(null); // The $3Dmol viewer instance
  const containerRef = useRef<HTMLDivElement | null>(null); // Where the viewer mounts
  const modelRef = useRef<Model3D | null>(null); // The $3Dmol model instance
  const mountTimerRef = useRef<ReturnType<typeof setTimeout> | null>(null); // Delayed mounting

  // Generate a 3D structure from SMILES with the selected toolkit
  const generateStructure = useCallback(
    async (smilesStr: string) => {
      if (!smilesStr) return;
      if (!convertService || typeof convertService.generate3DCoordinates !== "function") {
        console.error("convertService.generate3DCoordinates is not available.");
        setError("3D Conversion service is not configured correctly.");
        return;
      }
      setLoading(true);
      setError(null);
      setMolData(null);
      try {
        const molblockResult = await convertService.generate3DCoordinates(smilesStr, toolkit);
        if (typeof molblockResult === "string" && molblockResult.trim() !== "") {
          setMolData(molblockResult);
        } else {
          throw new Error("Invalid or empty molecule data received from API.");
        }
      } catch (err) {
        console.error("Error generating 3D structure:", err);
        setError(`Error generating 3D structure: ${errorText(err) || "Unknown error"}`);
        setMolData(null);
      } finally {
        setLoading(false);
      }
    },
    [toolkit]
  );

  const handleSubmit = async (e?: FormEvent) => {
    if (e) e.preventDefault();
    const trimmedSmiles = smiles.trim();
    if (trimmedSmiles) {
      await generateStructure(trimmedSmiles);
    } else {
      setError("Please enter a SMILES string.");
    }
  };

  // Destroy and clean up the $3Dmol viewer instance
  const destroyViewer = useCallback(() => {
    if (viewerRef.current) {
      try {
        viewerRef.current.removeAllModels();
        viewerRef.current.removeAllLabels();
        viewerRef.current.clear();
      } catch (e) {
        console.warn("Error cleaning up viewer:", e);
      } finally {
        viewerRef.current = null;
        modelRef.current = null;
        setViewerInitialized(false);
      }
    }
  }, []);

  // Create a new $3Dmol viewer instance
  const createViewer = useCallback(() => {
    const lib = window.$3Dmol;
    if (!containerRef.current || !lib) return null;
    const containerWidth = containerRef.current.clientWidth;
    const containerHeight = containerRef.current.clientHeight;
    if (containerWidth <= 0 || containerHeight <= 0) return null;
    try {
      const bgColorHex =
        BACKGROUND_COLORS.find((c) => c.id === backgroundColor)?.value || "#000000";
      destroyViewer();
      const viewer = lib.createViewer(containerRef.current, {
        backgroundColor: bgColorHex,
        width: "100%",
        height: "100%",
        antialias: true,
        defaultcolors: lib.elementColors.rasmol,
      });
      if (!viewer) throw new Error("$3Dmol.createViewer returned null or undefined.");
      viewer.resize();
      return viewer;
    } catch (err) {
      console.error("Failed to create $3Dmol viewer:", err);
      setError(`Failed to initialize 3D viewer: ${errorText(err)}`);
      return null;
    }
  }, [backgroundColor, destroyViewer]);

  // Apply the visual style settings
  const applyStyle = useCallback(() => {
    if (!viewerRef.current || !modelRef.current) return false;
    const viewer = viewerRef.current;

    try {
      // Reset existing styles and labels before applying new ones
      viewer.setStyle({}, {});
      viewer.removeAllLabels();

      // undefined keeps $3Dmol's default element colouring
      const styleColorScheme = colorScheme === "default" ? undefined : colorScheme;
      const styleConfig = { colorscheme: styleColorScheme };

      switch (style) {
        case "stick":
          viewer.setStyle({}, { stick: { radius: 0.15, ...styleConfig } });
          break;
        case "line":
          viewer.setStyle({}, { line: { linewidth: 2, ...styleConfig } });
          break;
        case "sphere":
          viewer.setStyle({}, { sphere: { scale: 0.3, ...styleConfig } });
          break;
        default:
          viewer.setStyle({}, { stick: { ...styleConfig } });
      }

      // Atom labels for non-hydrogen atoms
      if (showLabels && modelRef.current) {
        try {
          const model = modelRef.current;
          const atoms = model.atoms || [];

          // Label colour chosen for contrast with the viewer background (white is the only light one)
          const labelFontColor = backgroundColor === "white" ? "#000000" : "white";

          // Small offset so labels sit beside the atom centre
          const labelOffset = { x: 0.1, y: 0.1, z: 0 };

          for (let i = 0; i < atoms.length; i++) {
            const atom = atoms[i];
            if (
              atom &&
              atom.elem &&
              atom.elem !== "H" &&
              typeof atom.x === "number" &&
              typeof atom.y === "number" &&
              typeof atom.z === "number"
            ) {
              viewer.addLabel(atom.elem, {
                position: {
                  x: atom.x + labelOffset.x,
                  y: atom.y + labelOffset.y,
                  z: atom.z + labelOffset.z,
                },
                useScreen: false,
                fontColor: labelFontColor,
                fontSize: 14,
                showBackground: false,
                borderThickness: 1.0,
                borderColor: "cyan",
                alignment: "center",
              });
            }
          }
        } catch (labelError) {
          // Labels failing must not block rendering
          console.warn("Could not add atom labels:", labelError);
        }
      }

      viewer.spin(spin);
      viewer.render();
      return true;
    } catch (e) {
      console.error("Error applying style:", e);
      setError(`Failed to apply style: ${errorText(e)}`);
      return false;
    }
  }, [style, colorScheme, showLabels, spin, backgroundColor]);

  // Render the molecule data in the viewer
  const renderMolecule = useCallback(() => {
    if (!viewerRef.current || !molData) return;
    const viewer = viewerRef.current;
    try {
      viewer.removeAllModels();
      modelRef.current = null;
      if (typeof molData !== "string" || !molData.includes("M  END"))
        throw new Error("Invalid Molblock data format.");
      const model = viewer.addModel(molData, "mol");
      if (!model) throw new Error("$3Dmol.addModel failed to return a model.");
      modelRef.current = model;
      const styleApplied = applyStyle();
      if (styleApplied) {
        viewer.zoomTo();
        viewer.render();
      } else {
        viewer.render();
      }
    } catch (e) {
      console.error("Error rendering molecule:", e);
      setError(`Failed to render molecule: ${errorText(e)}. Try reinitializing or refreshing.`);
      destroyViewer();
    }
  }, [molData, applyStyle, destroyViewer]);

  // Initialize the viewer instance
  const initializeViewer = useCallback(() => {
    if (viewerInitialized || !containerRef.current) return false;
    const viewer = createViewer();
    if (viewer) {
      viewerRef.current = viewer;
      setViewerInitialized(true);
      if (molData) setTimeout(() => renderMolecule(), 100);
      return true;
    } else {
      setViewerInitialized(false);
      return false;
    }
  }, [createViewer, molData, renderMolecule, viewerInitialized]);

  // Mount and initialize the viewer
  useEffect(() => {
    if (isActive && !viewerMounted && containerRef.current) {
      mountTimerRef.current = setTimeout(() => {
        setViewerMounted(true);
        initializeViewer();
      }, 150);
    }
    return () => {
      if (mountTimerRef.current) clearTimeout(mountTimerRef.current);
    };
  }, [isActive, viewerMounted, initializeViewer]);

  // Re-render when the molecule data changes
  useEffect(() => {
    if (molData && viewerInitialized) renderMolecule();
    else if (!molData && viewerRef.current) {
      viewerRef.current.removeAllModels();
      modelRef.current = null;
      viewerRef.current.render();
    }
  }, [molData, viewerInitialized, renderMolecule]);

  // Apply style changes
  useEffect(() => {
    if (viewerInitialized && viewerRef.current && modelRef.current) applyStyle();
  }, [style, colorScheme, showLabels, spin, viewerInitialized, applyStyle]);

  // Background colour changes
  useEffect(() => {
    if (viewerRef.current && viewerInitialized) {
      try {
        const bgColorHex =
          BACKGROUND_COLORS.find((c) => c.id === backgroundColor)?.value || "#000000";
        viewerRef.current.setBackgroundColor(bgColorHex);
        viewerRef.current.render();
      } catch (e) {
        console.warn("Error changing background color:", e);
      }
    }
  }, [backgroundColor, viewerInitialized]);

  // Clean up on unmount
  useEffect(() => {
    return () => {
      destroyViewer();
    };
  }, [destroyViewer]);

  // Window resize
  useEffect(() => {
    const handleResize = () => {
      if (viewerRef.current && containerRef.current && viewerInitialized) {
        try {
          viewerRef.current.resize();
          viewerRef.current.render();
        } catch (e) {
          console.warn("Error handling resize:", e);
        }
      }
    };
    window.addEventListener("resize", handleResize);
    return () => window.removeEventListener("resize", handleResize);
  }, [viewerInitialized]);

  const handleUseExampleMolecule = () => {
    setSmiles(EXAMPLE_MOLECULE);
    generateStructure(EXAMPLE_MOLECULE);
  };

  // Save the current view as a PNG
  const handleTakeScreenshot = () => {
    if (!viewerRef.current || !viewerInitialized) {
      setError("Viewer not ready.");
      return;
    }
    try {
      const pngDataUri = viewerRef.current.pngURI();
      if (!pngDataUri) throw new Error("Failed to generate PNG data URI.");
      const a = document.createElement("a");
      a.href = pngDataUri;
      a.download = "molecule-3d-screenshot.png";
      document.body.appendChild(a);
      a.click();
      document.body.removeChild(a);
    } catch (err) {
      console.error("Error taking screenshot:", err);
      setError("Failed to take screenshot.");
    }
  };

  // Reset the camera view
  const handleResetView = () => {
    if (!viewerRef.current || !viewerInitialized || !modelRef.current) return;
    try {
      if (spin) {
        viewerRef.current.spin(false);
      }
      if (viewerRef.current.getModel() && modelRef.current) {
        // Identity matrix
        viewerRef.current.setView(
          new Float32Array([1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1]),
          100
        );
      }
      viewerRef.current.zoomTo();
      if (spin) {
        viewerRef.current.spin(spin);
      }
      viewerRef.current.render();
    } catch (e) {
      console.error("Error resetting view:", e);
      setError(`Failed to reset view: ${errorText(e)}`);
    }
  };

  // Try the viewer again after a failed start. Unmounting is enough: the mount effect
  // schedules a fresh initialization. (Setting `viewerMounted` back to true here as well
  // cancelled that scheduled initialization, so the retry never did anything.)
  const handleReinitialize = () => {
    setError(null);
    setViewerMounted(false);
    destroyViewer();
  };

  const handleCopy = () => {
    if (!molData || !navigator.clipboard) return;
    navigator.clipboard
      .writeText(molData)
      .then(() => {
        setCopied(true);
        setTimeout(() => setCopied(false), 2000);
      })
      .catch((err) => {
        console.error("Failed to copy 3D molblock:", err);
        setError("Failed to copy 3D Molblock.");
      });
  };

  // Download the molblock as a .mol file
  const downloadMolblock = () => {
    if (!molData) return;
    try {
      const blob = new Blob([molData], {
        type: "chemical/x-mdl-molfile;charset=utf-8",
      });
      const url = URL.createObjectURL(blob);
      const a = document.createElement("a");
      a.href = url;
      const filenameBase = smiles
        ? smiles.replace(/[^a-z0-9]/gi, "_").substring(0, 30)
        : "molecule";
      a.download = `${filenameBase}_3d.mol`;
      document.body.appendChild(a);
      a.click();
      document.body.removeChild(a);
      URL.revokeObjectURL(url);
    } catch (err) {
      console.error("Error creating download link:", err);
      setError("Could not create file.");
    }
  };

  const viewerBusy = loading || !molData || !viewerInitialized;
  const done = !!molData && !loading && !error;

  useRevealOnResult("d3d-viewer", loading, !!molData, !!error);

  return (
    <div className="pc-bento">
      <Cell id="d3d-input" title="Structure input" className="pc-t-input">
        <form onSubmit={handleSubmit} className="pc-tool-form">
          <SMILESInput value={smiles} onChange={setSmiles} label="Input SMILES" required />
          <div className="pc-row pc-actions">
            <Button type="submit" disabled={!smiles.trim() || loading}>
              {loading ? <Loader2 className="animate-spin" /> : <Box />}
              {loading ? "Generating…" : "Generate 3D view"}
            </Button>
            {!smiles && !molData && (
              <Button variant="outline" type="button" onClick={handleUseExampleMolecule}>
                Use example molecule
              </Button>
            )}
          </div>
        </form>
      </Cell>

      <Cell id="d3d-opts" title="Options" className="pc-t-side">
        <div className="pc-tool-form">
          <div className="pc-opt-group">
            <span className="pc-lbl" id="d3d-toolkit">
              Toolkit <small>For 3D generation</small>
            </span>
            <Seg
              labelId="d3d-toolkit"
              full
              value={toolkit}
              onChange={setToolkit}
              options={TOOLKIT_OPTIONS_3D}
            />
          </div>
          <div className="pc-opt-group">
            <span className="pc-lbl" id="d3d-style">
              Style
            </span>
            <Seg
              labelId="d3d-style"
              full
              value={style}
              onChange={setStyle}
              options={VISUALIZATION_STYLES}
            />
          </div>
          <div className="pc-field-row">
            <Select
              label="Background"
              value={backgroundColor}
              onChange={setBackgroundColor}
              options={BACKGROUND_COLORS.map((c) => [c.id, c.label])}
            />
            <Select
              label="Colour scheme"
              value={colorScheme}
              onChange={setColorScheme}
              options={COLOR_SCHEMES.map((c) => [c.id, c.label])}
            />
          </div>
          <div className="pc-opt-list">
            <Switch
              checked={showLabels}
              onChange={setShowLabels}
              label="Atom labels"
              hint="Non-hydrogen atoms"
            />
            <Switch
              checked={spin}
              onChange={setSpin}
              label="Spin"
              hint="Rotate the model continuously"
              disabled={viewerBusy}
            />
          </div>
          <p className="pc-tool-note">
            Generates a 3D conformer from SMILES and shows it in an interactive 3Dmol.js viewer.
          </p>
        </div>
      </Cell>

      <Cell
        id="d3d-viewer"
        title="3D viewer"
        className="pc-t-input"
        meta={
          <>
            <Button
              variant="outline"
              size="sm"
              type="button"
              onClick={handleTakeScreenshot}
              disabled={viewerBusy}
            >
              <Camera />
              Screenshot
            </Button>
            <Button
              variant="outline"
              size="sm"
              type="button"
              onClick={handleResetView}
              disabled={viewerBusy}
            >
              <LocateFixed />
              Reset view
            </Button>
          </>
        }
      >
        <div
          className="relative h-[500px] overflow-hidden rounded-xl border border-border"
          style={{
            backgroundColor:
              BACKGROUND_COLORS.find((c) => c.id === backgroundColor)?.value || "#000000",
          }}
        >
          <div
            className="relative h-full w-full"
            ref={containerRef}
            data-testid="molecule-3d-container"
          >
            {/* $3Dmol viewer mounts here */}
          </div>
          {isActive && (!viewerInitialized || !viewerMounted) && (
            <div className="absolute inset-0 z-10 grid place-items-center bg-card/80 backdrop-blur-sm">
              <p className="m-0 grid justify-items-center gap-2 text-center text-foreground">
                <Loader2 className="size-8 animate-spin text-primary" aria-hidden="true" />
                <span className="font-medium">Initializing 3D viewer…</span>
                <span className="text-sm text-muted-foreground">Please wait a moment</span>
              </p>
            </div>
          )}
        </div>
        <p className="pc-note">
          Drag to rotate, scroll or pinch to zoom, right-drag or two-finger drag to pan.
          Double-click or Reset view to recentre.
        </p>
      </Cell>

      <Cell
        id="d3d-results"
        title={done ? "3D molblock" : "Results"}
        className="pc-t-side"
        meta={
          done && (
            <div className="pc-row">
              <Button variant="outline" onClick={handleCopy}>
                {copied ? <Check /> : <Clipboard />}
                {copied ? "Copied" : "Copy molblock"}
              </Button>
              <Button variant="outline" onClick={downloadMolblock}>
                <Download />
                Download (.mol)
              </Button>
            </div>
          )
        }
      >
        {loading && <ToolSkeleton variant="conversion" />}
        {error && !loading && (
          <GlassErrorCard
            message={error}
            onRetry={error.includes("initialize") ? handleReinitialize : () => setError(null)}
          />
        )}
        {!molData && !loading && !error && (
          <p className="pc-empty">Enter a SMILES string to generate and view its 3D structure.</p>
        )}
        {done && (
          <pre className="max-h-[500px] overflow-auto rounded-xl border border-border bg-muted p-3 font-mono text-xs whitespace-pre text-foreground">
            {molData}
          </pre>
        )}
        {molData && !loading && (
          <p className="pc-note">
            Computed with {toolkit === "rdkit" ? "RDKit (ETKDG/MMFF94)" : "OpenBabel (MMFF94)"}; it
            may differ from the conformation found in nature or in experiments. For accurate
            structures use X-ray, NMR or high-level quantum-chemical data.
          </p>
        )}
      </Cell>
    </div>
  );
};

export default Depict3DView;
