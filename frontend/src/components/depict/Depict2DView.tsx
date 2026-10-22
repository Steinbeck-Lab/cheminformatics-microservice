import { useEffect, useMemo, useState } from "react";
import { Download, Play, ClipboardPaste, X, AlertCircle, Layers } from "lucide-react";
import { Button } from "@/components/ui/button";
import {
  Cell,
  CopyButton,
  RotationSlider,
  Select,
  Seg,
  Switch,
} from "@/components/common/porcelain-ui";
import { get2DDepictionUrl, get2DDepictionUrlEnhanced } from "@/services/depictService";
import { calculateNPLikeness } from "@/services/chemService";
import { useMoleculeInsights } from "@/hooks/useMoleculeInsights";
import { useRevealOnResult } from "@/hooks/useRevealOnResult";
import { useDebounced } from "@/hooks/useDebounced";
import { Formula } from "@/components/common/Formula";
import { errorText } from "@/lib/error-messages";

type Toolkit = "cdk" | "rdkit";
type Format = "svg" | "png";
interface Mol {
  smiles: string;
  title: string;
}

const EXAMPLES: readonly [Mol, ...Mol[]] = [
  { title: "Caffeine", smiles: "CN1C=NC2=C1C(=O)N(C(=O)N2C)C" },
  { title: "Aspirin", smiles: "CC(=O)OC1=CC=CC=C1C(=O)O" },
  { title: "Ibuprofen", smiles: "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O" },
  { title: "Paracetamol", smiles: "CC(=O)NC1=CC=C(C=C1)O" },
];
const MAX_MOLS = 50;
/** "Load examples": a batch to draw at once, one "SMILES title" per line (CXSMILES included). */
const BATCH_EXAMPLES = [
  ...EXAMPLES.map((e) => `${e.smiles} ${e.title}`),
  "Cl*.Cl*.c1ccccc1-c1ccccc1 |m:1:4.5.6.7.8.9,3:10.11.12.13.14.15| Dichlorobiphenyl",
  "COc1cc2c(cc1OC)NC(=O)C3=C2C=CN3C Emetine",
  "C1=CC=C(C=C1)C(=O)C(=O)O Phenylglyoxylic acid",
  "C1=CC=C(C=C1)C2=CC=C(C=C2)C(=O)O Biphenyl-4-carboxylic acid",
  "C1=CC=C2C(=C1)C=CC=C2 Naphthalene",
  "CC1=CC=C(C=C1)C Toluene",
].join("\n");

/** A ?smiles= link (from the home page, say) opens that structure; otherwise caffeine. */
function startMol(): Mol {
  const smiles = new URLSearchParams(window.location.search).get("smiles")?.trim();
  return smiles ? { smiles, title: "Your structure" } : EXAMPLES[0];
}
const TOOLKIT_NAME: Record<Toolkit, string> = { cdk: "CDK", rdkit: "RDKit" };
const TOOLKIT_OPTIONS = (["cdk", "rdkit"] as const).map((value) => ({
  value,
  label: TOOLKIT_NAME[value],
}));
const FORMAT_OPTIONS: { value: Format; label: string }[] = [
  { value: "svg", label: "SVG" },
  { value: "png", label: "PNG" },
];

/** One structure per line: "SMILES title". A CXSMILES "|...|" block stays with the SMILES. */
export function parseLines(text: string): Mol[] {
  return text
    .split(/[\r\n]+/)
    .map((l) => l.trim())
    .filter((l) => l && !l.startsWith("#"))
    .slice(0, MAX_MOLS)
    .map((line, i) => {
      const cx = line.match(/^(\S+\s+\|[^|]*\|)\s*(.*)$/);
      const [smiles, rest] = cx
        ? [cx[1] ?? "", cx[2] ?? ""]
        : [line.split(/\s+/)[0] ?? "", line.split(/\s+/).slice(1).join(" ")];
      return { smiles, title: rest || `Molecule ${i + 1}` };
    });
}

interface Options {
  cip: boolean;
  unicolor: boolean;
  hydrogen: string;
  atomNumbers: boolean;
  highlight: string;
  width: number;
  height: number;
  style: string;
  abbreviate: string;
  annotate: string;
  dative: string;
  multicenter: string;
  donuts: boolean;
  flip: boolean;
  showTitle: boolean;
  radicals: boolean;
  zoom: number;
  ratio: number;
}

const DEFAULTS: Options = {
  cip: false,
  unicolor: false,
  hydrogen: "Smart",
  atomNumbers: false,
  highlight: "",
  width: 512,
  height: 512,
  style: "cow",
  abbreviate: "off",
  annotate: "none",
  dative: "metals",
  multicenter: "provided",
  donuts: false,
  flip: false,
  showTitle: false,
  radicals: false,
  zoom: 1,
  ratio: 1,
};

function depictUrl(smiles: string, tk: Toolkit, o: Options, rotate: number, size?: number) {
  const width = size ?? o.width;
  const height = size ?? o.height;
  const common = {
    width,
    height,
    rotate,
    CIP: o.cip,
    unicolor: o.unicolor,
    highlight: o.highlight || undefined,
    showAtomNumbers: o.atomNumbers,
  };
  if (tk === "rdkit") return get2DDepictionUrl(smiles, { ...common, toolkit: "rdkit" });
  return get2DDepictionUrlEnhanced(smiles, {
    ...common,
    hydrogen_display: o.hydrogen,
    abbreviate: o.abbreviate,
    dative: o.dative,
    multicenter: o.multicenter,
    style: o.style,
    annotate: o.annotate,
    donuts: o.donuts,
    zoom: o.zoom,
    ratio: o.ratio,
    flip: o.flip,
    showtitle: o.showTitle,
    perceive_radicals: o.radicals,
  });
}

const safeName = (s: string) => s.replace(/[^a-z0-9]+/gi, "_").toLowerCase() || "molecule";

function saveBlob(blob: Blob, name: string) {
  const href = URL.createObjectURL(blob);
  const a = Object.assign(document.createElement("a"), { href, download: name });
  document.body.appendChild(a);
  a.click();
  a.remove();
  URL.revokeObjectURL(href);
}

/** Fetch the SVG; for PNG draw it onto a canvas (the CDK endpoint only serves SVG). */
async function fetchImage(url: string, format: Format, w: number, h: number): Promise<Blob> {
  const res = await fetch(url);
  if (!res.ok) throw new Error(`The depiction service answered ${res.status}.`);
  const svg = await res.blob();
  if (format === "svg") return new Blob([svg], { type: "image/svg+xml" });
  const src = URL.createObjectURL(new Blob([svg], { type: "image/svg+xml" }));
  try {
    const img = await new Promise<HTMLImageElement>((ok, fail) => {
      const i = new Image();
      i.onload = () => ok(i);
      i.onerror = () => fail(new Error("The SVG could not be turned into a PNG."));
      i.src = src;
    });
    const canvas = Object.assign(document.createElement("canvas"), { width: w * 2, height: h * 2 });
    const ctx = canvas.getContext("2d")!;
    ctx.fillStyle = "#fff";
    ctx.fillRect(0, 0, canvas.width, canvas.height);
    ctx.drawImage(img, 0, 0, canvas.width, canvas.height);
    return await new Promise<Blob>((ok, fail) =>
      canvas.toBlob((b) => (b ? ok(b) : fail(new Error("PNG export failed."))), "image/png")
    );
  } finally {
    URL.revokeObjectURL(src);
  }
}

const minus = (n: number, d = 2) => n.toFixed(d).replace("-", "−");

// NP-likeness of the example chips, fetched once per session for the comparison list.
let exampleScores: Promise<(number | null)[]> | null = null;
const loadExampleScores = () =>
  (exampleScores ??= Promise.all(
    EXAMPLES.map((e) => calculateNPLikeness(e.smiles).then(Number, () => null))
  ));

const TILES: { key: string; label: string; group: string; unit?: string; of?: string }[] = [
  { key: "molecular_weight", label: "Molecular weight", group: "size", unit: "g/mol" },
  { key: "alogp", label: "ALogP", group: "polarity" },
  { key: "topological_polar_surface_area", label: "TPSA", group: "polarity", unit: "Å²" },
  { key: "hbond", label: "H-bond acceptors / donors", group: "polarity" },
  { key: "rotatable_bond_count", label: "Rotatable bonds", group: "size" },
  { key: "aromatic_rings_count", label: "Aromatic rings", group: "size" },
  { key: "qed_drug_likeliness", label: "QED", group: "drug", of: "of 1" },
  {
    key: "lipinski_rule_of_five_violations",
    label: "Lipinski violations",
    group: "drug",
    of: "of 4 rules",
  },
];

const fmt = (v: unknown) =>
  typeof v === "number"
    ? Number.isInteger(v)
      ? String(v)
      : minus(v)
    : v == null
      ? "–"
      : String(v);

/**
 * 2D Depiction: one or many SMILES rendered by CDK or RDKit, compared side by side, with the
 * structure's identifiers, NP-likeness and key descriptors. Opens with caffeine already drawn.
 */
const Depict2DView = () => {
  const [start] = useState(startMol);
  const [input, setInput] = useState(start.smiles);
  const [mols, setMols] = useState<Mol[]>([start]);
  const [toolkit, setToolkit] = useState<Toolkit>("cdk");
  const [opts, setOpts] = useState<Options>(DEFAULTS);
  const [rotation, setRotation] = useState<Record<number, number>>({});
  const [error, setError] = useState<string | null>(null);
  const [busy, setBusy] = useState(false);
  const [broken, setBroken] = useState(false);
  const [exporting, setExporting] = useState(false);
  const [npExamples, setNpExamples] = useState<(number | null)[]>([]);
  const [format, setFormat] = useState<Format>("svg");

  const set = <K extends keyof Options>(k: K, v: Options[K]) => setOpts((o) => ({ ...o, [k]: v }));
  // Typed and dragged values stay live in the controls; the server only draws them once they
  // settle, so a keystroke or a slider step does not each cost a render (or a grid of them).
  const [highlight] = useDebounced(opts.highlight, 250);
  const [width] = useDebounced(opts.width, 250);
  const [height] = useDebounced(opts.height, 250);
  const [zoom] = useDebounced(opts.zoom, 200);
  const [ratio] = useDebounced(opts.ratio, 200);
  const drawn = useMemo<Options>(
    () => ({ ...opts, highlight, width, height, zoom, ratio }),
    [opts, highlight, width, height, zoom, ratio]
  );
  const [drawnRotation, commitRotation] = useDebounced(rotation, 200);
  const resetRotation = () => {
    setRotation({});
    commitRotation({});
  };
  const batch = mols.length > 1; // several structures: a grid of drawings instead of one-molecule panels
  const mol = mols[0] as Mol | undefined;
  const rot = rotation[0] ?? 0;
  const drawnRot = drawnRotation[0] ?? 0;
  const mainSrc = mol && !batch ? depictUrl(mol.smiles, toolkit, drawn, drawnRot) : "";
  const isCdk = toolkit === "cdk";

  // a new src starts loading: clear the old error and show the rendering tag
  useEffect(() => {
    setBroken(false);
    setBusy(!!mainSrc);
  }, [mainSrc]);

  useEffect(() => {
    let alive = true;
    loadExampleScores().then((s) => alive && setNpExamples(s));
    return () => {
      alive = false;
    };
  }, []);

  // the one-molecule panels are hidden in batch mode, so their data is not fetched
  const smiles = batch ? undefined : mol?.smiles;
  const insights = useMoleculeInsights(smiles);

  const depict = (text = input) => {
    const list = parseLines(text);
    if (!list.length) {
      setError("Enter at least one SMILES string.");
      setMols([]);
      return;
    }
    const total = text.split(/[\r\n]+/).filter((l) => l.trim() && !l.trim().startsWith("#")).length;
    setError(total > MAX_MOLS ? `Only the first ${MAX_MOLS} structures are drawn.` : null);
    setMols(list);
    resetRotation();
  };

  const pickExample = (e: Mol) => {
    setInput(e.smiles);
    setMols([e]);
    resetRotation();
    setError(null);
  };

  const loadExamples = () => {
    setInput(BATCH_EXAMPLES);
    depict(BATCH_EXAMPLES);
  };

  const paste = async () => {
    try {
      const text = await navigator.clipboard.readText();
      if (text) setInput(text);
    } catch {
      setError("The browser did not allow reading the clipboard. Paste with Ctrl+V or ⌘V instead.");
    }
  };

  const clear = () => {
    setInput("");
    setMols([]);
    setError(null);
  };

  const downloadMol = async (i: number, fmt: Format) => {
    const m = mols[i];
    if (!m) return;
    try {
      setExporting(true);
      const url = depictUrl(m.smiles, toolkit, drawn, drawnRotation[i] ?? 0);
      const blob = await fetchImage(url, fmt, drawn.width, drawn.height);
      saveBlob(blob, `${safeName(m.title)}_${toolkit}.${fmt}`);
    } catch (e) {
      setError(errorText(e, "The download failed."));
    } finally {
      setExporting(false);
    }
  };

  const downloadAll = async () => {
    try {
      setExporting(true);
      const JSZip = (await import("jszip")).default;
      const zip = new JSZip();
      await Promise.all(
        mols.map(async (m, i) => {
          const url = depictUrl(m.smiles, toolkit, drawn, drawnRotation[i] ?? 0);
          try {
            zip.file(
              `${String(i + 1).padStart(2, "0")}_${safeName(m.title)}.${format}`,
              await fetchImage(url, format, drawn.width, drawn.height)
            );
          } catch (e) {
            zip.file(
              `${String(i + 1).padStart(2, "0")}_${safeName(m.title)}.error.txt`,
              `${m.smiles}\n${e}`
            );
          }
        })
      );
      saveBlob(await zip.generateAsync({ type: "blob" }), `depictions_${toolkit}_${format}.zip`);
    } catch (e) {
      setError(errorText(e, "The download failed."));
    } finally {
      setExporting(false);
    }
  };

  const activeExample =
    mols.length === 1 ? EXAMPLES.find((e) => e.smiles === mols[0]!.smiles) : undefined;
  const d = insights?.descriptors;
  const np = typeof d?.nplikeness === "number" ? d.nplikeness : null;
  const npPos = (v: number) => `${Math.min(100, Math.max(0, ((v + 5) / 10) * 100))}%`;
  const npRows = useMemo(() => {
    const rows = EXAMPLES.map((e, i) => ({
      name: e.title,
      score: npExamples[i] ?? null,
      cur: e.smiles === smiles,
    }));
    if (np !== null && !rows.some((r) => r.cur))
      rows.push({ name: "This structure", score: np, cur: true });
    return rows;
  }, [npExamples, smiles, np]);

  const hydrogenId = "h-display-label";

  useRevealOnResult("h-all", batch);

  return (
    <div className="pc-bento">
      <Cell id="h-input" title="Structure input" meta="SMILES, one per line" className="pc-c-input">
        <label className="pc-lbl" htmlFor="smiles-input">
          SMILES
        </label>
        <textarea
          id="smiles-input"
          className="pc-well"
          rows={mols.length > 1 || input.includes("\n") ? 4 : 2}
          spellCheck={false}
          autoComplete="off"
          autoCapitalize="off"
          value={input}
          onChange={(e) => setInput(e.target.value)}
          onKeyDown={(e) => {
            if (e.key === "Enter" && (e.metaKey || e.ctrlKey)) depict();
          }}
          aria-describedby={error ? "smiles-error" : undefined}
        />
        <div className="pc-row pc-actions">
          <Button onClick={() => depict()}>
            <Play />
            Depict
          </Button>
          <Button variant="outline" onClick={paste}>
            <ClipboardPaste />
            Paste
          </Button>
          <Button variant="outline" onClick={loadExamples}>
            <Layers />
            Load examples
          </Button>
          <Button variant="destructive" onClick={clear}>
            <X />
            Clear
          </Button>
        </div>
        {error && (
          <p className="pc-err on" id="smiles-error" role="alert">
            <AlertCircle className="pc-i" aria-hidden="true" />
            <span>{error}</span>
          </p>
        )}
        <div className="pc-examples">
          <p className="pc-lbl2" id="ex-label">
            Examples
          </p>
          <div className="pc-row" role="group" aria-labelledby="ex-label">
            {EXAMPLES.map((e) => (
              <button
                key={e.title}
                type="button"
                className="clay clay-outline pc-chip"
                aria-pressed={activeExample === e}
                onClick={() => pickExample(e)}
              >
                {e.title}
              </button>
            ))}
          </div>
        </div>
      </Cell>

      <Cell id="h-opts" title="Depiction options" className="pc-c-opts">
        <div className="pc-opt-list">
          <Switch
            checked={opts.cip}
            onChange={(v) => set("cip", v)}
            label="Show CIP labels"
            hint="R/S and E/Z stereo descriptors on atoms and bonds"
          />
          <hr />
          <Switch
            checked={opts.unicolor}
            onChange={(v) => set("unicolor", v)}
            label="Unicolor"
            hint="Draw every atom in one colour instead of by element"
          />
          <hr />
          <div className="pc-opt-group">
            <span className="pc-lbl" id={hydrogenId}>
              Hydrogen display
              <small>
                {isCdk
                  ? "Which hydrogens are drawn explicitly"
                  : "CDK only. Switch the toolkit to CDK."}
              </small>
            </span>
            <Seg
              labelId={hydrogenId}
              full
              value={opts.hydrogen}
              disabled={!isCdk}
              onChange={(v) => set("hydrogen", v)}
              options={["Smart", "Minimal", "Provided", "Explicit", "Stereo"].map((v) => ({
                value: v,
                label: v,
              }))}
            />
          </div>
          <details className="pc-more">
            <summary>More options</summary>
            <div className="pc-more-body">
              <Switch
                checked={opts.atomNumbers}
                onChange={(v) => set("atomNumbers", v)}
                label="Atom numbers"
              />
              <label className="pc-field">
                <span>Highlight substructure (SMARTS)</span>
                <input
                  className="pc-input"
                  value={opts.highlight}
                  placeholder="e.g. c1ccccc1"
                  spellCheck={false}
                  onChange={(e) => set("highlight", e.target.value)}
                />
              </label>
              <div className="pc-field-row">
                {(["width", "height"] as const).map((k) => (
                  <label key={k} className="pc-field">
                    <span>{k === "width" ? "Width (px)" : "Height (px)"}</span>
                    <input
                      className="pc-input"
                      type="number"
                      min={50}
                      max={2000}
                      step={10}
                      value={opts[k]}
                      onChange={(e) =>
                        set(k, Math.min(2000, Math.max(50, Number(e.target.value) || 50)))
                      }
                    />
                  </label>
                ))}
              </div>
              <p className="pc-note">{isCdk ? "CDK drawing options" : "These options need CDK."}</p>
              <div className="pc-field-row">
                <Select
                  label="Style"
                  value={opts.style}
                  disabled={!isCdk}
                  onChange={(v) => set("style", v)}
                  options={[
                    ["cow", "Colour on white"],
                    ["cob", "Colour on black"],
                    ["cot", "Colour on transparent"],
                    ["bow", "Black on white"],
                    ["bot", "Black on transparent"],
                    ["wob", "White on black"],
                    ["nob", "Neon on black"],
                  ]}
                />
                <Select
                  label="Abbreviations"
                  value={opts.abbreviate}
                  disabled={!isCdk}
                  onChange={(v) => set("abbreviate", v)}
                  options={[
                    ["off", "None"],
                    ["groups", "Functional groups"],
                    ["reagents", "Reagents (THF, DMF)"],
                    ["on", "Groups and reagents"],
                  ]}
                />
                <Select
                  label="Annotation"
                  value={opts.annotate}
                  disabled={!isCdk}
                  onChange={(v) => set("annotate", v)}
                  options={[
                    ["none", "None"],
                    ["number", "Atom numbers"],
                    ["cip", "CIP labels"],
                    ["mapidx", "Atom mapping"],
                    ["colmap", "Colour-coded mapping"],
                    ["bondnumber", "Bond numbers"],
                    ["atomvalue", "Atom values"],
                  ]}
                />
                <Select
                  label="Dative bonds"
                  value={opts.dative}
                  disabled={!isCdk}
                  onChange={(v) => set("dative", v)}
                  options={[
                    ["never", "Never"],
                    ["metals", "Metal complexes"],
                    ["always", "Always"],
                  ]}
                />
                <Select
                  label="Multicentre bonds"
                  value={opts.multicenter}
                  disabled={!isCdk}
                  onChange={(v) => set("multicenter", v)}
                  options={[
                    ["provided", "As provided"],
                    ["dative", "Dative arrows"],
                    ["dashed", "Dashed lines"],
                    ["dashed_neutral", "Dashed (neutral)"],
                    ["hidden", "Hidden"],
                    ["hidden_neutral", "Hidden (neutral)"],
                  ]}
                />
              </div>
              <Switch
                checked={opts.donuts}
                disabled={!isCdk}
                onChange={(v) => set("donuts", v)}
                label="Aromatic rings as circles"
              />
              <Switch
                checked={opts.flip}
                disabled={!isCdk}
                onChange={(v) => set("flip", v)}
                label="Flip horizontally"
              />
              <Switch
                checked={opts.showTitle}
                disabled={!isCdk}
                onChange={(v) => set("showTitle", v)}
                label="Show title"
              />
              <Switch
                checked={opts.radicals}
                disabled={!isCdk}
                onChange={(v) => set("radicals", v)}
                label="Perceive radicals"
              />
              {(["zoom", "ratio"] as const).map((k) => (
                <div key={k} className="pc-rot">
                  <label htmlFor={`opt-${k}`}>{k === "zoom" ? "Zoom" : "Bond ratio"}</label>
                  <input
                    id={`opt-${k}`}
                    type="range"
                    min={0.5}
                    max={k === "zoom" ? 3 : 2}
                    step={0.1}
                    value={opts[k]}
                    disabled={!isCdk}
                    onChange={(e) => set(k, Number(e.target.value))}
                  />
                  <output htmlFor={`opt-${k}`}>{opts[k].toFixed(1)}×</output>
                </div>
              ))}
            </div>
          </details>
        </div>
      </Cell>

      {batch ? (
        <Depict2DBatch
          mols={mols}
          toolkit={toolkit}
          onToolkit={setToolkit}
          format={format}
          onFormat={setFormat}
          opts={drawn}
          rotation={rotation}
          drawnRotation={drawnRotation}
          onRotate={(i, deg) => setRotation((o) => ({ ...o, [i]: deg }))}
          onRotateEnd={() => commitRotation()}
          exporting={exporting}
          onDownload={downloadMol}
          onDownloadAll={downloadAll}
        />
      ) : (
        <>
          <Cell
            id="h-depict"
            title="Depiction"
            className="pc-c-depict"
            meta={
              <Seg
                label="Toolkit"
                value={toolkit}
                onChange={setToolkit}
                options={TOOLKIT_OPTIONS}
              />
            }
          >
            <div className={busy ? "pc-stage busy" : "pc-stage"}>
              {mol ? (
                <img
                  key={mainSrc}
                  src={mainSrc}
                  alt={`2D depiction of ${mol.title} rendered by ${TOOLKIT_NAME[toolkit]}`}
                  width={512}
                  height={512}
                  onLoad={() => setBusy(false)}
                  onError={() => {
                    setBusy(false);
                    setBroken(true);
                  }}
                />
              ) : (
                <p className="pc-stage-empty">Enter a SMILES string and press Depict.</p>
              )}
              <span className="pc-loading" aria-hidden="true">
                Rendering…
              </span>
              {mol && (
                <span className="pc-stage-tag" aria-hidden="true">
                  SVG · {opts.width} × {opts.height}
                </span>
              )}
            </div>
            {broken && (
              <p className="pc-err on" role="alert">
                <AlertCircle className="pc-i" aria-hidden="true" />
                <span>
                  The depiction could not be rendered. Check the SMILES or try the other toolkit.
                </span>
              </p>
            )}
            <RotationSlider
              id="rotation"
              label="Rotation"
              value={rot}
              disabled={!mol}
              onChange={(deg) => setRotation((r) => ({ ...r, 0: deg }))}
              onCommit={() => commitRotation()}
            />
            <div className="pc-row pc-actions">
              <Button onClick={() => downloadMol(0, "svg")} disabled={!mol || exporting}>
                <Download />
                Download SVG
              </Button>
              <Button
                variant="outline"
                onClick={() => downloadMol(0, "png")}
                disabled={!mol || exporting}
              >
                <Download />
                PNG
              </Button>
              <CopyButton text={mol?.smiles ?? ""} label="Copy SMILES" />
            </div>
          </Cell>

          <Cell
            id="h-cmp"
            title="Compare toolkits"
            meta="Same SMILES, two renderers"
            className="pc-c-cmp"
          >
            <div className="pc-cmp-grid">
              {(["cdk", "rdkit"] as const).map((tk) => (
                <button
                  key={tk}
                  type="button"
                  className="pc-cmp"
                  aria-pressed={toolkit === tk}
                  aria-label={`${TOOLKIT_NAME[tk]} depiction. Use ${TOOLKIT_NAME[tk]} for the main depiction`}
                  onClick={() => setToolkit(tk)}
                  disabled={!mol}
                >
                  <span className="pc-stage">
                    {mol && (
                      <img
                        alt=""
                        width={256}
                        height={256}
                        src={depictUrl(mol.smiles, tk, drawn, drawnRot, 256)}
                      />
                    )}
                  </span>
                  <span className="pc-cap">
                    {TOOLKIT_NAME[tk]} <small>{toolkit === tk ? "In use" : "Use"}</small>
                  </span>
                </button>
              ))}
            </div>
          </Cell>

          <Cell
            id="h-ids"
            title="Identifiers"
            meta={<span className="pc-tag">RDKit</span>}
            className="pc-c-ids"
          >
            <dl className="pc-ids">
              <div className="pc-id">
                <dt>Molecular formula</dt>
                <dd>
                  {d?.molecular_formula ? (
                    <Formula as="span" value={String(d.molecular_formula)} />
                  ) : (
                    <span className="pc-val pc-wait">{mol ? "Calculating…" : "–"}</span>
                  )}
                  <CopyButton
                    icon
                    text={String(d?.molecular_formula ?? "")}
                    label="Copy molecular formula"
                  />
                </dd>
              </div>
              {(["inchi", "inchikey"] as const).map((k) => (
                <div key={k} className="pc-id">
                  <dt>{k === "inchi" ? "InChI" : "InChIKey"}</dt>
                  <dd>
                    <span className="pc-val">
                      {insights?.[k] ? insights[k] : mol ? "Calculating…" : "–"}
                    </span>
                    <CopyButton
                      icon
                      text={insights?.[k] ?? ""}
                      label={`Copy ${k === "inchi" ? "InChI" : "InChIKey"}`}
                    />
                  </dd>
                </div>
              ))}
            </dl>
            {insights?.failed.ids && (
              <p className="pc-err on">
                <AlertCircle className="pc-i" aria-hidden="true" />
                <span>The identifiers could not be calculated for this structure.</span>
              </p>
            )}
          </Cell>

          <Cell
            id="h-np"
            title="NP-likeness"
            meta={<span className="pc-tag">RDKit</span>}
            className="pc-c-np"
          >
            <p className="pc-np-score">
              <span className="pc-big pc-mono pc-val">{np === null ? "–" : minus(np)}</span>
              <span className="pc-of">on a −5 to +5 scale</span>
            </p>
            <div
              className="pc-scale"
              role="img"
              aria-label={
                np === null
                  ? "NP-likeness score not available"
                  : `NP-likeness score ${minus(np)} on a scale from −5 to +5`
              }
            >
              <div className="pc-bar" />
              <div className="pc-half neg" />
              <div className="pc-half" />
              {Array.from({ length: 11 }, (_, i) => i - 5).map((t) => (
                <span key={t}>
                  <span
                    className={t % 5 === 0 ? "pc-tick maj" : "pc-tick"}
                    style={{ left: npPos(t) }}
                  />
                  {t % 5 === 0 && (
                    <span className="pc-tl" style={{ left: npPos(t) }}>
                      {t > 0 ? `+${t}` : minus(t, 0)}
                    </span>
                  )}
                </span>
              ))}
              {np !== null && (
                <div className="pc-mk" style={{ left: npPos(np) }}>
                  <b>{minus(np)}</b>
                  <i />
                </div>
              )}
            </div>
            <div className="pc-np-ends">
              <span>Synthetic-like</span>
              <span>Natural-product-like</span>
            </div>
            <p className="pc-note">Higher values are more natural-product-like.</p>
            <div className="pc-np-ex">
              <h3 id="np-ex-h">The examples on the same scale</h3>
              <ul aria-labelledby="np-ex-h">
                {npRows.map((r) => (
                  <li key={r.name} className={r.cur ? "cur" : undefined}>
                    <span>{r.name}</span>
                    <span className="pc-tr" aria-hidden="true">
                      {r.score !== null && <i style={{ left: npPos(r.score) }} />}
                    </span>
                    <span className="pc-mono">{r.score === null ? "–" : minus(r.score)}</span>
                  </li>
                ))}
              </ul>
            </div>
            {insights?.failed.descriptors && (
              <p className="pc-err on">
                <AlertCircle className="pc-i" aria-hidden="true" />
                <span>NP-likeness could not be calculated for this structure.</span>
              </p>
            )}
          </Cell>

          <Cell
            id="h-desc"
            title="Descriptors"
            meta="Descriptors from RDKit"
            className="pc-c-legend"
          >
            <ul className="pc-legend" aria-label="Descriptor groups">
              <li data-group="polarity">
                <span className="pc-atom" aria-hidden="true" />
                Polarity <span className="pc-d">lipophilicity, polar surface, H-bonding</span>
              </li>
              <li data-group="size">
                <span className="pc-atom" aria-hidden="true" />
                Size and shape <span className="pc-d">mass, flexibility, rings</span>
              </li>
              <li data-group="drug">
                <span className="pc-atom" aria-hidden="true" />
                Drug-likeness <span className="pc-d">QED, rule of five</span>
              </li>
            </ul>
            {insights?.failed.descriptors && (
              <p className="pc-err on">
                <AlertCircle className="pc-i" aria-hidden="true" />
                <span>The descriptors could not be calculated for this structure.</span>
              </p>
            )}
          </Cell>

          {TILES.map((t) => (
            <section
              key={t.key}
              className="pc-cell glass pc-tile"
              data-group={t.group}
              aria-label={t.label}
            >
              <p className="pc-tl">
                {t.label}
                <span className="pc-atom" aria-hidden="true" />
              </p>
              <p className="pc-num">
                {t.key === "hbond" ? (
                  <>
                    <span className="pc-val">{fmt(d?.hydrogen_bond_acceptors)}</span>
                    <span className="pc-sl">/</span>
                    <span className="pc-val">{fmt(d?.hydrogen_bond_donors)}</span>
                  </>
                ) : (
                  <span className="pc-val">{fmt(d?.[t.key])}</span>
                )}
                {(t.unit || t.of) && <span className="pc-u">{t.unit ?? t.of}</span>}
              </p>
              {t.key === "qed_drug_likeliness" && (
                <div className="pc-meter" aria-hidden="true">
                  <i
                    style={{ width: `${Math.round(Number(d?.qed_drug_likeliness ?? 0) * 100)}%` }}
                  />
                </div>
              )}
            </section>
          ))}
        </>
      )}
    </div>
  );
};

/** Batch mode: every structure drawn in a grid, each with its own rotation and download. */
function Depict2DBatch({
  mols,
  toolkit,
  onToolkit,
  format,
  onFormat,
  opts,
  rotation,
  drawnRotation,
  onRotate,
  onRotateEnd,
  exporting,
  onDownload,
  onDownloadAll,
}: {
  mols: Mol[];
  toolkit: Toolkit;
  onToolkit: (t: Toolkit) => void;
  format: Format;
  onFormat: (f: Format) => void;
  /** The settled options the drawings are made with. */
  opts: Options;
  rotation: Record<number, number>;
  drawnRotation: Record<number, number>;
  onRotate: (i: number, deg: number) => void;
  onRotateEnd: () => void;
  exporting: boolean;
  onDownload: (i: number, fmt: Format) => void;
  onDownloadAll: () => void;
}) {
  return (
    <Cell
      id="h-all"
      title="All structures"
      className="pc-c-all"
      meta={
        <>
          <span>{mols.length} structures</span>
          <Seg label="Toolkit" value={toolkit} onChange={onToolkit} options={TOOLKIT_OPTIONS} />
          <Seg
            label="Download format"
            value={format}
            onChange={onFormat}
            options={FORMAT_OPTIONS}
          />
          <CopyButton text={mols.map((m) => m.smiles).join("\n")} label="Copy all SMILES" />
          <Button onClick={onDownloadAll} disabled={exporting}>
            <Download />
            Download all (ZIP)
          </Button>
        </>
      }
    >
      <ul className="pc-thumbs">
        {mols.map((m, i) => {
          const r = rotation[i] ?? 0;
          return (
            <li key={`${i}-${m.smiles}`} className="pc-thumb">
              <span className="pc-stage">
                <img
                  alt={`2D depiction of ${m.title} rendered by ${TOOLKIT_NAME[toolkit]}`}
                  loading="lazy"
                  width={opts.width}
                  height={opts.height}
                  src={depictUrl(m.smiles, toolkit, opts, drawnRotation[i] ?? 0)}
                />
              </span>
              <span className="pc-thumb-title">{m.title}</span>
              <RotationSlider
                label={`Rotation for ${m.title}`}
                value={r}
                onChange={(deg) => onRotate(i, deg)}
                onCommit={onRotateEnd}
              />
              <div className="pc-row">
                <Button
                  variant="outline"
                  size="sm"
                  onClick={() => onDownload(i, format)}
                  disabled={exporting}
                  aria-label={`Download ${m.title} as ${format.toUpperCase()}`}
                >
                  <Download />
                  {format.toUpperCase()}
                </Button>
                <CopyButton text={m.smiles} label="Copy SMILES" />
              </div>
            </li>
          );
        })}
      </ul>
    </Cell>
  );
}

export default Depict2DView;
