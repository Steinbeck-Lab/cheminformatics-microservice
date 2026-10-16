/// <reference types="vite/client" />

// 3Dmol.js loaded via CDN script tag in index.html. Only the parts of its API the
// 3D viewers use; members a viewer checks for before calling are optional.
interface Atom3D {
  elem?: string;
  x?: number;
  y?: number;
  z?: number;
}
interface Model3D {
  atoms?: Atom3D[];
}
interface Point3D {
  x: number;
  y: number;
  z: number;
}
interface Viewer3D {
  removeAllModels: () => void;
  removeAllLabels: () => void;
  clear: () => void;
  resize: () => void;
  render: () => void;
  zoomTo: () => void;
  setStyle: (sel: object, style: object) => void;
  addLabel: (text: string, opts: object) => void;
  spin: (on: boolean) => void;
  setBackgroundColor: (color: string) => void;
  getView?: () => number[];
  setView: (view: number[] | Float32Array, ms?: number) => void;
  resetView?: () => void;
  setRotationCenter: (center: Point3D) => void;
  setCenter: (center: Point3D) => void;
  getModel: () => unknown;
  addModel: (data: string, format: string) => Model3D | null;
  pngURI: () => string;
}
interface Lib3Dmol {
  createViewer: (el: HTMLElement, config: object) => Viewer3D | null;
  elementColors: { rasmol: unknown };
}
/** Undefined until the 3Dmol.js script has loaded. */
declare const $3Dmol: Lib3Dmol | undefined;

// Vite environment variables
interface ImportMetaEnv {
  readonly VITE_API_URL: string;
}

interface ImportMeta {
  readonly env: ImportMetaEnv;
}

// InChI and RInChI WebAssembly builds (Emscripten), loaded by script tag from the
// IUPAC InChI-Web-Demo site. Each script puts a module factory on window.
interface InchiWasmModule {
  ccall(
    ident: string,
    returnType: "number",
    argTypes: string[],
    args: (string | number | boolean)[]
  ): number;
  ccall(
    ident: string,
    returnType: "string",
    argTypes: string[],
    args: (string | number | boolean)[]
  ): string;
  UTF8ToString(ptr: number): string;
  getValue(ptr: number, type: string): number;
  _malloc(size: number): number;
  _free(ptr: number): void;
}

type InchiWasmModuleFactory = () => Promise<InchiWasmModule>;

interface Window {
  $3Dmol?: Lib3Dmol;
  [k: `inchiModule${string}`]: InchiWasmModuleFactory | undefined;
  [k: `rinchiModule${string}`]: InchiWasmModuleFactory | undefined;
}
