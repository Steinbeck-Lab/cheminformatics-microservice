/// <reference types="vite/client" />

// 3Dmol.js loaded via CDN script tag in index.html
declare const $3Dmol: {
  createViewer: (element: HTMLElement, config?: Record<string, unknown>) => unknown;
  [key: string]: unknown;
};

// Vite environment variables
interface ImportMetaEnv {
  readonly VITE_API_URL: string;
}

interface ImportMeta {
  readonly env: ImportMetaEnv;
}
