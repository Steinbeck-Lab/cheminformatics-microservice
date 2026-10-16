/// <reference types="vitest/config" />
import path from "path";
import { defineConfig } from "vite";
import react from "@vitejs/plugin-react-swc";
import tailwindcss from "@tailwindcss/vite";

const VENDOR_REACT = new Set([
  "react",
  "react-dom",
  "scheduler",
  "react-router",
  "react-router-dom",
]);
const VENDOR_MOTION = new Set(["motion", "framer-motion", "motion-dom", "motion-utils"]);
const VENDOR_UI = new Set(["lucide-react", "class-variance-authority", "clsx", "tailwind-merge"]);

export default defineConfig({
  plugins: [
    tailwindcss(),
    react({
      jsxRuntime: "automatic",
    }),
  ],
  resolve: {
    alias: {
      "@": path.resolve(__dirname, "./src"),
    },
  },
  server: {
    port: 3000,
    proxy: {
      "/v1": {
        target: "http://localhost:8000",
        changeOrigin: true,
      },
      "/latest": {
        // VITE_PROXY_TARGET=https://dev.api.naturalproducts.net lets a local UI use the dev API
        // without CORS (pair it with VITE_API_URL=/latest)
        target: process.env.VITE_PROXY_TARGET || "http://localhost:8000",
        changeOrigin: true,
      },
    },
  },
  // optimizeDeps.esbuildOptions.loader REMOVED -- no longer needed
  build: {
    sourcemap: "hidden",
    outDir: "dist",
    rollupOptions: {
      output: {
        // by node_modules path, so deep imports (react-dom/client, motion/react) join their group
        manualChunks(id) {
          const pkg = id.match(/.*node_modules\/((?:@[^/]+\/)?[^/]+)\//)?.[1];
          if (!pkg) return undefined;
          if (VENDOR_REACT.has(pkg)) return "vendor-react";
          if (VENDOR_MOTION.has(pkg)) return "vendor-motion";
          if (VENDOR_UI.has(pkg)) return "vendor-ui";
          return undefined;
        },
      },
    },
  },
  test: {
    globals: true,
    environment: "jsdom",
    // Node 25+ ships a global localStorage that shadows jsdom's and lacks its methods
    execArgv: ["--no-experimental-webstorage"],
    setupFiles: "./src/__tests__/setup.ts",
    css: true,
    alias: {
      "@": path.resolve(__dirname, "./src"),
    },
    coverage: {
      provider: "v8",
      reporter: ["text", "lcov"],
    },
  },
});
