/// <reference types="vitest/config" />
import { defineConfig } from "vite";
import react from "@vitejs/plugin-react-swc";
import tailwindcss from "@tailwindcss/vite";

export default defineConfig({
  plugins: [
    tailwindcss(),
    react({
      jsxRuntime: "automatic",
    }),
  ],
  server: {
    port: 3000,
    proxy: {
      "/v1": {
        target: "http://localhost:8000",
        changeOrigin: true,
      },
      "/latest": {
        target: "http://localhost:8000",
        changeOrigin: true,
      },
    },
  },
  // optimizeDeps.esbuildOptions.loader REMOVED -- no longer needed
  build: {
    sourcemap: "hidden",
    outDir: "dist",
  },
  test: {
    globals: true,
    environment: "jsdom",
    setupFiles: "./src/__tests__/setup.ts",
    css: true,
    coverage: {
      provider: "v8",
      reporter: ["text", "lcov"],
    },
  },
});
