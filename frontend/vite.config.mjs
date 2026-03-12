import { defineConfig } from "vite";
import react from "@vitejs/plugin-react-swc";

export default defineConfig({
  plugins: [
    // Enable JSX handling for .js files (CRA convention)
    react({
      jsxRuntime: "automatic",
    }),
    // Custom plugin to let Vite know .js files may contain JSX
    {
      name: "treat-js-as-jsx",
      async transform(code, id) {
        if (!id.match(/src\/.*\.js$/)) return null;
        // Use esbuild to transform JSX in .js files
        const { transformWithEsbuild } = await import("vite");
        return transformWithEsbuild(code, id + "x", {
          loader: "jsx",
          jsx: "automatic",
        });
      },
    },
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
  optimizeDeps: {
    esbuildOptions: {
      loader: {
        ".js": "jsx",
      },
    },
  },
  build: {
    sourcemap: "hidden",
    outDir: "dist",
  },
});
