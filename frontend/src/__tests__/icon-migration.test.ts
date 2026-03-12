import { describe, it, expect } from "vitest";
import * as fs from "fs";
import * as path from "path";

describe("Icon migration", () => {
  const packageJsonPath = path.resolve(__dirname, "../../package.json");
  const srcDir = path.resolve(__dirname, "..");

  it("no react-icons in package.json dependencies", () => {
    const packageJson = JSON.parse(fs.readFileSync(packageJsonPath, "utf8"));
    const allDeps = {
      ...packageJson.dependencies,
      ...packageJson.devDependencies,
    };
    expect(allDeps).not.toHaveProperty("react-icons");
  });

  it("no @fortawesome in package.json dependencies", () => {
    const packageJson = JSON.parse(fs.readFileSync(packageJsonPath, "utf8"));
    const allDeps = {
      ...packageJson.dependencies,
      ...packageJson.devDependencies,
    };
    expect(allDeps).not.toHaveProperty("@fortawesome/react-fontawesome");
    expect(allDeps).not.toHaveProperty("@fortawesome/free-solid-svg-icons");
  });

  it("lucide-react is installed as a dependency", () => {
    const packageJson = JSON.parse(fs.readFileSync(packageJsonPath, "utf8"));
    expect(packageJson.dependencies).toHaveProperty("lucide-react");
  });

  it("Header.tsx imports from lucide-react", () => {
    const headerPath = path.join(srcDir, "components", "common", "Header.tsx");
    const content = fs.readFileSync(headerPath, "utf8");
    expect(content).toContain('from "lucide-react"');
    expect(content).not.toContain("react-icons");
  });

  it("HomePage.tsx imports from lucide-react", () => {
    const homePath = path.join(srcDir, "pages", "HomePage.tsx");
    const content = fs.readFileSync(homePath, "utf8");
    expect(content).toContain('from "lucide-react"');
    expect(content).not.toContain("react-icons");
  });

  it("Footer.tsx imports from lucide-react", () => {
    const footerPath = path.join(srcDir, "components", "common", "Footer.tsx");
    const content = fs.readFileSync(footerPath, "utf8");
    expect(content).toContain('from "lucide-react"');
    expect(content).not.toContain("react-icons");
  });
});
