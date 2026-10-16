import { describe, it, expect } from "vitest";
import { DEFAULT_VERSION, INCHI_VERSIONS } from "@/utils/inchiUtils";

function cfg(key: string) {
  const c = INCHI_VERSIONS[key];
  if (!c) throw new Error(`No InChI version ${key}`);
  return c;
}

// The IUPAC InChI-Web-Demo site renamed its modules (pages/inchi_versions.json):
// "latest" is now 1.07.5 (inchi-web-1075.js) and the molecular-inorganics build is
// inchi-web-dev-moin.js. The old file names return 404.
describe("InChI module config", () => {
  it("defaults to InChI 1.07.5 from the live file", () => {
    const def = cfg(DEFAULT_VERSION);
    expect(def.default).toBe(true);
    expect(def.scriptSrc).toMatch(/\/inchi\/inchi-web-1075\.js$/);
    expect(def.moduleName).toBe("inchiModule1075");
    expect(def.label).toContain("1.07.5");
  });

  it("loads the molecular-inorganics build from its new file", () => {
    expect(cfg("Latest-MoIn").scriptSrc).toMatch(/\/inchi\/inchi-web-dev-moin\.js$/);
    expect(cfg("Latest-MoIn").moduleName).toBe("inchiModuleDevMoIn");
  });

  it("no longer points at the removed latest files", () => {
    for (const x of Object.values(INCHI_VERSIONS))
      expect(x.scriptSrc).not.toMatch(/inchi-web-latest/);
  });

  it("offers the enhanced-stereochemistry and polymer builds from their live files", () => {
    expect(cfg("Dev-EnhancedStereo").scriptSrc).toMatch(
      /\/inchi\/inchi-web-dev-enhanced-stereo\.js$/
    );
    expect(cfg("Dev-EnhancedStereo").moduleName).toBe("inchiModuleDevEnhancedStereo");
    expect(cfg("Polymer").scriptSrc).toMatch(/\/inchi\/inchi-web-polymer-support\.js$/);
    expect(cfg("Polymer").moduleName).toBe("inchiModulePolymerSupport");
  });

  it.each([
    ["Dev", "inchi-web-dev", "inchiModuleDev"],
    ["Dev-NoMetalH", "inchi-web-no-metal-h", "inchiModuleNoMetalH"],
    [
      "Dev-ExplicitZeroValence",
      "inchi-web-explicit-zero-valence",
      "inchiModuleExplicitZeroValence",
    ],
  ])("offers the %s build from its live file", (key, file, mod) => {
    expect(cfg(key).scriptSrc).toMatch(new RegExp(`/inchi/${file}\\.js$`));
    expect(cfg(key).moduleName).toBe(mod);
  });
});
