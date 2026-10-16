import { beforeEach, describe, it, expect, vi } from "vitest";
import { fireEvent, render, screen, waitFor, within } from "@testing-library/react";
import { AppProvider } from "@/context/AppContext";
import InChIView from "@/components/tools/InChIView";
import { convertMolfileToInchi } from "@/utils/inchiUtils";

const MOLFILE = `
  Ketcher  1016261212D 1   1.00000     0.00000     0

  2  1  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0     0  0
M  END
`;

const ketcher = vi.hoisted(() => ({
  executeCommand: vi.fn(),
  resetEditor: vi.fn(),
}));

vi.mock("@/hooks/useKetcherEditor", () => ({
  useKetcherEditor: () => ({
    ketcherFrame: { current: null },
    isEditorReady: true,
    executeCommand: ketcher.executeCommand,
    resetEditor: ketcher.resetEditor,
  }),
}));

vi.mock("@/utils/inchiUtils", () => ({
  VERSION_OPTIONS: [
    ["1.06", "InChI 1.06"],
    ["Latest", "InChI Latest"],
    ["Latest-MoIn", "InChI latest with Molecular inorganics"],
    ["Dev-EnhancedStereo", "InChI dev with Enhanced Stereochemistry"],
    ["Polymer", "InChI with Polymer Support"],
    ["Dev", "InChI dev"],
    ["Dev-NoMetalH", "InChI dev, no metal hydrogens"],
    ["Dev-ExplicitZeroValence", "InChI dev, explicit zero valence"],
  ],
  DEFAULT_VERSION: "Latest",
  loadInchiModule: vi.fn(async () => undefined),
  convertMolfileToInchi: vi.fn(async () => ({
    return_code: 0,
    inchi: "InChI=1S/CH4O/c1-2/h2H,1H3",
    auxinfo: "AuxInfo=1/0/N:1,2/rA:2nCO/rB:s1;/rC:0,0,0;1,0,0;",
    message: "",
  })),
  generateInchiKey: vi.fn(async () => ({
    return_code: 0,
    inchikey: "OKKJLVBELUTLKV-UHFFFAOYSA-N",
  })),
  convertInchiToMolfile: vi.fn(async () => ({ return_code: 0, molfile: MOLFILE })),
  convertAuxinfoToMolfile: vi.fn(async () => ({ return_code: 0, molfile: MOLFILE })),
}));

const renderView = () =>
  render(
    <AppProvider>
      <InChIView />
    </AppProvider>
  );

const generateButton = () => screen.getByRole("button", { name: /Generate InChI/i });

describe("InChIView (bento)", () => {
  beforeEach(() => {
    vi.clearAllMocks();
    ketcher.executeCommand.mockImplementation(async (cmd: string) =>
      cmd === "getMolfile" ? MOLFILE : undefined
    );
  });

  it("lays out input, options, editor and results panels", async () => {
    renderView();
    const input = screen.getByRole("region", { name: "Structure input" });
    expect(within(input).getByRole("radiogroup", { name: "Input method" })).toBeInTheDocument();
    const options = screen.getByRole("region", { name: "Options" });
    expect(within(options).getByRole("combobox", { name: "InChI version" })).toHaveValue("Latest");
    expect(within(options).getByRole("switch", { name: /Mobile H perception/ })).toHaveAttribute(
      "aria-checked",
      "true"
    );
    const editor = screen.getByRole("region", { name: "Structure editor" });
    expect(within(editor).getByTitle("Ketcher Editor")).toHaveAttribute(
      "src",
      "/standalone/index.html"
    );
    expect(screen.getByRole("region", { name: "Results" })).toHaveTextContent(/generate its InChI/);
    await waitFor(() => expect(generateButton()).toBeEnabled());
  });

  it("generates InChI, key and AuxInfo from the drawn structure", async () => {
    renderView();
    await waitFor(() => expect(generateButton()).toBeEnabled());
    fireEvent.click(generateButton());
    const results = await screen.findByRole("region", { name: "InChI" });
    expect(ketcher.executeCommand).toHaveBeenCalledWith("getMolfile");
    const [sent, opts, version] = vi.mocked(convertMolfileToInchi).mock.calls[0]!;
    expect(sent).toContain("M  END");
    expect(opts).toBe("");
    expect(version).toBe("Latest");
    expect(within(results).getByText("InChI=1S/CH4O/c1-2/h2H,1H3")).toBeInTheDocument();
    expect(await within(results).findByText("OKKJLVBELUTLKV-UHFFFAOYSA-N")).toBeInTheDocument();
    expect(within(results).getByRole("button", { name: "Copy InChIKey" })).toBeInTheDocument();
  });

  it("passes changed options to the converter", async () => {
    renderView();
    await waitFor(() => expect(generateButton()).toBeEnabled());
    fireEvent.click(screen.getByRole("switch", { name: /Mobile H perception/ }));
    fireEvent.click(generateButton());
    await screen.findByRole("region", { name: "InChI" });
    expect(vi.mocked(convertMolfileToInchi).mock.calls[0]![1]).toBe("-FixedH");
  });

  it("loads a SMILES string into the editor", async () => {
    renderView();
    fireEvent.click(screen.getByRole("radio", { name: "SMILES" }));
    fireEvent.change(screen.getByLabelText(/SMILES to edit/), { target: { value: "CO" } });
    fireEvent.click(screen.getByRole("button", { name: /Load structure/i }));
    await waitFor(() => expect(ketcher.executeCommand).toHaveBeenCalledWith("setMolecule", ["CO"]));
  });

  it("keeps an error visible when nothing is drawn", async () => {
    ketcher.executeCommand.mockImplementation(async () => "");
    renderView();
    await waitFor(() => expect(generateButton()).toBeEnabled());
    fireEvent.click(generateButton());
    const results = screen.getByRole("region", { name: "Results" });
    expect(await within(results).findByText(/No structure drawn/)).toBeInTheDocument();
  });

  it("says the editor is empty when it returns a molfile with no atoms", async () => {
    // what Ketcher really returns for an empty canvas: a header and a 0-atom counts line
    const empty = "\n  Ketcher  0101010000\n\n  0  0  0  0  0  0  0  0  0  0999 V2000\nM  END\n";
    ketcher.executeCommand.mockImplementation(async (cmd: string) =>
      cmd === "getMolfile" ? empty : undefined
    );
    renderView();
    await waitFor(() => expect(generateButton()).toBeEnabled());
    fireEvent.click(generateButton());
    const results = screen.getByRole("region", { name: "Results" });
    expect(await within(results).findByText(/No structure drawn/)).toBeInTheDocument();
    expect(vi.mocked(convertMolfileToInchi)).not.toHaveBeenCalled();
  });

  it("sends a V3000 molfile and the enhanced-stereo option to the enhanced-stereo build", async () => {
    renderView();
    await waitFor(() => expect(generateButton()).toBeEnabled());
    fireEvent.change(screen.getByRole("combobox", { name: "InChI version" }), {
      target: { value: "Dev-EnhancedStereo" },
    });
    fireEvent.click(await screen.findByRole("switch", { name: /Enhanced stereochemistry/ }));
    fireEvent.click(generateButton());
    await screen.findByRole("region", { name: "InChI" });
    expect(ketcher.executeCommand).toHaveBeenCalledWith("getMolfile", ["v3000"]);
    const [, opts, version] = vi.mocked(convertMolfileToInchi).mock.calls[0]!;
    expect(version).toBe("Dev-EnhancedStereo");
    expect(opts).toContain("-EnhancedStereochemistry");
  });

  it("switches polymer handling and NPZz on for the polymer build", async () => {
    renderView();
    await waitFor(() => expect(generateButton()).toBeEnabled());
    fireEvent.change(screen.getByRole("combobox", { name: "InChI version" }), {
      target: { value: "Polymer" },
    });
    // the new build loads first; the button is disabled meanwhile
    await waitFor(() => expect(generateButton()).toBeEnabled());
    fireEvent.click(generateButton());
    await screen.findByRole("region", { name: "InChI" });
    const [, opts, version] = vi.mocked(convertMolfileToInchi).mock.calls[0]!;
    expect(version).toBe("Polymer");
    expect(opts).toContain("-Polymers");
    expect(opts).toContain("-NPZz");
  });

  it.each([
    ["Dev", false],
    ["Dev-NoMetalH", true],
    ["Dev-ExplicitZeroValence", true],
  ])("the %s build sends MolecularInorganics by default: %s", async (version, moin) => {
    renderView();
    await waitFor(() => expect(generateButton()).toBeEnabled());
    fireEvent.change(screen.getByRole("combobox", { name: "InChI version" }), {
      target: { value: version },
    });
    await waitFor(() => expect(generateButton()).toBeEnabled());
    fireEvent.click(generateButton());
    await screen.findByRole("region", { name: "InChI" });
    const [, opts, sent] = vi.mocked(convertMolfileToInchi).mock.calls[0]!;
    expect(sent).toBe(version);
    expect(opts.includes("-MolecularInorganics")).toBe(moin);
    expect(opts).not.toContain("-NPZz");
    expect(screen.queryByRole("switch", { name: /Molecular inorganics/ }) !== null).toBe(moin);
  });
});
