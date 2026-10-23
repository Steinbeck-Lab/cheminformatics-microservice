import { describe, it, expect, vi, beforeEach } from "vitest";
import { act, fireEvent, render, screen, within } from "@testing-library/react";
import { AppProvider } from "@/context/AppContext";
import FormatConversionView from "@/components/convert/FormatConversionView";
import convertService from "@/services/convertService";
import type { XYZBatchConversionResult } from "@/types/api";

const MOLBLOCK = `
  CDK     08302311362D

  2  1  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
M  END`;

const frame = (index: number, title: string, smiles: string) => ({
  index,
  title,
  success: true,
  error: null,
  canonicalsmiles: smiles,
  inchi: `InChI=1S/${title}`,
  inchikey: `${title.toUpperCase()}-KEY`,
  molblock: "BLOCK\nM  END",
  method: "bond_orders",
  bond_orders_perceived: true,
  warnings: [],
});

const XYZ_BATCH: XYZBatchConversionResult = {
  structures: [frame(0, "water", "O"), frame(1, "methane", "C")],
  sdf: "BLOCK\nM  END\n$$$$\n",
  summary: {
    total: 2,
    successful: 2,
    failed: 0,
    bond_orders_count: 2,
    connectivity_only_count: 0,
  },
} as XYZBatchConversionResult;

vi.mock("@/services/convertService", () => {
  const svc = {
    detectInputFormat: vi.fn(async () => ({ detected_format: "smiles", confidence: "high" })),
    generate2DCoordinates: vi.fn(async () => "MOL"),
    generateSMILES: vi.fn(async () => "CCO"),
    generateCanonicalSMILES: vi.fn(async () => "CCO"),
    generateCXSMILES: vi.fn(async () => "CCO |cx|"),
    generateInChI: vi.fn(async () => "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3"),
    generateInChIKey: vi.fn(async () => '"LFQSCWFLJHTTHZ-UHFFFAOYSA-N"'),
    generateSELFIES: vi.fn(async () => "[C][C][O]"),
    generateSMARTS: vi.fn(async () => "[#6][#6][#8]"),
    generateIUPACName: vi.fn(async () => ({
      name: "ethanol",
      tier: "pin_verified",
      is_pin: true,
      verified: "opsin",
    })),
    molblockToSMILES: vi.fn(async () => "CO"),
    convertCDXToMol: vi.fn(async () => "MOL"),
    convertXYZ: vi.fn(async () => XYZ_BATCH),
  };
  return { default: svc, ...svc };
});

vi.mock("@/components/depict/MoleculeDepiction2D", () => ({
  default: ({ smiles }: { smiles: string }) => <div data-testid="depiction">{smiles}</div>,
}));

const renderView = () =>
  render(
    <AppProvider>
      <FormatConversionView />
    </AppProvider>
  );

describe("FormatConversionView (bento)", () => {
  beforeEach(() => vi.clearAllMocks());

  it("lays out input, options and results panels", () => {
    renderView();
    const input = screen.getByRole("region", { name: "Structure input" });
    expect(within(input).getByRole("combobox", { name: "Input format" })).toHaveValue("smiles");
    expect(within(input).getByLabelText(/^SMILES/)).toBeInTheDocument();
    const options = screen.getByRole("region", { name: "Options" });
    expect(within(options).getByRole("combobox", { name: "Output format" })).toHaveValue(
      "canonicalsmiles"
    );
    expect(within(options).getByRole("radiogroup", { name: "Toolkit" })).toBeInTheDocument();
    expect(screen.getByRole("region", { name: "Results" })).toHaveTextContent(/Enter a structure/);
  });

  it("converts SMILES to an InChI Key and shows it with a copy button", async () => {
    renderView();
    fireEvent.change(screen.getByLabelText(/^SMILES/), { target: { value: "CCO" } });
    fireEvent.change(screen.getByRole("combobox", { name: "Output format" }), {
      target: { value: "inchikey" },
    });
    fireEvent.click(screen.getByRole("button", { name: /Convert format/i }));
    const results = await screen.findByRole("region", { name: "InChI Key" });
    expect(convertService.generateInChIKey).toHaveBeenCalledWith("CCO", "cdk", {
      autoDetect: false,
    });
    // surrounding quotes from the API are stripped
    expect(within(results).getByText("LFQSCWFLJHTTHZ-UHFFFAOYSA-N")).toBeInTheDocument();
    expect(within(results).getByRole("button", { name: /Copy result/i })).toBeInTheDocument();
    expect(within(results).getByTestId("depiction")).toHaveTextContent("CCO");
  });

  it("names SMILES with Orthonym at the chosen tier and shows the tier", async () => {
    renderView();
    fireEvent.change(screen.getByLabelText(/^SMILES/), { target: { value: "CCO" } });
    fireEvent.change(screen.getByRole("combobox", { name: "Output format" }), {
      target: { value: "iupac" },
    });
    expect(screen.queryByRole("radiogroup", { name: "Toolkit" })).not.toBeInTheDocument();
    fireEvent.click(screen.getByRole("radio", { name: "Valid" }));
    fireEvent.click(screen.getByRole("button", { name: /Convert format/i }));
    const results = await screen.findByRole("region", { name: "IUPAC Name" });
    expect(convertService.generateIUPACName).toHaveBeenCalledWith("CCO", "valid", {
      autoDetect: false,
    });
    expect(within(results).getByText("ethanol", { selector: "pre" })).toBeInTheDocument();
    expect(results).toHaveTextContent("with Orthonym (pin verified)");
  });

  it("converts a pasted MOL block to SMILES, keeping the two-space terminator", async () => {
    renderView();
    fireEvent.change(screen.getByRole("combobox", { name: "Input format" }), {
      target: { value: "molsdf" },
    });
    expect(screen.getByRole("combobox", { name: "Output format" })).toBeDisabled();
    expect(screen.getByRole("button", { name: /Upload MOL\/SDF file/i })).toBeInTheDocument();
    fireEvent.change(screen.getByLabelText(/MOL\/SDF block/), { target: { value: MOLBLOCK } });
    fireEvent.click(screen.getByRole("button", { name: /Convert format/i }));
    const results = await screen.findByRole("region", { name: "SMILES" });
    const [sent, toolkit] = vi.mocked(convertService.molblockToSMILES).mock.calls[0]!;
    expect(sent).toMatch(/\nM {2}END\n$/);
    expect(toolkit).toBe("cdk");
    expect(within(results).getByText("CO", { selector: "pre" })).toBeInTheDocument();
  });

  it("shows a multi-frame XYZ batch as a grid", async () => {
    renderView();
    fireEvent.change(screen.getByRole("combobox", { name: "Input format" }), {
      target: { value: "xyz" },
    });
    const options = screen.getByRole("region", { name: "Options" });
    // CDK cannot perceive bonds from XYZ, so the toolkit moves to RDKit
    expect(within(options).getByRole("radio", { name: "RDKit" })).toHaveAttribute(
      "aria-checked",
      "true"
    );
    expect(within(options).queryByRole("radio", { name: "CDK" })).not.toBeInTheDocument();
    fireEvent.click(within(options).getByRole("radio", { name: "OpenBabel" }));
    expect(within(options).getByLabelText("Net charge")).toBeDisabled();
    fireEvent.click(within(options).getByRole("radio", { name: "RDKit" }));
    expect(within(options).getByLabelText("Net charge")).toBeEnabled();
    fireEvent.change(screen.getByLabelText(/XYZ block/), {
      target: { value: "3\nwater\nO 0 0 0" },
    });
    await act(async () => {
      fireEvent.click(screen.getByRole("button", { name: /Convert format/i }));
    });
    const results = await screen.findByRole("region", { name: "2 structures" });
    expect(convertService.convertXYZ).toHaveBeenCalledWith("3\nwater\nO 0 0 0", {
      charge: 0,
      useHueckel: false,
      toolkit: "rdkit",
    });
    expect(within(results).getByTestId("xyz-summary")).toHaveTextContent(/2\s+of\s+2/);
    expect(within(results).queryByText(/Enter a structure/)).not.toBeInTheDocument();
  });
});
