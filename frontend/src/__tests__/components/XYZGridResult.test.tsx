import { describe, it, expect, vi, beforeEach } from "vitest";
import { render, screen } from "@testing-library/react";
import userEvent from "@testing-library/user-event";
import XYZGridResult from "@/components/convert/XYZGridResult";
import type { XYZBatchConversionResult } from "@/types/api";

// Stub MoleculeDepiction2D — its real implementation hits an external API.
vi.mock("@/components/depict/MoleculeDepiction2D", () => ({
  default: ({ smiles }: { smiles: string }) => <div data-testid="depiction">{smiles}</div>,
}));

const fixture: XYZBatchConversionResult = {
  structures: [
    {
      index: 0,
      title: "water",
      success: true,
      error: null,
      canonicalsmiles: "[H]O[H]",
      inchi: "InChI=1S/H2O/h1H2",
      inchikey: "XLYOFNOQVPJJNP-UHFFFAOYSA-N",
      molblock: "BLOCK1\nM  END",
      method: "bond_orders",
      bond_orders_perceived: true,
      warnings: [],
    },
    {
      index: 1,
      title: "ferrocene",
      success: true,
      error: null,
      canonicalsmiles: "[Fe]",
      inchi: "InChI=1S/Fe",
      inchikey: "XEEYBQQBJWHFJM-UHFFFAOYSA-N",
      molblock: "BLOCK2\nM  END",
      method: "connectivity_only",
      bond_orders_perceived: false,
      warnings: ["RDKit DetermineBonds failed; fell back to connectivity-only perception"],
    },
    {
      index: 2,
      title: "bad",
      success: false,
      error: "unknown element Xx",
      canonicalsmiles: "",
      inchi: "",
      inchikey: "",
      molblock: "",
      method: "",
      bond_orders_perceived: false,
      warnings: [],
    },
  ],
  sdf: "BLOCK1\nM  END\n$$$$\nBLOCK2\nM  END\n$$$$\n",
  summary: {
    total: 3,
    successful: 2,
    failed: 1,
    bond_orders_count: 1,
    connectivity_only_count: 1,
  },
};

beforeEach(() => {
  // Stub URL.createObjectURL / revokeObjectURL — jsdom doesn't ship them.
  Object.assign(URL, {
    createObjectURL: vi.fn(() => "blob:fake"),
    revokeObjectURL: vi.fn(),
  });
});

describe("XYZGridResult", () => {
  it("renders one card per structure", () => {
    render(<XYZGridResult result={fixture} />);
    // 2 successful frames have depictions; the failed one does not.
    expect(screen.getAllByTestId("depiction")).toHaveLength(2);
    expect(screen.getByText("water")).toBeInTheDocument();
    expect(screen.getByText("ferrocene")).toBeInTheDocument();
    expect(screen.getByText("bad")).toBeInTheDocument();
  });

  it("badges connectivity-only frames", () => {
    render(<XYZGridResult result={fixture} />);
    // Multiple elements will contain "connectivity-only" (badge + summary span) — just confirm at least one exists.
    const els = screen.getAllByText(/connectivity-only/i);
    expect(els.length).toBeGreaterThanOrEqual(1);
  });

  it("shows error string for failed frames", () => {
    render(<XYZGridResult result={fixture} />);
    expect(screen.getByText(/unknown element Xx/)).toBeInTheDocument();
  });

  it("download SDF (all) creates a Blob with the full SDF", async () => {
    const user = userEvent.setup();
    render(<XYZGridResult result={fixture} />);
    const btn = screen.getByRole("button", { name: /download sdf/i });
    await user.click(btn);
    const calls = (URL.createObjectURL as unknown as ReturnType<typeof vi.fn>).mock.calls;
    expect(calls).toHaveLength(1);
    expect(calls[0][0]).toBeInstanceOf(Blob);
  });

  it("displays summary line", () => {
    render(<XYZGridResult result={fixture} />);
    // Summary text is spread across child spans; match against the container's full textContent.
    const summaryDiv = screen.getByTestId("xyz-summary");
    expect(summaryDiv.textContent).toMatch(/2\s+of\s+3/i);
    expect(summaryDiv.textContent).toMatch(/1\s+connectivity-only/i);
    expect(summaryDiv.textContent).toMatch(/1\s+failed/i);
  });
});
