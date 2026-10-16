import { afterEach, beforeEach, describe, it, expect, vi } from "vitest";
import { fireEvent, render, screen, within } from "@testing-library/react";
import { AppProvider } from "@/context/AppContext";
import StandardizeView from "@/components/chem/StandardizeView";
import { generate2DCoordinates } from "@/services/convertService";

const MOLBLOCK = `
  RDKit          2D

  3  2  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0
    1.2990    0.7500    0.0000 C   0  0
    2.5981   -0.0000    0.0000 O   0  0
  1  2  1  0
  2  3  1  0
M  END
`;

vi.mock("@/services/convertService", () => ({
  generate2DCoordinates: vi.fn(async () => MOLBLOCK),
}));

const RESULT = {
  standardized_mol: MOLBLOCK,
  canonical_smiles: "CCO",
  inchi: "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3",
  inchikey: "LFQSCWFLJHTTHZ-UHFFFAOYSA-N",
  formula: "C2H6O",
};

const renderView = () =>
  render(
    <AppProvider>
      <StandardizeView />
    </AppProvider>
  );

describe("StandardizeView (bento)", () => {
  beforeEach(() => {
    vi.stubGlobal(
      "fetch",
      vi.fn(async () => new Response(JSON.stringify(RESULT), { status: 200 }))
    );
  });
  afterEach(() => vi.unstubAllGlobals());

  it("lays out input, about and results panels", () => {
    renderView();
    const input = screen.getByRole("region", { name: "Structure input" });
    expect(within(input).getByRole("radiogroup", { name: "Input method" })).toBeInTheDocument();
    expect(within(input).getByLabelText(/Molblock/)).toBeInTheDocument();
    expect(
      within(input).getByRole("button", { name: /Upload MOL\/SDF file/i })
    ).toBeInTheDocument();
    expect(screen.getByRole("region", { name: "About standardization" })).toBeInTheDocument();
    expect(screen.getByRole("region", { name: "Results" })).toHaveTextContent(
      /molblock or a SMILES/
    );
  });

  it("standardizes a pasted molblock and shows the identifiers", async () => {
    renderView();
    fireEvent.change(screen.getByLabelText(/Molblock/), { target: { value: MOLBLOCK } });
    fireEvent.click(screen.getByRole("button", { name: /Standardize molecule/i }));
    const results = await screen.findByRole("region", { name: "Standardized structure" });
    expect(within(results).getByText("LFQSCWFLJHTTHZ-UHFFFAOYSA-N")).toBeInTheDocument();
    expect(within(results).getByRole("button", { name: /Copy molblock/i })).toBeInTheDocument();
    const extra = screen.getByRole("region", { name: "Additional properties" });
    expect(within(extra).getByText("C2H6O")).toBeInTheDocument();
  });

  it("converts a SMILES string to a molblock first", async () => {
    renderView();
    fireEvent.click(screen.getByRole("radio", { name: "SMILES" }));
    fireEvent.change(screen.getByLabelText(/SMILES/), { target: { value: "CCO" } });
    fireEvent.click(screen.getByRole("button", { name: /Standardize molecule/i }));
    await screen.findByRole("region", { name: "Standardized structure" });
    expect(generate2DCoordinates).toHaveBeenCalledWith("CCO", "rdkit");
  });
});
