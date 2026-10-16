import { beforeEach, describe, it, expect, vi } from "vitest";
import { fireEvent, render, screen, within } from "@testing-library/react";
import { AppProvider } from "@/context/AppContext";
import StructureVisualizerView from "@/components/depict/StructureVisualizerView";
import { lookupPubChem } from "@/services/chemService";

vi.mock("@/services/chemService", () => ({
  lookupPubChem: vi.fn(async (identifier: string) => ({
    success: true,
    input: identifier,
    input_type: "name",
    name: "Caffeine",
    canonical_smiles: "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
  })),
}));

// The viewers draw with images and 3Dmol; stand-ins show what they were given.
vi.mock("@/components/depict/MoleculeDepiction2D", () => ({
  default: ({ smiles, title, toolkit }: { smiles: string; title: string; toolkit: string }) => (
    <div data-testid="depiction-2d" data-toolkit={toolkit} data-smiles={smiles}>
      {title}
    </div>
  ),
}));
vi.mock("@/components/depict/MoleculeDepiction3D", () => ({
  default: ({ smiles, title, toolkit }: { smiles: string; title: string; toolkit: string }) => (
    <div data-testid="depiction-3d" data-toolkit={toolkit} data-smiles={smiles}>
      {title}
    </div>
  ),
}));

const renderView = () =>
  render(
    <AppProvider>
      <StructureVisualizerView />
    </AppProvider>
  );

const search = (value: string) => {
  fireEvent.change(screen.getByLabelText("Chemical identifier"), { target: { value } });
  fireEvent.click(screen.getByRole("button", { name: /Search/ }));
};

describe("StructureVisualizerView (bento)", () => {
  beforeEach(() => vi.clearAllMocks());

  it("lays out input, options and results panels", () => {
    renderView();
    const input = screen.getByRole("region", { name: "Structure input" });
    expect(within(input).getByLabelText("Chemical identifier")).toBeInTheDocument();
    const options = screen.getByRole("region", { name: "Options" });
    const tk2d = within(options).getByRole("radiogroup", { name: "2D toolkit" });
    expect(within(tk2d).getByRole("radio", { name: "RDKit" })).toHaveAttribute(
      "aria-checked",
      "true"
    );
    expect(within(options).getByRole("radiogroup", { name: "3D toolkit" })).toBeInTheDocument();
    expect(screen.getByRole("region", { name: "Results" })).toHaveTextContent(/Search PubChem/);
  });

  it("fills the field from an example", () => {
    renderView();
    fireEvent.click(screen.getByRole("button", { name: "CAS: 50-78-2" }));
    expect(screen.getByLabelText("Chemical identifier")).toHaveValue("50-78-2");
  });

  it("looks the identifier up and shows SMILES with 2D and 3D views", async () => {
    renderView();
    fireEvent.click(screen.getByRole("radio", { name: "CDK" }));
    search(" caffeine ");
    const results = await screen.findByRole("region", { name: "Structure" });
    expect(lookupPubChem).toHaveBeenCalledWith("caffeine");
    expect(within(results).getByText("Found: Caffeine")).toBeInTheDocument();
    expect(within(results).getAllByText("CN1C=NC2=C1C(=O)N(C(=O)N2C)C").length).toBeGreaterThan(0);
    expect(within(results).getByRole("button", { name: /Copy SMILES/ })).toBeInTheDocument();
    expect(within(results).getByTestId("depiction-2d")).toHaveAttribute("data-toolkit", "cdk");
    expect(within(results).getByTestId("depiction-3d")).toHaveAttribute(
      "data-toolkit",
      "openbabel"
    );
  });

  it("says so when PubChem has no structure", async () => {
    vi.mocked(lookupPubChem).mockResolvedValueOnce({
      success: false,
      input: "nothing",
      input_type: "name",
      canonical_smiles: "",
    });
    renderView();
    search("nothing");
    expect(await screen.findByText("No structure found in PubChem")).toBeInTheDocument();
    expect(screen.queryByTestId("depiction-2d")).toBeNull();
  });

  it("keeps a lookup error in the results panel", async () => {
    vi.mocked(lookupPubChem).mockRejectedValueOnce(new Error("Network down"));
    renderView();
    search("caffeine");
    const results = screen.getByRole("region", { name: "Results" });
    expect(await within(results).findByRole("alert")).toBeInTheDocument();
  });
});
