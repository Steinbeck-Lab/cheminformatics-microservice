import { describe, it, expect, vi } from "vitest";
import { fireEvent, render, screen, within } from "@testing-library/react";
import { AppProvider } from "@/context/AppContext";
import PubChemLookupView from "@/components/chem/PubChemLookupView";

vi.mock("@/services/chemService", () => ({
  lookupPubChem: vi.fn(async () => ({
    success: true,
    input: "aspirin",
    input_type: "name",
    canonical_smiles: "CC(=O)OC1=CC=CC=C1C(=O)O",
    cids: [2244],
    pubchem_links: ["https://pubchem.ncbi.nlm.nih.gov/compound/2244"],
  })),
}));

const renderView = () =>
  render(
    <AppProvider>
      <PubChemLookupView />
    </AppProvider>
  );

describe("PubChemLookupView (bento)", () => {
  it("lays out input, about and results panels", () => {
    renderView();
    expect(screen.getByRole("region", { name: "Structure input" })).toBeInTheDocument();
    expect(screen.getByRole("region", { name: "About PubChem lookup" })).toBeInTheDocument();
    expect(screen.getByRole("region", { name: "Results" })).toHaveTextContent(/Enter a name/);
  });

  it("shows the PubChem match in the results panel", async () => {
    renderView();
    fireEvent.click(screen.getByRole("button", { name: "Chemical Name" }));
    expect(screen.getByLabelText(/Chemical identifier/)).toHaveValue("aspirin");
    fireEvent.click(screen.getByRole("button", { name: /Find structure/i }));
    const results = await screen.findByRole("region", { name: "PubChem match" });
    expect(within(results).getByText(/PubChem CID: 2244/)).toBeInTheDocument();
    expect(within(results).getByRole("link", { name: /CID 2244/ })).toHaveAttribute(
      "href",
      "https://pubchem.ncbi.nlm.nih.gov/compound/2244"
    );
  });
});
