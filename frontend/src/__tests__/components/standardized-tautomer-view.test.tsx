import { describe, it, expect, vi } from "vitest";
import { fireEvent, render, screen, within } from "@testing-library/react";
import { AppProvider } from "@/context/AppContext";
import StandardizedTautomerView from "@/components/chem/StandardizedTautomerView";

vi.mock("@/services/chemService", () => ({
  generateStandardizedTautomer: vi.fn(async () => "CCC(C)=O"),
}));

const renderView = () =>
  render(
    <AppProvider>
      <StandardizedTautomerView />
    </AppProvider>
  );

describe("StandardizedTautomerView (bento)", () => {
  it("lays out input, about and results panels", () => {
    renderView();
    expect(screen.getByRole("region", { name: "Structure input" })).toBeInTheDocument();
    expect(screen.getByRole("region", { name: "About tautomers" })).toBeInTheDocument();
    expect(screen.getByRole("region", { name: "Results" })).toHaveTextContent(/Enter a SMILES/);
  });

  it("shows the original and the standardized tautomer in the results panel", async () => {
    renderView();
    fireEvent.change(screen.getByLabelText(/SMILES/), { target: { value: "CC(O)=CC" } });
    fireEvent.click(screen.getByRole("button", { name: /Standardize tautomer/i }));
    const results = await screen.findByRole("region", { name: "Standardized tautomer" });
    expect(within(results).getAllByText("CCC(C)=O").length).toBeGreaterThan(0);
    expect(within(results).getAllByText("CC(O)=CC").length).toBeGreaterThan(0);
    expect(within(results).getByText("Standardized SMILES")).toBeInTheDocument();
  });
});
