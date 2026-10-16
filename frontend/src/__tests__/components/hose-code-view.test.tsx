import { describe, it, expect, vi } from "vitest";
import { fireEvent, render, screen, within } from "@testing-library/react";
import { AppProvider } from "@/context/AppContext";
import HOSECodeView from "@/components/chem/HOSECodeView";

vi.mock("@/services/chemService", () => ({
  generateHOSECodes: vi.fn(async () => ["C-4;CC(//)", "C-4;C(//)"]),
}));

const renderView = () =>
  render(
    <AppProvider>
      <HOSECodeView />
    </AppProvider>
  );

describe("HOSECodeView (bento)", () => {
  it("lays out input, about and results panels", () => {
    renderView();
    expect(screen.getByRole("region", { name: "Structure input" })).toBeInTheDocument();
    expect(screen.getByRole("region", { name: "About HOSE codes" })).toBeInTheDocument();
    expect(screen.getByRole("region", { name: "Results" })).toHaveTextContent(/Enter a SMILES/);
  });

  it("lists one HOSE code per atom in the results panel", async () => {
    renderView();
    fireEvent.change(screen.getByLabelText(/SMILES/), { target: { value: "CCO" } });
    fireEvent.click(screen.getByRole("button", { name: /Generate HOSE codes/i }));
    const results = await screen.findByRole("region", { name: "2 HOSE codes" });
    expect(within(results).getByText("C-4;CC(//)")).toBeInTheDocument();
    expect(within(results).getByText("C-4;C(//)")).toBeInTheDocument();
  });
});
