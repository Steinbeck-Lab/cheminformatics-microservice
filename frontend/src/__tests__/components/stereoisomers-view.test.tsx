import { describe, it, expect, vi } from "vitest";
import { fireEvent, render, screen, within } from "@testing-library/react";
import { AppProvider } from "@/context/AppContext";
import StereoisomersView from "@/components/chem/StereoisomersView";

vi.mock("@/services/chemService", () => ({
  generateStereoisomers: vi.fn(async () => ["C[C@H](N)O", "C[C@@H](N)O"]),
}));

const renderView = () =>
  render(
    <AppProvider>
      <StereoisomersView />
    </AppProvider>
  );

describe("StereoisomersView (bento)", () => {
  it("lays out input, about and results panels", () => {
    renderView();
    expect(screen.getByRole("region", { name: "Structure input" })).toBeInTheDocument();
    expect(screen.getByRole("region", { name: "About stereoisomers" })).toBeInTheDocument();
    expect(screen.getByRole("region", { name: "Results" })).toHaveTextContent(/Enter a SMILES/);
  });

  it("shows each stereoisomer in the results panel", async () => {
    renderView();
    fireEvent.change(screen.getByLabelText(/SMILES/), { target: { value: "CC(N)O" } });
    fireEvent.click(screen.getByRole("button", { name: /Generate stereoisomers/i }));
    const results = await screen.findByRole("region", { name: "2 stereoisomers" });
    expect(within(results).getByText("Stereoisomer 1")).toBeInTheDocument();
    expect(within(results).getByText("Stereoisomer 2")).toBeInTheDocument();
  });
});
