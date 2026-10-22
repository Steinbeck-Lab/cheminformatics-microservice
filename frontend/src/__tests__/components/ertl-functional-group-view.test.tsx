import { describe, it, expect, vi } from "vitest";
import { fireEvent, render, screen, within } from "@testing-library/react";
import { AppProvider } from "@/context/AppContext";
import ErtlFunctionalGroupView from "@/components/chem/ErtlFunctionalGroupView";

vi.mock("@/services/chemService", () => ({
  generateFunctionalGroups: vi.fn(async () => [
    { description: "Carboxylic acid", atomIds: [1, 2, 3] },
    { description: "Ester", atomIds: [4, 5] },
  ]),
}));

const renderView = () =>
  render(
    <AppProvider>
      <ErtlFunctionalGroupView />
    </AppProvider>
  );

describe("ErtlFunctionalGroupView (bento)", () => {
  it("lays out input, about and results panels", () => {
    renderView();
    expect(screen.getByRole("region", { name: "Structure input" })).toBeInTheDocument();
    expect(screen.getByRole("region", { name: "About functional groups" })).toBeInTheDocument();
    expect(screen.getByRole("region", { name: "Results" })).toHaveTextContent(/Enter a SMILES/);
  });

  it("lists the groups and highlights one on click", async () => {
    renderView();
    fireEvent.change(screen.getByLabelText(/SMILES/), {
      target: { value: "CC(=O)OC1=CC=CC=C1C(=O)O" },
    });
    fireEvent.click(screen.getByRole("button", { name: /Detect functional groups/i }));
    const results = await screen.findByRole("region", { name: "2 functional groups" });
    const acid = within(results).getByRole("button", { name: /Carboxylic acid/ });
    expect(within(results).getByRole("button", { name: /Ester/ })).toBeInTheDocument();
    fireEvent.click(acid);
    expect(acid).toHaveAttribute("aria-pressed", "true");
  });
});
