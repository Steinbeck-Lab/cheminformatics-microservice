import { describe, it, expect, vi } from "vitest";
import { fireEvent, render, screen, within } from "@testing-library/react";
import { AppProvider } from "@/context/AppContext";
import NPlikenessView from "@/components/chem/NPlikenessView";

vi.mock("@/services/chemService", () => ({
  calculateNPLikeness: vi.fn(async () => 1.234),
}));

const renderView = () =>
  render(
    <AppProvider>
      <NPlikenessView />
    </AppProvider>
  );

describe("NPlikenessView (bento)", () => {
  it("lays out input, about and results panels", () => {
    renderView();
    expect(screen.getByRole("region", { name: "Structure input" })).toBeInTheDocument();
    expect(screen.getByRole("region", { name: "About NP-likeness" })).toBeInTheDocument();
    expect(screen.getByRole("region", { name: "Results" })).toHaveTextContent(/Enter a SMILES/);
  });

  it("shows the score and its reading in the results panel", async () => {
    renderView();
    fireEvent.change(screen.getByLabelText(/SMILES/), { target: { value: "CCO" } });
    fireEvent.click(screen.getByRole("button", { name: /Calculate NP score/i }));
    const results = await screen.findByRole("region", { name: "NP-likeness score" });
    expect(within(results).getByText("1.23")).toBeInTheDocument();
    expect(within(results).getByText(/Very high natural product likeness/)).toBeInTheDocument();
  });
});
