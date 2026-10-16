import { describe, it, expect, vi } from "vitest";
import { fireEvent, render, screen, within } from "@testing-library/react";
import { AppProvider } from "@/context/AppContext";
import StructureErrorView from "@/components/chem/StructureErrorView";

vi.mock("@/services/chemService", () => ({
  checkStructureErrors: vi.fn(async () => ({
    original: { smi: "C1=CC=CC=C1", messages: ["Valence problem on atom 2"] },
    standardized: { smi: "c1ccccc1", messages: ["No Errors Found"] },
  })),
}));

const renderView = () =>
  render(
    <AppProvider>
      <StructureErrorView />
    </AppProvider>
  );

describe("StructureErrorView (bento)", () => {
  it("lays out input, about and results panels", () => {
    renderView();
    expect(screen.getByRole("region", { name: "Structure input" })).toBeInTheDocument();
    expect(screen.getByRole("region", { name: "About structure checks" })).toBeInTheDocument();
    expect(screen.getByRole("region", { name: "Results" })).toHaveTextContent(/Enter a SMILES/);
    expect(screen.getByRole("switch", { name: /Fix issues/i })).toHaveAttribute(
      "aria-checked",
      "true"
    );
  });

  it("reports the original and the standardized structure", async () => {
    renderView();
    fireEvent.change(screen.getByLabelText(/SMILES/), { target: { value: "C1=CC=CC=C1" } });
    fireEvent.click(screen.getByRole("button", { name: /Check structure/i }));
    const original = await screen.findByRole("region", { name: "Original structure" });
    expect(within(original).getByText("Valence problem on atom 2")).toBeInTheDocument();
    const fixed = screen.getByRole("region", { name: "Standardized structure" });
    expect(within(fixed).getAllByText("c1ccccc1").length).toBeGreaterThan(0);
    expect(within(fixed).getByText("No Errors Found")).toBeInTheDocument();
  });
});
