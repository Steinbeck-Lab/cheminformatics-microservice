import { describe, it, expect, vi } from "vitest";
import { fireEvent, render, screen, within } from "@testing-library/react";
import { AppProvider } from "@/context/AppContext";
import FixRadicalsView from "@/components/chem/FixRadicalsView";

vi.mock("@/services/chemService", () => ({
  fixRadicals: vi.fn(async () => ({
    fixed_smiles: "CC",
    radicals_detected: 1,
    radicals_fixed: 1,
  })),
}));

vi.mock("@/services/depictService", () => ({
  generate2DDepictionEnhanced: vi.fn(async () => "<svg></svg>"),
  generate2DDepiction: vi.fn(async () => "<svg></svg>"),
}));

const renderView = () =>
  render(
    <AppProvider>
      <FixRadicalsView />
    </AppProvider>
  );

describe("FixRadicalsView (bento)", () => {
  it("lays out input, about and results panels", () => {
    renderView();
    expect(screen.getByRole("region", { name: "Structure input" })).toBeInTheDocument();
    expect(screen.getByRole("region", { name: "About radical fixing" })).toBeInTheDocument();
    expect(screen.getByRole("region", { name: "Results" })).toHaveTextContent(/Enter a SMILES/);
  });

  it("shows the counts and the fixed structure", async () => {
    renderView();
    fireEvent.change(screen.getByLabelText(/SMILES/), { target: { value: "C[CH2]" } });
    fireEvent.click(screen.getByRole("button", { name: /Fix radicals/i }));
    const stats = await screen.findByRole("region", { name: "1 radical fixed" });
    expect(within(stats).getByText("Radicals detected")).toBeInTheDocument();
    expect(within(stats).getByText("100%")).toBeInTheDocument();
    const cmp = screen.getByRole("region", { name: "Structure comparison" });
    expect(await within(cmp).findByText("CC")).toBeInTheDocument();
    expect(within(cmp).getByRole("button", { name: /Copy fixed SMILES/i })).toBeInTheDocument();
  });
});
