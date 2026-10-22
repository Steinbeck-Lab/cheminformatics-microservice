import { describe, it, expect, vi } from "vitest";
import { fireEvent, render, screen, within } from "@testing-library/react";
import { AppProvider } from "@/context/AppContext";
import AllFiltersView from "@/components/chem/AllFiltersView";
import api from "@/services/api";

vi.mock("@/services/api", () => ({
  default: {
    post: vi.fn(async () => ({
      data: ["CCO : T, T, T, F, T, T, 0.41, 1.98, -0.12"],
    })),
  },
}));

vi.mock("@/services/toolsService", () => ({
  applyChemicalFiltersDetailed: vi.fn(async () => ({
    total_molecules: 1,
    passing_molecules: 1,
    filter_operator: "OR",
    results: [
      {
        smiles: "CCO",
        valid: true,
        overall_pass: true,
        filters: {
          lipinski: { passes: true, details: "0 violations", violations: [] },
        },
      },
    ],
  })),
}));

const renderView = () =>
  render(
    <AppProvider>
      <AllFiltersView />
    </AppProvider>
  );

describe("AllFiltersView (bento)", () => {
  it("lays out input, about, options and results panels", () => {
    renderView();
    expect(screen.getByRole("region", { name: "Structure input" })).toBeInTheDocument();
    expect(screen.getByRole("region", { name: "About chemical filters" })).toBeInTheDocument();
    const options = screen.getByRole("region", { name: "Options" });
    expect(within(options).getByRole("switch", { name: /^PAINS/ })).toBeInTheDocument();
    expect(within(options).getByRole("switch", { name: /Rule of 3/ })).toBeInTheDocument();
    expect(within(options).getByLabelText(/QED/)).toBeInTheDocument();
    expect(
      within(options).getByRole("radiogroup", { name: "Filter match logic" })
    ).toBeInTheDocument();
    expect(screen.getByRole("region", { name: "Results" })).toHaveTextContent(/Enter SMILES/);
  });

  it("tabulates the filter results", async () => {
    renderView();
    fireEvent.change(screen.getByLabelText(/Molecule list/), { target: { value: "CCO" } });
    fireEvent.click(screen.getByRole("button", { name: /Apply filters/i }));
    const results = await screen.findByRole("region", { name: "Filter results" });
    expect(api.post).toHaveBeenCalled();
    expect(within(results).getByText("CCO")).toBeInTheDocument();
    expect(within(results).getByText("0.41")).toBeInTheDocument();
    expect(within(results).getAllByLabelText("Fail")).toHaveLength(1);
    expect(within(results).getByRole("button", { name: /Download/ })).toBeInTheDocument();
  });

  it("shows the detailed breakdown when detailed mode is on", async () => {
    renderView();
    fireEvent.click(screen.getByRole("switch", { name: /detailed violation/i }));
    fireEvent.change(screen.getByLabelText(/Molecule list/), { target: { value: "CCO" } });
    fireEvent.click(screen.getByRole("button", { name: /Apply filters/i }));
    const detail = await screen.findByRole("region", { name: "Detailed filter results" });
    expect(within(detail).getByText("0 violations")).toBeInTheDocument();
  });
});
