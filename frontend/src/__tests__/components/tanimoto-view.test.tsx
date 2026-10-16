import { describe, it, expect, vi } from "vitest";
import { fireEvent, render, screen, within } from "@testing-library/react";
import { AppProvider } from "@/context/AppContext";
import TanimotoView from "@/components/chem/TanimotoView";

const get = vi.fn(async (..._args: unknown[]) => ({ data: 0.75 }));
vi.mock("@/services/api", () => ({ default: { get: (...a: unknown[]) => get(...a) } }));

const renderView = () =>
  render(
    <AppProvider>
      <TanimotoView />
    </AppProvider>
  );

describe("TanimotoView (bento)", () => {
  it("lays out input, about and results panels", () => {
    renderView();
    expect(screen.getByRole("region", { name: "Structures to compare" })).toBeInTheDocument();
    expect(screen.getByRole("region", { name: "About Tanimoto similarity" })).toBeInTheDocument();
    expect(screen.getByRole("region", { name: "Results" })).toHaveTextContent(/Enter two SMILES/);
  });

  it("shows the similarity score in the results panel", async () => {
    renderView();
    fireEvent.change(screen.getByLabelText(/Molecule 1 SMILES/), { target: { value: "CCO" } });
    fireEvent.change(screen.getByLabelText(/Molecule 2 SMILES/), { target: { value: "CCN" } });
    fireEvent.click(screen.getByRole("button", { name: /Calculate similarity/i }));
    const results = await screen.findByRole("region", { name: "Tanimoto similarity" });
    expect(within(results).getByText("0.7500")).toBeInTheDocument();
    expect(within(results).getByText(/Good similarity/)).toBeInTheDocument();
    expect(get).toHaveBeenCalledWith("/chem/tanimoto", {
      params: {
        smiles: "CCO,CCN",
        toolkit: "rdkit",
        fingerprinter: "ECFP",
        nBits: 2048,
        radius: 2,
      },
    });
  });
});
