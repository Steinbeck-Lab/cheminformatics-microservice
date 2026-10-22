import { describe, it, expect, vi } from "vitest";
import { fireEvent, render, screen, within } from "@testing-library/react";
import { AppProvider } from "@/context/AppContext";
import CoconutPreProcessingView from "@/components/chem/CoconutPreProcessingView";
import { coconutPreprocessing } from "@/services/chemService";

vi.mock("@/services/chemService", () => ({
  coconutPreprocessing: vi.fn(async () => ({
    original: {
      representations: { "2D_MOL": "orig molblock", canonical_smiles: "OCC" },
      has_stereo: false,
      errors: [],
    },
    standardized: {
      representations: { "2D_MOL": "std molblock", canonical_smiles: "CCO" },
      has_stereo: true,
    },
    parent: { representations: { canonical_smiles: "CCO-parent" } },
  })),
}));

const renderView = () =>
  render(
    <AppProvider>
      <CoconutPreProcessingView />
    </AppProvider>
  );

describe("CoconutPreProcessingView (bento)", () => {
  it("lays out input, about and results panels", () => {
    renderView();
    const input = screen.getByRole("region", { name: "Structure input" });
    expect(within(input).getByRole("switch", { name: /3D coordinates/i })).toBeInTheDocument();
    expect(within(input).getByRole("switch", { name: /descriptors/i })).toBeInTheDocument();
    expect(
      screen.getByRole("region", { name: "About COCONUT pre-processing" })
    ).toBeInTheDocument();
    expect(screen.getByRole("region", { name: "Results" })).toHaveTextContent(/Enter a SMILES/);
  });

  it("shows each structure version in the results panel", async () => {
    renderView();
    fireEvent.click(screen.getByRole("switch", { name: /3D coordinates/i }));
    fireEvent.change(screen.getByLabelText(/SMILES/), { target: { value: "OCC" } });
    fireEvent.click(screen.getByRole("button", { name: /Process for COCONUT/i }));
    const results = await screen.findByRole("region", { name: "Original structure" });
    expect(coconutPreprocessing).toHaveBeenCalledWith("OCC", true, false);
    expect(within(results).getByText("orig molblock")).toBeInTheDocument();
    fireEvent.click(within(results).getByRole("radio", { name: "Parent" }));
    const parent = screen.getByRole("region", { name: "Parent structure" });
    expect(within(parent).getByText("CCO-parent")).toBeInTheDocument();
  });
});
