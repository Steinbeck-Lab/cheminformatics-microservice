import { describe, it, expect, vi } from "vitest";
import { fireEvent, render, screen, within } from "@testing-library/react";
import { AppProvider } from "@/context/AppContext";
import DescriptorsView from "@/components/chem/DescriptorsView";

vi.mock("@/services/chemService", () => ({
  getDescriptors: vi.fn(async () => ({ MolecularWeight: 180.16, TPSA: 63.6 })),
}));

const renderView = () =>
  render(
    <AppProvider>
      <DescriptorsView />
    </AppProvider>
  );

describe("DescriptorsView (bento)", () => {
  it("lays out input, options and results panels", () => {
    renderView();
    expect(screen.getByRole("region", { name: "Structure input" })).toBeInTheDocument();
    const options = screen.getByRole("region", { name: "Options" });
    expect(within(options).getByRole("radiogroup", { name: "Toolkit" })).toBeInTheDocument();
    expect(within(options).getByRole("radiogroup", { name: "Output format" })).toBeInTheDocument();
    expect(screen.getByRole("region", { name: "Results" })).toHaveTextContent(/Enter a SMILES/);
  });

  it("shows JSON output as copyable, pretty-printed JSON", async () => {
    renderView();
    fireEvent.change(screen.getByLabelText(/SMILES/), { target: { value: "CCO" } });
    fireEvent.click(screen.getByRole("button", { name: /Calculate descriptors/i }));
    const results = await screen.findByRole("region", { name: "2 descriptors" });
    expect(results.querySelector("pre")?.textContent).toBe(
      JSON.stringify({ MolecularWeight: 180.16, TPSA: 63.6 }, null, 2)
    );
    expect(within(results).queryByRole("table")).toBeNull();
    expect(within(results).getByRole("button", { name: /Copy JSON/ })).toBeInTheDocument();
  });

  it("keeps showing a result in the format it was fetched with after the toggle changes", async () => {
    renderView();
    fireEvent.change(screen.getByLabelText(/SMILES/), { target: { value: "CCO" } });
    fireEvent.click(screen.getByRole("button", { name: /Calculate descriptors/i }));
    const results = await screen.findByRole("region", { name: "2 descriptors" });
    fireEvent.click(screen.getByRole("radio", { name: "HTML" }));
    expect(results.querySelector("pre")).not.toBeNull();
    expect(within(results).queryByText(/Could not render descriptors/)).toBeNull();
    expect(within(results).getByText(/RDKIT · JSON/)).toBeInTheDocument();
  });
});
