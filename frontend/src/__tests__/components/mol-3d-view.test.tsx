import { describe, it, expect, vi, beforeEach } from "vitest";
import { fireEvent, render, screen, within } from "@testing-library/react";
import { AppProvider } from "@/context/AppContext";
import Mol3DView from "@/components/convert/Mol3DView";
import { generate3DCoordinates } from "@/services/convertService";

const MOLBLOCK = `
  CDK     08302311363D

  3  2  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2990    0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.5981    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  1  0  0  0  0
M  END
`;

vi.mock("@/services/convertService", () => ({
  generate3DCoordinates: vi.fn(async () => MOLBLOCK),
}));

const renderView = () =>
  render(
    <AppProvider>
      <Mol3DView />
    </AppProvider>
  );

describe("Mol3DView (bento)", () => {
  beforeEach(() => vi.clearAllMocks());

  it("lays out input, options and results panels", () => {
    renderView();
    const input = screen.getByRole("region", { name: "Structure input" });
    expect(within(input).getByLabelText(/SMILES/)).toBeInTheDocument();
    const options = screen.getByRole("region", { name: "Options" });
    expect(within(options).getByRole("radiogroup", { name: "Toolkit" })).toBeInTheDocument();
    expect(within(options).getByRole("radio", { name: "Open Babel" })).toHaveAttribute(
      "aria-checked",
      "true"
    );
    expect(screen.getByRole("region", { name: "Results" })).toHaveTextContent(/Enter a SMILES/);
  });

  it("shows the 3D molfile verbatim with copy and download", async () => {
    renderView();
    fireEvent.click(screen.getByRole("radio", { name: "RDKit" }));
    fireEvent.change(screen.getByLabelText(/SMILES/), { target: { value: "CCO" } });
    fireEvent.click(screen.getByRole("button", { name: /Generate 3D coordinates/i }));
    const results = await screen.findByRole("region", { name: "3D molfile" });
    expect(generate3DCoordinates).toHaveBeenCalledWith("CCO", "rdkit");
    const pre = within(results).getByText((_, el) => el?.tagName === "PRE");
    expect(pre.textContent).toBe(MOLBLOCK);
    expect(within(results).getByRole("button", { name: /Copy molblock/i })).toBeInTheDocument();
    expect(within(results).getByRole("button", { name: /Download/i })).toBeInTheDocument();
  });
});
