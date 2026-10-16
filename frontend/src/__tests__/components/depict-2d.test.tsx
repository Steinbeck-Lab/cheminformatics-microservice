import { describe, it, expect, vi, beforeEach } from "vitest";
import { act, render, screen, fireEvent, within } from "@testing-library/react";
import Depict2DView from "@/components/depict/Depict2DView";

vi.mock("@/services/convertService", () => ({
  generateInChI: vi.fn(async () => "InChI=1S/TEST"),
  generateInChIKey: vi.fn(async () => '"TESTKEY-UHFFFAOYSA-N"'), // the API quotes it
}));

vi.mock("@/services/chemService", () => ({
  calculateDescriptors: vi.fn(async () => ({
    molecular_formula: "C8H10N4O2",
    nplikeness: -1.09,
    molecular_weight: 194.19,
    alogp: -1.03,
    topological_polar_surface_area: 61.82,
    hydrogen_bond_acceptors: 3,
    hydrogen_bond_donors: 0,
    rotatable_bond_count: 0,
    aromatic_rings_count: 2,
    qed_drug_likeliness: 0.54,
    lipinski_rule_of_five_violations: 0,
  })),
  calculateNPLikeness: vi.fn(async () => 0.5),
}));

const CAFFEINE = "CN1C=NC2=C1C(=O)N(C(=O)N2C)C";
const mainImage = () => screen.getByRole("img", { name: /^2D depiction of/ });

describe("Depict2DView", () => {
  beforeEach(() => vi.clearAllMocks());

  it("opens with caffeine drawn by CDK", () => {
    render(<Depict2DView />);
    expect(screen.getByLabelText("SMILES")).toHaveValue(CAFFEINE);
    expect(mainImage()).toHaveAccessibleName("2D depiction of Caffeine rendered by CDK");
    expect(mainImage().getAttribute("src")).toContain("2D_enhanced");
    expect(mainImage().getAttribute("src")).toContain(encodeURIComponent("CN1C=NC2"));
  });

  it("switches the main depiction to RDKit", () => {
    render(<Depict2DView />);
    fireEvent.click(screen.getByRole("radio", { name: "RDKit" }));
    expect(mainImage().getAttribute("src")).toContain("toolkit=rdkit");
    expect(mainImage()).toHaveAccessibleName("2D depiction of Caffeine rendered by RDKit");
    expect(screen.getByRole("button", { name: /Use RDKit/ })).toHaveAttribute(
      "aria-pressed",
      "true"
    );
  });

  it("loads an example chip into the input and depicts it", () => {
    render(<Depict2DView />);
    fireEvent.click(screen.getByRole("button", { name: "Aspirin" }));
    expect(screen.getByLabelText("SMILES")).toHaveValue("CC(=O)OC1=CC=CC=C1C(=O)O");
    expect(screen.getByRole("button", { name: "Aspirin" })).toHaveAttribute("aria-pressed", "true");
    expect(mainImage()).toHaveAccessibleName(/Aspirin/);
  });

  it("clears, then asks for a SMILES when depicting nothing", () => {
    render(<Depict2DView />);
    fireEvent.click(screen.getByRole("button", { name: "Clear" }));
    expect(screen.getByLabelText("SMILES")).toHaveValue("");
    expect(screen.queryByRole("img", { name: /^2D depiction of/ })).toBeNull();
    fireEvent.click(screen.getByRole("button", { name: "Depict" }));
    expect(screen.getByText("Enter at least one SMILES string.")).toBeInTheDocument();
  });

  it("depicts several lines as a batch grid without one-molecule panels", () => {
    render(<Depict2DView />);
    fireEvent.change(screen.getByLabelText("SMILES"), {
      target: { value: `${CAFFEINE} Caffeine\nCC(=O)NC1=CC=C(C=C1)O Paracetamol` },
    });
    fireEvent.click(screen.getByRole("button", { name: "Depict" }));
    const all = screen.getByRole("region", { name: "All structures" });
    expect(within(all).getAllByRole("img", { name: /^2D depiction of / })).toHaveLength(2);
    expect(
      within(all).getByRole("slider", { name: "Rotation for Paracetamol" })
    ).toBeInTheDocument();
    expect(
      within(all).getByRole("button", { name: "Download Paracetamol as SVG" })
    ).toBeInTheDocument();
    for (const name of [
      "Depiction",
      "Compare toolkits",
      "Identifiers",
      "NP-likeness",
      "Descriptors",
    ])
      expect(screen.queryByRole("region", { name })).toBeNull();
  });

  it("loads the example batch and draws every structure", () => {
    render(<Depict2DView />);
    fireEvent.click(screen.getByRole("button", { name: "Load examples" }));
    const all = screen.getByRole("region", { name: "All structures" });
    expect(within(all).getAllByRole("img", { name: /^2D depiction of / })).toHaveLength(10);
    expect(within(all).getByRole("img", { name: /Dichlorobiphenyl/ })).toBeInTheDocument();
  });

  it("goes back to one-molecule panels when an example chip is picked", () => {
    render(<Depict2DView />);
    fireEvent.click(screen.getByRole("button", { name: "Load examples" }));
    fireEvent.click(screen.getByRole("button", { name: "Aspirin" }));
    expect(screen.queryByRole("region", { name: "All structures" })).toBeNull();
    expect(mainImage()).toHaveAccessibleName(/Aspirin/);
  });

  it("redraws the rotation once the slider settles, not on every step", () => {
    vi.useFakeTimers();
    try {
      render(<Depict2DView />);
      const slider = screen.getByRole("slider", { name: "Rotation" });
      fireEvent.change(slider, { target: { value: "45" } });
      fireEvent.change(slider, { target: { value: "90" } });
      expect(screen.getByText("90°")).toBeInTheDocument(); // the control is live
      expect(mainImage().getAttribute("src")).toContain("rotate=0&");
      act(() => vi.advanceTimersByTime(200));
      expect(mainImage().getAttribute("src")).toContain("rotate=90&");
    } finally {
      vi.useRealTimers();
    }
  });

  it("redraws at once when the rotation slider is released", () => {
    render(<Depict2DView />);
    const slider = screen.getByRole("slider", { name: "Rotation" });
    fireEvent.change(slider, { target: { value: "90" } });
    fireEvent.pointerUp(slider);
    expect(mainImage().getAttribute("src")).toContain("rotate=90&");
  });

  it("waits for typing to pause before drawing a highlight", () => {
    vi.useFakeTimers();
    try {
      render(<Depict2DView />);
      const input = screen.getByPlaceholderText("e.g. c1ccccc1");
      fireEvent.change(input, { target: { value: "c1" } });
      fireEvent.change(input, { target: { value: "c1ccccc1" } });
      expect(mainImage().getAttribute("src")).not.toContain("highlight=");
      act(() => vi.advanceTimersByTime(250));
      expect(mainImage().getAttribute("src")).toContain(
        `highlight=${encodeURIComponent("c1ccccc1")}`
      );
    } finally {
      vi.useRealTimers();
    }
  });

  it("shows identifiers, NP-likeness and descriptors for the structure", async () => {
    render(<Depict2DView />);
    expect(await screen.findByText("TESTKEY-UHFFFAOYSA-N")).toBeInTheDocument();
    expect(screen.getByText("InChI=1S/TEST")).toBeInTheDocument();
    const ids = screen.getByRole("region", { name: "Identifiers" });
    expect(within(ids).getByTestId("formula")).toHaveTextContent("C8H10N4O2");
    const np = screen.getByRole("region", { name: "NP-likeness" });
    expect(within(np).getByRole("img", { name: /NP-likeness score −1\.09/ })).toBeInTheDocument();
    expect(screen.getByText("194.19")).toBeInTheDocument();
  });

  it("marks CDK-only options as unavailable with RDKit", () => {
    render(<Depict2DView />);
    expect(screen.getByRole("radio", { name: "Minimal" })).toBeEnabled();
    fireEvent.click(screen.getByRole("radio", { name: "RDKit" }));
    expect(screen.getByRole("radio", { name: "Minimal" })).toBeDisabled();
  });
});

describe("Depict2DView deep link", () => {
  it("opens the structure from ?smiles=", () => {
    window.history.pushState({}, "", "/depict/2ddepiction?smiles=CCO");
    render(<Depict2DView />);
    expect(screen.getByLabelText("SMILES")).toHaveValue("CCO");
    expect(mainImage().getAttribute("src")).toContain("smiles=CCO");
    window.history.pushState({}, "", "/");
  });
});
