import { beforeEach, describe, it, expect, vi } from "vitest";
import { fireEvent, render, screen, within } from "@testing-library/react";
import { AppProvider } from "@/context/AppContext";
import SugarRemovalView from "@/components/tools/SugarRemovalView";
import {
  extractAglyconeAndSugars,
  getAglyconeAndSugarIndices,
  getSugarInfo,
  removeAllSugars,
  removeCircularSugars,
} from "@/services/toolsService";

vi.mock("@/services/toolsService", () => ({
  getSugarInfo: vi.fn(async () => "The molecule contains circular and linear sugars"),
  removeLinearSugars: vi.fn(async () => "LINEAR_FREE"),
  removeCircularSugars: vi.fn(async () => "CIRCULAR_FREE"),
  removeAllSugars: vi.fn(async () => "CCCC"),
  extractAglyconeAndSugars: vi.fn(async () => ["AGLY", "OCC1OC(O)C(O)C(O)C1O", "OCC(O)CO"]),
  getAglyconeAndSugarIndices: vi.fn(async () => [
    [0, 1],
    [2, 3, 4],
  ]),
}));

const GLYCOSIDE = "OCC1OC(OC2=CC=CC=C2)C(O)C(O)C1O";

const renderView = () =>
  render(
    <AppProvider>
      <SugarRemovalView />
    </AppProvider>
  );

const enterSmiles = () =>
  fireEvent.change(screen.getByLabelText(/SMILES String/), { target: { value: GLYCOSIDE } });

describe("SugarRemovalView (bento)", () => {
  beforeEach(() => vi.clearAllMocks());

  it("lays out input, about, options and results panels", () => {
    renderView();
    const input = screen.getByRole("region", { name: "Structure input" });
    expect(within(input).getByRole("radiogroup", { name: "Operation" })).toBeInTheDocument();
    expect(within(input).getByRole("radio", { name: "Detect" })).toHaveAttribute(
      "aria-checked",
      "true"
    );
    expect(within(input).getByRole("button", { name: /Detect sugars/i })).toBeDisabled();
    const about = screen.getByRole("region", { name: "About sugar removal" });
    expect(within(about).getByRole("link", { name: /MORTAR/ })).toBeInTheDocument();
    const options = screen.getByRole("region", { name: "Options" });
    // every option is visible without expanding anything
    expect(
      within(options).getByRole("group", { name: "Circular sugar detection options" })
    ).toBeInTheDocument();
    expect(
      within(options).getByRole("switch", { name: /Require O-glycosidic bonds/ })
    ).toBeInTheDocument();
    expect(screen.getByRole("region", { name: "Results" })).toHaveTextContent(/SMILES string/);
  });

  it("detects sugars and highlights them by type", async () => {
    renderView();
    enterSmiles();
    fireEvent.click(screen.getByRole("button", { name: /Detect sugars/i }));
    const results = await screen.findByRole("region", { name: "Sugar detection" });
    expect(getSugarInfo).toHaveBeenCalledWith(
      GLYCOSIDE,
      expect.objectContaining({ gly_bond: false })
    );
    expect(getAglyconeAndSugarIndices).toHaveBeenCalledTimes(2);
    expect(within(results).getByText(/contains circular and linear sugars/)).toBeInTheDocument();
    const highlight = within(results).getByRole("radiogroup", { name: "Highlight" });
    fireEvent.click(within(highlight).getByRole("radio", { name: "Both" }));
    expect(
      within(results).getByText("Both circular and linear sugars highlighted")
    ).toBeInTheDocument();
    expect(within(results).getByAltText("Molecule with highlighted sugars")).toBeInTheDocument();
  });

  it("removes both sugar types and shows before and after", async () => {
    renderView();
    enterSmiles();
    fireEvent.click(screen.getByRole("radio", { name: "Remove" }));
    fireEvent.click(screen.getByRole("radio", { name: "Both" }));
    // basic options are open in remove mode
    const options = screen.getByRole("region", { name: "Options" });
    expect(within(options).getByRole("combobox", { name: /Preservation mode/ })).toHaveValue("2");
    fireEvent.click(screen.getByRole("button", { name: /Remove sugars/i }));
    const results = await screen.findByRole("region", { name: "Sugar removal" });
    expect(removeAllSugars).toHaveBeenCalledWith(
      GLYCOSIDE,
      expect.objectContaining({ only_terminal: true, preservation_mode: 2 })
    );
    expect(removeCircularSugars).not.toHaveBeenCalled();
    expect(within(results).getByRole("figure", { name: "Before Removal" })).toBeInTheDocument();
    expect(within(results).getByRole("figure", { name: "After Removal" })).toBeInTheDocument();
  });

  it("extracts the aglycone and each sugar into a grid", async () => {
    renderView();
    enterSmiles();
    fireEvent.click(screen.getByRole("radio", { name: "Extract" }));
    fireEvent.click(screen.getByRole("button", { name: /Extract aglycone and sugars/i }));
    const results = await screen.findByRole("region", { name: "Aglycone and sugars" });
    expect(extractAglyconeAndSugars).toHaveBeenCalledWith(
      GLYCOSIDE,
      expect.objectContaining({ extract_circular_sugars: true, extract_linear_sugars: false })
    );
    expect(within(results).getByRole("figure", { name: "Aglycone" })).toBeInTheDocument();
    expect(within(results).getByRole("figure", { name: "Sugar 1" })).toBeInTheDocument();
    expect(within(results).getByRole("figure", { name: "Sugar 2" })).toBeInTheDocument();
  });
});
