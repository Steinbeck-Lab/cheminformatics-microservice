import { beforeEach, describe, it, expect, vi } from "vitest";
import { fireEvent, render, screen, waitFor, within } from "@testing-library/react";
import { AppProvider } from "@/context/AppContext";
import StructureDrawView from "@/components/depict/StructureDrawView";

const ketcher = vi.hoisted(() => ({
  executeCommand: vi.fn(),
  resetEditor: vi.fn(),
  ready: true,
}));

vi.mock("@/hooks/useKetcherEditor", () => ({
  useKetcherEditor: () => ({
    ketcherFrame: { current: null },
    isEditorReady: ketcher.ready,
    executeCommand: ketcher.executeCommand,
    resetEditor: ketcher.resetEditor,
  }),
}));

const renderView = () =>
  render(
    <AppProvider>
      <StructureDrawView />
    </AppProvider>
  );

describe("StructureDrawView (bento)", () => {
  beforeEach(() => {
    vi.clearAllMocks();
    ketcher.ready = true;
    ketcher.executeCommand.mockImplementation(async (cmd: string) =>
      cmd === "getSmiles" ? "CCO" : undefined
    );
  });

  it("lays out input, results, about and editor panels", () => {
    renderView();
    const input = screen.getByRole("region", { name: "Structure input" });
    expect(within(input).getByLabelText("SMILES to edit")).toBeInTheDocument();
    expect(screen.getByRole("region", { name: "About the editor" })).toBeInTheDocument();
    expect(screen.getByRole("region", { name: "Results" })).toHaveTextContent(
      /generate its SMILES/
    );
    const editor = screen.getByRole("region", { name: "Structure editor" });
    expect(within(editor).getByTitle("Ketcher Editor")).toHaveAttribute(
      "src",
      "/standalone/index.html"
    );
    expect(within(editor).getByText("Editor ready")).toBeInTheDocument();
  });

  it("generates SMILES from the drawing and copies it", async () => {
    const writeText = vi.fn(async () => undefined);
    Object.assign(navigator, { clipboard: { writeText } });
    renderView();
    fireEvent.click(screen.getByRole("button", { name: "Generate SMILES" }));
    const results = await screen.findByRole("region", { name: "SMILES" });
    expect(ketcher.executeCommand).toHaveBeenCalledWith("getSmiles");
    expect(within(results).getByText("CCO")).toBeInTheDocument();
    fireEvent.click(within(results).getByRole("button", { name: /Copy SMILES/ }));
    await waitFor(() => expect(writeText).toHaveBeenCalledWith("CCO"));
  });

  it("loads an example SMILES into the editor", async () => {
    renderView();
    fireEvent.click(screen.getByRole("button", { name: "Aspirin" }));
    expect(screen.getByLabelText("SMILES to edit")).toHaveValue("CC(=O)OC1=CC=CC=C1C(=O)O");
    fireEvent.click(screen.getByRole("button", { name: "Load structure" }));
    await waitFor(() =>
      expect(ketcher.executeCommand).toHaveBeenCalledWith("setMolecule", [
        "CC(=O)OC1=CC=CC=C1C(=O)O",
      ])
    );
  });

  it("keeps an error visible when nothing is drawn", async () => {
    ketcher.executeCommand.mockImplementation(async () => "");
    renderView();
    fireEvent.click(screen.getByRole("button", { name: "Generate SMILES" }));
    const results = screen.getByRole("region", { name: "Results" });
    expect(await within(results).findByText(/No structure drawn/)).toBeInTheDocument();
  });

  it("offers a retry while the editor starts", () => {
    ketcher.ready = false;
    renderView();
    const editor = screen.getByRole("region", { name: "Structure editor" });
    expect(within(editor).getByText(/Initializing editor/)).toBeInTheDocument();
    fireEvent.click(within(editor).getByRole("button", { name: /Retry/ }));
    expect(ketcher.resetEditor).toHaveBeenCalled();
  });
});
