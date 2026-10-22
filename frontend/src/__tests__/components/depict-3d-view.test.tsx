import { afterEach, beforeEach, describe, it, expect, vi } from "vitest";
import { fireEvent, render, screen, waitFor, within } from "@testing-library/react";
import { AppProvider } from "@/context/AppContext";
import Depict3DView from "@/components/depict/Depict3DView";
import convertService from "@/services/convertService";

const MOLBLOCK = `
  RDKit          3D

  2  1  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.4000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
M  END
`;

vi.mock("@/services/convertService", () => {
  const generate3DCoordinates = vi.fn(async () => MOLBLOCK);
  return { default: { generate3DCoordinates }, generate3DCoordinates };
});

const makeViewer = () => ({
  addModel: vi.fn(() => ({ atoms: [] })),
  removeAllModels: vi.fn(),
  removeAllLabels: vi.fn(),
  clear: vi.fn(),
  resize: vi.fn(),
  render: vi.fn(),
  zoomTo: vi.fn(),
  setStyle: vi.fn(),
  addLabel: vi.fn(),
  spin: vi.fn(),
  setBackgroundColor: vi.fn(),
  setView: vi.fn(),
  getModel: vi.fn(() => ({})),
  pngURI: vi.fn(() => "data:image/png;base64,AAAA"),
});

const renderView = () =>
  render(
    <AppProvider>
      <Depict3DView />
    </AppProvider>
  );

const generate = (smiles = "CO") => {
  fireEvent.change(screen.getByLabelText(/SMILES/), { target: { value: smiles } });
  fireEvent.click(screen.getByRole("button", { name: /Generate 3D view/i }));
};

describe("Depict3DView (bento)", () => {
  beforeEach(() => vi.clearAllMocks());

  it("lays out input, options, viewer and results without WebGL", () => {
    renderView();
    const input = screen.getByRole("region", { name: "Structure input" });
    expect(within(input).getByLabelText(/SMILES/)).toBeInTheDocument();
    const options = screen.getByRole("region", { name: "Options" });
    expect(within(options).getByRole("radio", { name: "Open Babel" })).toHaveAttribute(
      "aria-checked",
      "true"
    );
    expect(within(options).getByRole("radio", { name: "Stick" })).toHaveAttribute(
      "aria-checked",
      "true"
    );
    expect(within(options).getByRole("combobox", { name: "Background" })).toHaveValue("white");
    expect(within(options).getByRole("combobox", { name: "Colour scheme" })).toHaveValue("default");
    expect(within(options).getByRole("switch", { name: /Atom labels/ })).toBeInTheDocument();
    expect(within(options).getByRole("switch", { name: /Spin/ })).toBeInTheDocument();
    const viewer = screen.getByRole("region", { name: "3D viewer" });
    expect(within(viewer).getByTestId("molecule-3d-container")).toBeInTheDocument();
    expect(screen.getByRole("region", { name: "Results" })).toHaveTextContent(/Enter a SMILES/);
  });

  it("shows the generated molblock verbatim even without 3Dmol", async () => {
    renderView();
    fireEvent.click(screen.getByRole("radio", { name: "RDKit" }));
    generate("CO");
    const results = await screen.findByRole("region", { name: "3D molblock" });
    expect(convertService.generate3DCoordinates).toHaveBeenCalledWith("CO", "rdkit");
    const pre = within(results).getByText((_, el) => el?.tagName === "PRE");
    expect(pre.textContent).toBe(MOLBLOCK);
    expect(within(results).getByRole("button", { name: /Copy molblock/i })).toBeInTheDocument();
    expect(within(results).getByRole("button", { name: /Download/i })).toBeInTheDocument();
    expect(within(results).getByText(/RDKit \(ETKDG\/MMFF94\)/)).toBeInTheDocument();
  });

  it("keeps an API error visible in the results panel", async () => {
    vi.mocked(convertService.generate3DCoordinates).mockRejectedValueOnce(new Error("boom"));
    renderView();
    generate("CO");
    const results = screen.getByRole("region", { name: "Results" });
    expect(
      await within(results).findByText(/Error generating 3D structure: boom/)
    ).toBeInTheDocument();
  });
});

describe("Depict3DView with a 3Dmol viewer", () => {
  let viewer: ReturnType<typeof makeViewer>;
  const sizes = ["clientWidth", "clientHeight"] as const;

  beforeEach(() => {
    vi.clearAllMocks();
    viewer = makeViewer();
    (window as unknown as { $3Dmol: unknown }).$3Dmol = {
      createViewer: vi.fn(() => viewer),
      elementColors: { rasmol: {} },
    };
    for (const s of sizes)
      Object.defineProperty(HTMLElement.prototype, s, { configurable: true, get: () => 500 });
  });

  afterEach(() => {
    delete (window as unknown as { $3Dmol?: unknown }).$3Dmol;
    for (const s of sizes) delete (HTMLElement.prototype as unknown as Record<string, unknown>)[s];
  });

  it("loads the molblock into the viewer and applies style, colours and spin", async () => {
    renderView();
    await waitFor(() =>
      expect(screen.queryByText(/Initializing 3D viewer/)).not.toBeInTheDocument()
    );
    generate("CO");
    await waitFor(() => expect(viewer.addModel).toHaveBeenCalledWith(MOLBLOCK, "mol"));
    expect(viewer.setStyle).toHaveBeenCalledWith(
      {},
      { stick: { radius: 0.15, colorscheme: undefined } }
    );

    fireEvent.click(screen.getByRole("radio", { name: "Sphere" }));
    await waitFor(() =>
      expect(viewer.setStyle).toHaveBeenCalledWith(
        {},
        { sphere: { scale: 0.3, colorscheme: undefined } }
      )
    );
    fireEvent.change(screen.getByRole("combobox", { name: "Colour scheme" }), {
      target: { value: "greenCarbon" },
    });
    await waitFor(() =>
      expect(viewer.setStyle).toHaveBeenCalledWith(
        {},
        { sphere: { scale: 0.3, colorscheme: "greenCarbon" } }
      )
    );
    fireEvent.change(screen.getByRole("combobox", { name: "Background" }), {
      target: { value: "black" },
    });
    await waitFor(() => expect(viewer.setBackgroundColor).toHaveBeenCalledWith("#000000"));
    fireEvent.click(screen.getByRole("switch", { name: /Spin/ }));
    await waitFor(() => expect(viewer.spin).toHaveBeenCalledWith(true));
  });

  it("takes a screenshot of the viewer", async () => {
    renderView();
    generate("CO");
    await waitFor(() => expect(viewer.addModel).toHaveBeenCalled());
    const cell = screen.getByRole("region", { name: "3D viewer" });
    const click = vi.spyOn(HTMLAnchorElement.prototype, "click").mockImplementation(() => {});
    fireEvent.click(within(cell).getByRole("button", { name: /Screenshot/ }));
    expect(viewer.pngURI).toHaveBeenCalled();
    expect(click).toHaveBeenCalled();
    click.mockRestore();
  });

  it("retries a viewer that failed to start", async () => {
    const lib = (window as unknown as { $3Dmol: { createViewer: ReturnType<typeof vi.fn> } })
      .$3Dmol;
    lib.createViewer.mockImplementationOnce(() => {
      throw new Error("no WebGL");
    });
    renderView();
    const results = screen.getByRole("region", { name: "Results" });
    expect(
      await within(results).findByText(/Failed to initialize 3D viewer: no WebGL/)
    ).toBeInTheDocument();
    fireEvent.click(within(results).getByRole("button", { name: /Try again/ }));
    await waitFor(() => expect(lib.createViewer).toHaveBeenCalledTimes(2));
    await waitFor(() =>
      expect(screen.queryByText(/Initializing 3D viewer/)).not.toBeInTheDocument()
    );
  });
});
