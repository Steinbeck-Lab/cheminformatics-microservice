import { beforeEach, describe, it, expect, vi } from "vitest";
import { fireEvent, render, screen, waitFor, within } from "@testing-library/react";
import { AppProvider } from "@/context/AppContext";
import RInChIView from "@/components/tools/RInChIView";
import { convertRinchiToFileText, convertRxnfileToRinchi } from "@/utils/rinchiUtils";

const RXN = `$RXN

  Ketcher 10162612122D 1   1.00000     0.00000     0

  1  1
$MOL

  Ketcher 10162612122D 1   1.00000     0.00000     0

  1  0  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
M  END
$MOL

  Ketcher 10162612122D 1   1.00000     0.00000     0

  1  0  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
M  END
`;

const RINCHI = "RInChI=1.00.1S/CH4/h1H4<>H2O/h1H2/d+";

const ketcher = vi.hoisted(() => ({
  executeCommand: vi.fn(),
  resetEditor: vi.fn(),
}));

vi.mock("@/hooks/useKetcherEditor", () => ({
  useKetcherEditor: () => ({
    ketcherFrame: { current: null },
    isEditorReady: true,
    executeCommand: ketcher.executeCommand,
    resetEditor: ketcher.resetEditor,
  }),
}));

vi.mock("@/utils/rinchiUtils", () => ({
  VERSION_OPTIONS: [["1.1", "RInChI 1.1-dev with Latest InChI"]],
  DEFAULT_VERSION: "1.1",
  loadRinchiModule: vi.fn(async () => undefined),
  convertRxnfileToRinchi: vi.fn(async () => ({
    return_code: 0,
    rinchi: "RInChI=1.00.1S/CH4/h1H4<>H2O/h1H2/d+",
    rauxinfo: "RAuxInfo=1.00.1/0/N:1/rA:1C/rB:/rC:0,0,0;<>0/N:1/rA:1O/rB:/rC:0,0,0;",
  })),
  generateRinchiKey: vi.fn(async (_r: string, type: string) => ({
    return_code: 0,
    rinchikey: `${type}-RInChIKey=SA-FUHFF`,
  })),
  convertRinchiToFileText: vi.fn(async () => ({ return_code: 0, fileText: RXN })),
}));

const renderView = () =>
  render(
    <AppProvider>
      <RInChIView />
    </AppProvider>
  );

const generateButton = () => screen.getByRole("button", { name: /Generate RInChI/i });

describe("RInChIView (bento)", () => {
  beforeEach(() => {
    vi.clearAllMocks();
    ketcher.executeCommand.mockImplementation(async (cmd: string) =>
      cmd === "getRxn" ? RXN : undefined
    );
  });

  it("lays out input, options, editor and results panels", async () => {
    renderView();
    const input = screen.getByRole("region", { name: "Reaction input" });
    expect(within(input).getByRole("radiogroup", { name: "Input method" })).toBeInTheDocument();
    const options = screen.getByRole("region", { name: "Options" });
    expect(within(options).getByRole("combobox", { name: "RInChI version" })).toHaveValue("1.1");
    expect(within(options).getByRole("switch", { name: /Force equilibrium/ })).toBeInTheDocument();
    expect(within(options).getByRole("radiogroup", { name: "File format" })).toBeInTheDocument();
    const editor = screen.getByRole("region", { name: "Reaction editor" });
    expect(within(editor).getByTitle("Ketcher Editor")).toBeInTheDocument();
    expect(screen.getByRole("region", { name: "Results" })).toHaveTextContent(
      /generate its RInChI/
    );
    await waitFor(() => expect(generateButton()).toBeEnabled());
  });

  it("generates RInChI and all three keys from the drawn reaction", async () => {
    renderView();
    await waitFor(() => expect(generateButton()).toBeEnabled());
    fireEvent.click(screen.getByRole("switch", { name: /Force equilibrium/ }));
    fireEvent.click(generateButton());
    const results = await screen.findByRole("region", { name: "RInChI" });
    const [sent, equilibrium, version] = vi.mocked(convertRxnfileToRinchi).mock.calls[0]!;
    expect(sent).toBe(RXN);
    expect(equilibrium).toBe(true);
    expect(version).toBe("1.1");
    expect(within(results).getByText(RINCHI)).toBeInTheDocument();
    expect(await within(results).findByText("Long-RInChIKey=SA-FUHFF")).toBeInTheDocument();
    expect(within(results).getByText("Short-RInChIKey=SA-FUHFF")).toBeInTheDocument();
    expect(within(results).getByText("Web-RInChIKey=SA-FUHFF")).toBeInTheDocument();
  });

  it("converts a RInChI to an RD file without the editor, keeping column spacing", async () => {
    renderView();
    fireEvent.click(screen.getByRole("radio", { name: "RInChI → file" }));
    expect(screen.queryByRole("region", { name: "Reaction editor" })).not.toBeInTheDocument();
    fireEvent.click(screen.getByRole("radio", { name: "RD" }));
    fireEvent.change(screen.getByLabelText(/^RInChI to convert/), { target: { value: RINCHI } });
    const convert = await screen.findByRole("button", { name: /Convert to RD file/i });
    await waitFor(() => expect(convert).toBeEnabled());
    fireEvent.click(convert);
    const results = await screen.findByRole("region", { name: "RD file" });
    expect(convertRinchiToFileText).toHaveBeenCalledWith(RINCHI, "", "RD", "1.1");
    const pre = within(results).getByTestId("rinchi-file-output");
    expect(pre.tagName).toBe("PRE");
    expect(pre.textContent).toContain("M  END");
    expect(within(results).getByRole("button", { name: "Copy RD file" })).toBeInTheDocument();
  });

  it("keeps an error visible when no reaction is drawn", async () => {
    ketcher.executeCommand.mockImplementation(async () => "");
    renderView();
    await waitFor(() => expect(generateButton()).toBeEnabled());
    fireEvent.click(generateButton());
    const results = screen.getByRole("region", { name: "Results" });
    expect(await within(results).findByRole("alert")).toHaveTextContent(/No reaction found/);
  });
});
