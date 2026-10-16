import { afterEach, beforeEach, describe, it, expect, vi } from "vitest";
import { fireEvent, render, screen, waitFor, within } from "@testing-library/react";
import userEvent from "@testing-library/user-event";
import { AppProvider } from "@/context/AppContext";
import OCRView from "@/components/ocsr/OCRView";
import OCSRPage from "@/pages/OCSRPage";
import ocsrService from "@/services/ocsrService";

vi.mock("@/services/ocsrService", () => {
  const processImage = vi.fn(async () => ({ smiles: ["CCO", "c1ccccc1"] }));
  const processImageUrl = vi.fn(async () => ({ smiles: "CN1C=NC2=C1C(=O)N(C(=O)N2C)C" }));
  return { default: { processImage, processImageUrl }, processImage, processImageUrl };
});

vi.mock("@/components/common/MoleculeCard", () => ({
  default: ({ smiles, title }: { smiles: string; title: string }) => (
    <figure data-testid="molecule">
      <figcaption>{title}</figcaption>
      <code>{smiles}</code>
    </figure>
  ),
}));

const renderView = () =>
  render(
    <AppProvider>
      <OCRView />
    </AppProvider>
  );

const urlApi = URL as unknown as Record<string, unknown>;

describe("OCRView (bento)", () => {
  beforeEach(() => {
    vi.clearAllMocks();
    urlApi.createObjectURL = vi.fn(() => "blob:preview");
    urlApi.revokeObjectURL = vi.fn();
  });
  afterEach(() => {
    delete urlApi.createObjectURL;
    delete urlApi.revokeObjectURL;
  });

  it("lays out input, options and results panels", () => {
    renderView();
    const input = screen.getByRole("region", { name: "Image input" });
    expect(within(input).getByRole("radio", { name: "Upload image" })).toHaveAttribute(
      "aria-checked",
      "true"
    );
    expect(within(input).getByText(/drag and drop/)).toBeInTheDocument();
    const options = screen.getByRole("region", { name: "Options" });
    expect(within(options).getByLabelText(/Reference name/)).toBeInTheDocument();
    expect(within(options).getByRole("switch", { name: /Hand-drawn model/ })).toHaveAttribute(
      "aria-checked",
      "false"
    );
    expect(screen.getByRole("region", { name: "Results" })).toHaveTextContent(/Upload an image/);
    expect(screen.getByRole("button", { name: /Extract structures/ })).toBeDisabled();
  });

  it("reads an uploaded image with the hand-drawn model", async () => {
    const user = userEvent.setup();
    const { container } = renderView();
    const file = new File(["png"], "structure.png", { type: "image/png" });
    await user.upload(container.querySelector<HTMLInputElement>('input[type="file"]')!, file);
    expect(await screen.findByText(/structure\.png/)).toBeInTheDocument();
    expect(screen.getByRole("img", { name: "Structure preview" })).toHaveAttribute(
      "src",
      "blob:preview"
    );
    fireEvent.change(screen.getByLabelText(/Reference name/), { target: { value: "Fig 1" } });
    fireEvent.click(screen.getByRole("switch", { name: /Hand-drawn model/ }));
    fireEvent.click(screen.getByRole("button", { name: /Extract structures/ }));
    const results = await screen.findByRole("region", { name: "Detected structures" });
    expect(ocsrService.processImage).toHaveBeenCalledWith(file, "Fig 1", true);
    expect(within(results).getAllByTestId("molecule")).toHaveLength(2);
    expect(within(results).getByText("c1ccccc1")).toBeInTheDocument();
    expect(within(results).getByText("2 structures")).toBeInTheDocument();
  });

  it("reads an image from a URL", async () => {
    renderView();
    fireEvent.click(screen.getByRole("radio", { name: "Image URL" }));
    fireEvent.change(screen.getByLabelText("Image URL"), {
      target: { value: " https://example.com/s.png " },
    });
    fireEvent.click(screen.getByRole("button", { name: /Extract structures/ }));
    const results = await screen.findByRole("region", { name: "Detected structures" });
    expect(ocsrService.processImageUrl).toHaveBeenCalledWith(
      "https://example.com/s.png",
      undefined,
      false
    );
    expect(within(results).getByText("Structure 1")).toBeInTheDocument();
  });

  it("keeps a processing error in the results panel", async () => {
    vi.mocked(ocsrService.processImageUrl).mockRejectedValueOnce(new Error("Server busy"));
    renderView();
    fireEvent.click(screen.getByRole("radio", { name: "Image URL" }));
    fireEvent.change(screen.getByLabelText("Image URL"), {
      target: { value: "https://example.com/s.png" },
    });
    fireEvent.click(screen.getByRole("button", { name: /Extract structures/ }));
    const results = screen.getByRole("region", { name: "Results" });
    expect(await within(results).findByText(/Server busy/)).toBeInTheDocument();
  });

  it("says so when no structure is found", async () => {
    vi.mocked(ocsrService.processImageUrl).mockResolvedValueOnce({ smiles: [] });
    renderView();
    fireEvent.click(screen.getByRole("radio", { name: "Image URL" }));
    fireEvent.change(screen.getByLabelText("Image URL"), {
      target: { value: "https://example.com/s.png" },
    });
    fireEvent.click(screen.getByRole("button", { name: /Extract structures/ }));
    await waitFor(() =>
      expect(screen.getByRole("region", { name: "Results" })).toHaveTextContent(
        /No chemical structures were detected/
      )
    );
  });
});

describe("OCSRPage", () => {
  it("has a title cell, the tool, and the DECIMER and MARCUS references", () => {
    render(
      <AppProvider>
        <OCSRPage />
      </AppProvider>
    );
    expect(
      screen.getByRole("heading", { level: 1, name: "Optical Chemical Structure Recognition" })
    ).toBeInTheDocument();
    expect(screen.getByRole("region", { name: "Image input" })).toBeInTheDocument();
    expect(screen.getByRole("link", { name: /Visit DECIMER/ })).toHaveAttribute(
      "href",
      "https://decimer.ai"
    );
    expect(screen.getByRole("link", { name: /Visit MARCUS/ })).toHaveAttribute(
      "href",
      "https://marcus.decimer.ai"
    );
    const refs = screen.getByRole("region", { name: "References" });
    expect(
      within(refs).getByRole("link", { name: "doi.org/10.1038/s41467-023-40782-0" })
    ).toBeInTheDocument();
    expect(
      within(refs).getByRole("link", { name: "doi.org/10.1039/d5dd00313j" })
    ).toBeInTheDocument();
  });
});
