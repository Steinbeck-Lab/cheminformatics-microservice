import { beforeEach, describe, it, expect, vi } from "vitest";
import { fireEvent, render, screen, within } from "@testing-library/react";
import { AppProvider } from "@/context/AppContext";
import StructureGenView from "@/components/tools/StructureGenView";
import toolsService from "@/services/toolsService";

vi.mock("@/services/toolsService", () => {
  const svc = { generateStructures: vi.fn() };
  return { default: svc, ...svc };
});

const OUTPUT = {
  output: {
    structures: ["CCO", "COC"],
    total_count: 2,
    generated_count: 2,
    limit_applied: false,
    settings: { "-P": "Require planarity" },
  },
};

const renderView = () =>
  render(
    <AppProvider>
      <StructureGenView />
    </AppProvider>
  );

const generate = (formula = "C2H6O") => {
  fireEvent.change(screen.getByLabelText(/Molecular formula/), { target: { value: formula } });
  fireEvent.click(screen.getByRole("button", { name: /Generate structures/i }));
};

describe("StructureGenView (bento)", () => {
  beforeEach(() => vi.clearAllMocks());

  it("lays out formula, about and results panels", () => {
    renderView();
    const input = screen.getByRole("region", { name: "Formula" });
    expect(within(input).getByLabelText(/Molecular formula/)).toBeInTheDocument();
    expect(within(input).getByRole("button", { name: /Generate structures/i })).toBeDisabled();
    expect(screen.getByRole("region", { name: "About structure generation" })).toHaveTextContent(
      /SURGE/
    );
    expect(screen.getByRole("region", { name: "Results" })).toHaveTextContent(/molecular formula/);
  });

  it("lists the generated isomers, with summary and selection", async () => {
    vi.mocked(toolsService.generateStructures).mockResolvedValueOnce(OUTPUT as unknown as string[]);
    renderView();
    generate();
    const results = await screen.findByRole("region", { name: "2 structures" });
    expect(toolsService.generateStructures).toHaveBeenCalledWith("C2H6O");
    expect(within(results).getByText("Require planarity")).toBeInTheDocument();
    expect(within(results).getByText("Complete")).toBeInTheDocument();
    const download = within(results).getByRole("button", { name: /Download selected \(0\)/ });
    expect(download).toBeDisabled();
    const first = within(results).getByRole("button", { name: "Select Isomer 1" });
    fireEvent.click(first);
    expect(first).toHaveAttribute("aria-pressed", "true");
    expect(within(results).getByRole("button", { name: /Download selected \(1\)/ })).toBeEnabled();
    fireEvent.click(within(results).getByRole("button", { name: /Select all/ }));
    expect(within(results).getByRole("button", { name: /Download selected \(2\)/ })).toBeEnabled();
    fireEvent.click(within(results).getByRole("button", { name: /Deselect all/ }));
    expect(first).toHaveAttribute("aria-pressed", "false");
  });

  it("notes when the result list was limited", async () => {
    vi.mocked(toolsService.generateStructures).mockResolvedValueOnce({
      output: { ...OUTPUT.output, total_count: 5400, generated_count: 2, limit_applied: true },
    } as unknown as string[]);
    renderView();
    generate("C7H8O");
    const results = await screen.findByRole("region", { name: "2 structures" });
    expect(within(results).getByText("Limited")).toBeInTheDocument();
    expect(within(results).getByText(/first 2 structures out of 5,400/)).toBeInTheDocument();
  });

  it("shows a message returned instead of structures as an error", async () => {
    vi.mocked(toolsService.generateStructures).mockResolvedValueOnce(
      "Too many heavy atoms" as unknown as string[]
    );
    renderView();
    generate("C20H42");
    const results = screen.getByRole("region", { name: "Results" });
    expect(await within(results).findByText("Too many heavy atoms")).toBeInTheDocument();
  });
});
