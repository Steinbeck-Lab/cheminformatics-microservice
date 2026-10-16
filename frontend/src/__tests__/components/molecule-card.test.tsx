import { describe, it, expect, vi, beforeEach } from "vitest";
import { render, screen } from "@testing-library/react";
import userEvent from "@testing-library/user-event";
import MoleculeCard from "@/components/common/MoleculeCard";
import { AppProvider } from "@/context/AppContext";

const renderCard = (props: Partial<Parameters<typeof MoleculeCard>[0]> = {}) =>
  render(
    <AppProvider>
      <MoleculeCard smiles="CCO" title="Ethanol" {...props} />
    </AppProvider>
  );

describe("MoleculeCard (porcelain)", () => {
  let writeText: ReturnType<typeof vi.fn>;
  beforeEach(() => {
    writeText = vi.fn(async () => {});
    Object.defineProperty(navigator, "clipboard", {
      configurable: true,
      value: { writeText },
    });
  });

  it("copy button writes the SMILES to the clipboard", async () => {
    renderCard();
    // userEvent installs its own clipboard stub, so click without setup()
    screen.getByRole("button", { name: "Copy SMILES" }).click();
    expect(writeText).toHaveBeenCalledWith("CCO");
    expect(await screen.findByRole("button", { name: "SMILES copied" })).toBeInTheDocument();
  });

  it("shows the depiction on a white stage and the SMILES in mono", () => {
    renderCard();
    const img = screen.getByRole("img", { name: /chemical structure of ethanol/i });
    expect(img.closest(".pc-stage")).not.toBeNull();
    expect(screen.getByText("CCO").closest(".pc-mono")).not.toBeNull();
  });

  it("labels every action", () => {
    renderCard();
    expect(screen.getByRole("button", { name: "Download structure as SDF" })).toBeInTheDocument();
    expect(screen.getByRole("button", { name: "View molecule details" })).toBeInTheDocument();
  });

  it("card click calls onClick with the SMILES", async () => {
    const onClick = vi.fn();
    const user = userEvent.setup();
    renderCard({ onClick });
    await user.click(screen.getByRole("button", { name: "Ethanol" }));
    expect(onClick).toHaveBeenCalledWith("CCO");
  });
});
