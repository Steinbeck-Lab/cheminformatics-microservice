import { describe, it, expect } from "vitest";
import { render, screen } from "@testing-library/react";
import userEvent from "@testing-library/user-event";
import MoleculeDepiction2D from "@/components/depict/MoleculeDepiction2D";

const img = () => screen.getByRole("img", { name: "Caffeine" });
const param = (name: string) => new URL(img().getAttribute("src") ?? "").searchParams.get(name);

describe("MoleculeDepiction2D (porcelain)", () => {
  it("renders the depiction on a white stage", () => {
    render(<MoleculeDepiction2D smiles="CN1C=NC2=C1C(=O)N(C)C(=O)N2C" title="Caffeine" />);
    expect(img().closest(".pc-stage")).not.toBeNull();
    expect(param("smiles")).toBe("CN1C=NC2=C1C(=O)N(C)C(=O)N2C");
  });

  it("switches change the image URL", async () => {
    const user = userEvent.setup();
    render(<MoleculeDepiction2D smiles="CCO" title="Caffeine" />);
    expect(param("unicolor")).toBe("false");
    expect(param("CIP")).toBe("false");

    await user.click(screen.getByRole("switch", { name: /black & white/i }));
    expect(param("unicolor")).toBe("true");

    await user.click(screen.getByRole("switch", { name: /show r\/s, e\/z/i }));
    expect(param("CIP")).toBe("true");
  });

  it("H display select changes the URL for CDK", async () => {
    const user = userEvent.setup();
    render(<MoleculeDepiction2D smiles="CCO" title="Caffeine" toolkit="cdk" />);
    expect(param("hydrogen_display")).toBe("Smart");
    await user.selectOptions(screen.getByLabelText(/h display/i), "Explicit");
    expect(param("hydrogen_display")).toBe("Explicit");
  });

  it("reset restores the initial options", async () => {
    const user = userEvent.setup();
    render(<MoleculeDepiction2D smiles="CCO" title="Caffeine" showCIP />);
    await user.click(screen.getByRole("switch", { name: /black & white/i }));
    await user.click(screen.getByRole("switch", { name: /show r\/s, e\/z/i }));
    expect(param("CIP")).toBe("false");
    await user.click(screen.getByRole("button", { name: /reset/i }));
    expect(param("unicolor")).toBe("false");
    expect(param("CIP")).toBe("true");
  });

  it("two instances do not share element ids", () => {
    const { container } = render(
      <>
        <MoleculeDepiction2D smiles="CCO" title="A" />
        <MoleculeDepiction2D smiles="CCN" title="B" />
      </>
    );
    const ids = [...container.querySelectorAll("[id]")].map((e) => e.id);
    expect(new Set(ids).size).toBe(ids.length);
  });
});
