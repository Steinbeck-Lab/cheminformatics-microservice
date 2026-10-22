import { describe, it, expect, vi } from "vitest";
import { render, screen } from "@testing-library/react";
import Molecule3D from "@/components/common/Molecule3D";

vi.mock("@/services/convertService", () => ({
  generate3DCoordinates: vi.fn(async () => "not a molblock"),
}));

describe("Molecule3D", () => {
  it("is a labelled image and survives a browser without WebGL", () => {
    render(<Molecule3D smiles="CCO" label="Ethanol in 3D" />);
    expect(screen.getByRole("img", { name: "Ethanol in 3D" })).toBeInTheDocument();
  });
});
