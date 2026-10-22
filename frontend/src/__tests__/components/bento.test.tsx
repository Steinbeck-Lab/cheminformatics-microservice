import { describe, it, expect } from "vitest";
import { render, screen } from "@testing-library/react";
import { Bento, BentoCell, ToolPageHead } from "@/components/common/Bento";

describe("Bento template", () => {
  it("each cell is a region named by its label", () => {
    render(
      <Bento>
        <BentoCell title="Structure input">body</BentoCell>
        <BentoCell title="Depiction options">opts</BentoCell>
      </Bento>
    );
    expect(screen.getByRole("region", { name: "Structure input" })).toHaveTextContent("body");
    expect(
      screen.getByRole("heading", { level: 2, name: "Depiction options" })
    ).toBeInTheDocument();
  });

  it("cells are glass, full width by default, and take page placement classes", () => {
    render(<BentoCell title="Depiction" className="lg:col-[1/8]" />);
    const cell = screen.getByRole("region", { name: "Depiction" });
    expect(cell).toHaveClass("glass", "col-span-full", "lg:col-[1/8]");
  });

  it("shows meta and the madder accent", () => {
    render(<BentoCell title="NP-likeness" meta="-5 to +5" accent="madder" />);
    expect(screen.getByText("-5 to +5")).toBeInTheDocument();
    expect(screen.getByRole("heading", { name: "NP-likeness" })).toHaveAttribute(
      "data-accent",
      "madder"
    );
  });

  it("page head renders the h1, description and status", () => {
    render(<ToolPageHead title="2D Depiction" description="CDK or RDKit" status="Ready" />);
    expect(screen.getByRole("heading", { level: 1, name: "2D Depiction" })).toBeInTheDocument();
    expect(screen.getByText("CDK or RDKit")).toBeInTheDocument();
    expect(screen.getByText("Ready")).toBeInTheDocument();
  });
});
