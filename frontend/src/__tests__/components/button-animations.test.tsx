import { render, screen } from "@testing-library/react";
import { describe, expect, it } from "vitest";
import { Button } from "@/components/ui/button";

describe("Button animation classes", () => {
  it("renders with btn-hover-lift class", () => {
    render(<Button>Click me</Button>);
    const button = screen.getByRole("button", { name: /click me/i });
    expect(button.className).toContain("btn-hover-lift");
  });

  it("renders with focus-ring-animate class", () => {
    render(<Button>Focus me</Button>);
    const button = screen.getByRole("button", { name: /focus me/i });
    expect(button.className).toContain("focus-ring-animate");
  });

  it("ghost variant has scale-100 overrides that neutralize btn-hover-lift transforms", () => {
    render(<Button variant="ghost">Ghost</Button>);
    const button = screen.getByRole("button", { name: /ghost/i });
    // Ghost has hover:scale-100 and active:scale-100 to override btn-hover-lift
    expect(button.className).toContain("hover:scale-100");
    expect(button.className).toContain("active:scale-100");
  });

  it("link variant has scale-100 overrides that neutralize btn-hover-lift transforms", () => {
    render(<Button variant="link">Link</Button>);
    const button = screen.getByRole("button", { name: /link/i });
    expect(button.className).toContain("hover:scale-100");
    expect(button.className).toContain("active:scale-100");
  });

  it("ghost variant still includes focus-ring-animate for accessible focus rings", () => {
    render(<Button variant="ghost">Ghost Focus</Button>);
    const button = screen.getByRole("button", { name: /ghost focus/i });
    expect(button.className).toContain("focus-ring-animate");
  });

  it("default variant does not include scale-100 overrides", () => {
    render(<Button>Default</Button>);
    const button = screen.getByRole("button", { name: /default/i });
    expect(button.className).not.toContain("hover:scale-100");
    expect(button.className).not.toContain("active:scale-100");
  });
});
