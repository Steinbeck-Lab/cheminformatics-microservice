import { describe, it, expect } from "vitest";
import { render, screen } from "@testing-library/react";
import { Button } from "@/components/ui/button";

describe("shadcn/ui Button", () => {
  it("renders with default variant", () => {
    render(<Button>Click me</Button>);
    expect(screen.getByRole("button", { name: "Click me" })).toBeInTheDocument();
  });

  it.todo("renders with destructive variant");
  it.todo("renders with size=icon for icon-only buttons");
  it.todo("forwards onClick handler");
});
