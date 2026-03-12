import { describe, it, expect } from "vitest";
import { render, screen } from "@testing-library/react";
import { Input } from "@/components/ui/input";

describe("shadcn/ui Input", () => {
  it("renders input with placeholder", () => {
    render(<Input placeholder="Enter text" />);
    expect(screen.getByPlaceholderText("Enter text")).toBeInTheDocument();
  });

  it.todo("handles value change");
});

describe("shadcn/ui Select", () => {
  it.todo("renders with options");
  it.todo("calls onValueChange on selection");
});

describe("shadcn/ui Textarea", () => {
  it.todo("renders textarea with rows prop");
});
