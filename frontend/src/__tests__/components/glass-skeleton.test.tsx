import { describe, it, expect } from "vitest";
import { render, screen } from "@testing-library/react";
import { GlassSkeleton } from "@/components/feedback/GlassSkeleton";

describe("GlassSkeleton", () => {
  it("carries the silk sheen (its animation is gated on reduced motion in CSS)", () => {
    render(<GlassSkeleton />);
    const skeleton = screen.getByTestId("glass-skeleton");
    expect(skeleton.className).toContain("pc-sheen");
    expect(skeleton.className).not.toContain("animate-pulse");
  });

  it("renders text variant by default", () => {
    render(<GlassSkeleton />);
    const skeleton = screen.getByTestId("glass-skeleton");
    expect(skeleton.className).toContain("h-4");
  });

  it("renders image variant", () => {
    render(<GlassSkeleton variant="image" />);
    const skeleton = screen.getByTestId("glass-skeleton");
    expect(skeleton.className).toContain("h-48");
  });

  it("accepts custom className", () => {
    render(<GlassSkeleton className="custom-test-class" />);
    const skeleton = screen.getByTestId("glass-skeleton");
    expect(skeleton.className).toContain("custom-test-class");
  });
});
