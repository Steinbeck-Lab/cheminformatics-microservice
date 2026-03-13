import { describe, it, expect } from "vitest";
import { render, screen } from "@testing-library/react";
import { GlassSkeleton } from "@/components/feedback/GlassSkeleton";

describe("GlassSkeleton", () => {
  it("renders with shimmer animation class", () => {
    render(<GlassSkeleton />);
    const skeleton = screen.getByTestId("glass-skeleton");
    const shimmer = skeleton.querySelector(".animate-shimmer");
    expect(shimmer).toBeInTheDocument();
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
