import { describe, it, expect } from "vitest";
import { render, screen } from "@testing-library/react";
import { ToolSkeleton } from "@/components/feedback/ToolSkeleton";
import { EmptyState } from "@/components/feedback/EmptyState";
import { RouteLoadingFallback } from "@/components/feedback/RouteLoadingFallback";

describe("ToolSkeleton", () => {
  it("renders descriptors variant", () => {
    render(<ToolSkeleton variant="descriptors" />);
    const container = screen.getByTestId("tool-skeleton");
    const skeletons = container.querySelectorAll('[data-testid="glass-skeleton"]');
    expect(skeletons.length).toBeGreaterThan(0);
  });

  it("renders molecule variant", () => {
    render(<ToolSkeleton variant="molecule" />);
    const container = screen.getByTestId("tool-skeleton");
    // molecule variant has a h-64 image placeholder
    const imageBlock = container.querySelector(".h-64");
    expect(imageBlock).toBeInTheDocument();
  });
});

describe("EmptyState", () => {
  it("renders default message", () => {
    render(<EmptyState />);
    expect(screen.getByText("Enter input above and submit to see results")).toBeInTheDocument();
  });
});

describe("RouteLoadingFallback", () => {
  it("renders skeleton blocks", () => {
    render(<RouteLoadingFallback />);
    const fallback = screen.getByTestId("route-loading-fallback");
    const skeletons = fallback.querySelectorAll('[data-testid="glass-skeleton"]');
    expect(skeletons.length).toBeGreaterThanOrEqual(3);
  });
});
