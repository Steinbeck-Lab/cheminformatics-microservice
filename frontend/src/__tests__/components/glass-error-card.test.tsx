import { describe, it, expect, vi } from "vitest";
import { render, screen, fireEvent } from "@testing-library/react";
import { GlassErrorCard } from "@/components/feedback/GlassErrorCard";

describe("GlassErrorCard", () => {
  it("renders error message", () => {
    render(<GlassErrorCard message="Something went wrong" />);
    expect(screen.getByText("Something went wrong")).toBeInTheDocument();
  });

  it("renders retry button when onRetry provided", () => {
    render(<GlassErrorCard message="Error" onRetry={() => {}} />);
    expect(screen.getByRole("button", { name: /try again/i })).toBeInTheDocument();
  });

  it("does not render retry button when onRetry omitted", () => {
    render(<GlassErrorCard message="Error" />);
    expect(screen.queryByRole("button", { name: /try again/i })).not.toBeInTheDocument();
  });

  it("calls onRetry when retry button clicked", () => {
    const onRetry = vi.fn();
    render(<GlassErrorCard message="Error" onRetry={onRetry} />);
    fireEvent.click(screen.getByRole("button", { name: /try again/i }));
    expect(onRetry).toHaveBeenCalledOnce();
  });
});
