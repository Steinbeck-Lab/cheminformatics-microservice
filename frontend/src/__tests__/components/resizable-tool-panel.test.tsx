import React from "react";
import { describe, it, expect, vi, beforeEach } from "vitest";
import { render, screen } from "@testing-library/react";
import { ResizableToolPanel } from "@/components/layout/ResizableToolPanel";

// Mock useMediaQuery to control wide/narrow screen behavior
const mockUseMediaQuery = vi.fn();
vi.mock("@/hooks/useMediaQuery", () => ({
  useMediaQuery: (...args: unknown[]) => mockUseMediaQuery(...args),
}));

// Mock the resizable primitives so tests don't need the full library
vi.mock("@/components/ui/resizable", () => ({
  ResizablePanelGroup: ({
    children,
    autoSaveId,
    ...props
  }: React.PropsWithChildren<{ autoSaveId?: string }>) => (
    <div data-testid="resizable-panel-group" data-autosave-id={autoSaveId} {...props}>
      {children}
    </div>
  ),
  ResizablePanel: ({ children, ...props }: React.PropsWithChildren<Record<string, unknown>>) => (
    <div data-testid="resizable-panel" {...props}>
      {children}
    </div>
  ),
  ResizableHandle: ({
    withHandle,
    ...props
  }: { withHandle?: boolean } & Record<string, unknown>) => (
    <div data-testid="resizable-handle" data-with-handle={withHandle} {...props} />
  ),
}));

describe("ResizableToolPanel", () => {
  beforeEach(() => {
    mockUseMediaQuery.mockReset();
  });

  it("renders horizontal split layout on wide screens (>1024px)", () => {
    mockUseMediaQuery.mockReturnValue(true); // wide screen

    render(
      <ResizableToolPanel
        toolId="descriptors"
        inputContent={<div>Input</div>}
        outputContent={<div>Output</div>}
      />
    );

    expect(screen.getByTestId("resizable-panel-group")).toBeInTheDocument();
    expect(screen.getAllByTestId("resizable-panel")).toHaveLength(2);
    expect(screen.getByTestId("resizable-handle")).toBeInTheDocument();
  });

  it("renders stacked vertical layout on narrow screens (<1024px)", () => {
    mockUseMediaQuery.mockReturnValue(false); // narrow screen

    render(
      <ResizableToolPanel
        toolId="descriptors"
        inputContent={<div>Input</div>}
        outputContent={<div>Output</div>}
      />
    );

    // No resizable group on narrow screen
    expect(screen.queryByTestId("resizable-panel-group")).not.toBeInTheDocument();
    expect(screen.getByText("Input")).toBeInTheDocument();
    expect(screen.getByText("Output")).toBeInTheDocument();
  });

  it("renders inputContent in left/top panel", () => {
    mockUseMediaQuery.mockReturnValue(true);

    render(
      <ResizableToolPanel
        toolId="test-tool"
        inputContent={<div>My Input Area</div>}
        outputContent={<div>My Output Area</div>}
      />
    );

    const panels = screen.getAllByTestId("resizable-panel");
    expect(panels[0]).toHaveTextContent("My Input Area");
  });

  it("renders outputContent in right/bottom panel", () => {
    mockUseMediaQuery.mockReturnValue(true);

    render(
      <ResizableToolPanel
        toolId="test-tool"
        inputContent={<div>My Input Area</div>}
        outputContent={<div>My Output Area</div>}
      />
    );

    const panels = screen.getAllByTestId("resizable-panel");
    expect(panels[1]).toHaveTextContent("My Output Area");
  });

  it("uses autoSaveId for localStorage persistence", () => {
    mockUseMediaQuery.mockReturnValue(true);

    render(
      <ResizableToolPanel
        toolId="tanimoto"
        inputContent={<div>Input</div>}
        outputContent={<div>Output</div>}
      />
    );

    const group = screen.getByTestId("resizable-panel-group");
    expect(group).toHaveAttribute("data-autosave-id", "panel-tanimoto");
  });
});
