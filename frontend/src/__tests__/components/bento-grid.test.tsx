import React from "react";
import { describe, it, expect, vi } from "vitest";
import { render, screen } from "@testing-library/react";
import { BentoGrid, BentoCard } from "@/components/common/BentoGrid";
import { Beaker } from "lucide-react";

// Mock motion/react for simpler testing
vi.mock("motion/react", () => ({
  motion: {
    div: ({
      children,
      className,
      animate,
      ...props
    }: React.PropsWithChildren<{
      className?: string;
      animate?: { opacity?: number; scale?: number };
    }>) => (
      <div
        className={className}
        data-opacity={animate?.opacity}
        data-scale={animate?.scale}
        {...props}
      >
        {children}
      </div>
    ),
  },
  AnimatePresence: ({ children }: React.PropsWithChildren) => <>{children}</>,
  LayoutGroup: ({ children }: React.PropsWithChildren) => (
    <div data-testid="layout-group">{children}</div>
  ),
}));

// Mock ResizableToolPanel
vi.mock("@/components/layout/ResizableToolPanel", () => ({
  ResizableToolPanel: ({
    toolId,
    inputContent,
    outputContent,
  }: {
    toolId: string;
    inputContent: React.ReactNode;
    outputContent: React.ReactNode;
  }) => (
    <div data-testid="resizable-tool-panel" data-tool-id={toolId}>
      <div data-testid="panel-input">{inputContent}</div>
      <div data-testid="panel-output">{outputContent}</div>
    </div>
  ),
}));

describe("BentoGrid", () => {
  it("renders a CSS grid container with correct column classes", () => {
    const { container } = render(
      <BentoGrid pageKey="chem" expandedId={null}>
        <div>Card 1</div>
      </BentoGrid>
    );
    const grid = container.querySelector(".grid");
    expect(grid?.className).toContain("sm:grid-cols-2");
    expect(grid?.className).toContain("lg:grid-cols-3");
    expect(grid?.className).toContain("xl:grid-cols-4");
  });

  it("switches to single column when expandedId is set", () => {
    const { container } = render(
      <BentoGrid pageKey="chem" expandedId="descriptors">
        <div>Card 1</div>
      </BentoGrid>
    );
    const grid = container.querySelector(".grid");
    expect(grid?.className).toContain("grid-cols-1");
    expect(grid?.className).not.toContain("sm:grid-cols-2");
  });

  it("wraps children in LayoutGroup", () => {
    render(
      <BentoGrid pageKey="chem" expandedId={null}>
        <div>Child</div>
      </BentoGrid>
    );
    expect(screen.getByTestId("layout-group")).toBeInTheDocument();
  });
});

describe("BentoCard", () => {
  const baseProps = {
    id: "descriptors",
    pageKey: "chem",
    name: "Descriptors",
    description: "Calculate molecular descriptors",
    icon: Beaker,
    size: "compact" as const,
    isExpanded: false,
    isOtherExpanded: false,
    onToggle: vi.fn(),
  };

  it("renders card with name, description, and icon", () => {
    render(<BentoCard {...baseProps} />);
    expect(screen.getByText("Descriptors")).toBeInTheDocument();
    expect(screen.getByText("Calculate molecular descriptors")).toBeInTheDocument();
  });

  it("calls onToggle when clicked", async () => {
    const onToggle = vi.fn();
    const { container } = render(<BentoCard {...baseProps} onToggle={onToggle} />);
    container.firstElementChild?.dispatchEvent(new MouseEvent("click", { bubbles: true }));
    expect(onToggle).toHaveBeenCalled();
  });

  it("renders children when isExpanded is true", () => {
    render(
      <BentoCard {...baseProps} isExpanded={true}>
        <div>Tool content</div>
      </BentoCard>
    );
    expect(screen.getByText("Tool content")).toBeInTheDocument();
  });

  it("shows close button when expanded", () => {
    render(
      <BentoCard {...baseProps} isExpanded={true}>
        <div>Content</div>
      </BentoCard>
    );
    expect(screen.getByLabelText("Close expanded tool")).toBeInTheDocument();
  });

  it("renders ResizableToolPanel when inputContent and outputContent are provided", () => {
    render(
      <BentoCard
        {...baseProps}
        isExpanded={true}
        inputContent={<div>Input area</div>}
        outputContent={<div>Output area</div>}
      />
    );
    const panel = screen.getByTestId("resizable-tool-panel");
    expect(panel).toBeInTheDocument();
    expect(panel).toHaveAttribute("data-tool-id", "descriptors");
    expect(screen.getByText("Input area")).toBeInTheDocument();
    expect(screen.getByText("Output area")).toBeInTheDocument();
  });

  it("renders children (not ResizableToolPanel) when only children are provided", () => {
    render(
      <BentoCard {...baseProps} isExpanded={true}>
        <div>Legacy children content</div>
      </BentoCard>
    );
    expect(screen.getByText("Legacy children content")).toBeInTheDocument();
    expect(screen.queryByTestId("resizable-tool-panel")).not.toBeInTheDocument();
  });
});
