import React from "react";
import { ResizablePanelGroup, ResizablePanel, ResizableHandle } from "@/components/ui/resizable";
import { useMediaQuery } from "@/hooks/useMediaQuery";

interface ResizableToolPanelProps {
  /** Used for localStorage key via autoSaveId (e.g. "descriptors", "tanimoto"). */
  toolId: string;
  /** Left/top panel content (typically the input area). */
  inputContent: React.ReactNode;
  /** Right/bottom panel content (typically the output/result area). */
  outputContent: React.ReactNode;
  /** Default width percentage for the input (left) panel. */
  defaultInputSize?: number;
  /** Default width percentage for the output (right) panel. */
  defaultOutputSize?: number;
  /** Minimum panel size in percent. */
  minSize?: number;
}

/**
 * Responsive resizable split panel for expanded tool views.
 *
 * - Wide screens (>=1024px): horizontal split with draggable divider.
 *   Panel sizes persist to localStorage via react-resizable-panels autoSaveId.
 * - Narrow screens (<1024px): simple vertical stack, no resizing.
 */
export function ResizableToolPanel({
  toolId,
  inputContent,
  outputContent,
  defaultInputSize = 45,
  defaultOutputSize = 55,
  minSize = 25,
}: ResizableToolPanelProps) {
  const isWide = useMediaQuery("(min-width: 1024px)");

  if (!isWide) {
    // Mobile / narrow: stacked layout
    return (
      <div className="space-y-4">
        <div>{inputContent}</div>
        <div>{outputContent}</div>
      </div>
    );
  }

  // Wide: horizontal resizable split
  return (
    <ResizablePanelGroup
      direction="horizontal"
      autoSaveId={`panel-${toolId}`}
      className="min-h-[300px] rounded-lg"
    >
      <ResizablePanel defaultSize={defaultInputSize} minSize={minSize}>
        <div className="h-full overflow-auto p-4">{inputContent}</div>
      </ResizablePanel>

      <ResizableHandle withHandle />

      <ResizablePanel defaultSize={defaultOutputSize} minSize={minSize}>
        <div className="h-full overflow-auto p-4">{outputContent}</div>
      </ResizablePanel>
    </ResizablePanelGroup>
  );
}
