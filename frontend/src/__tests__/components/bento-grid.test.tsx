import { describe, it, expect } from "vitest";

describe("BentoGrid", () => {
  it.todo("renders a CSS grid container with correct column classes");
  it.todo("switches to single column when expandedId is set");
  it.todo("wraps children in LayoutGroup");
});

describe("BentoCard", () => {
  it.todo("renders card with name, description, and icon");
  it.todo("applies hero size classes (col-span-2 row-span-2) on lg screens");
  it.todo("applies medium size classes (col-span-2) on sm screens");
  it.todo("applies compact size classes (col-span-1)");
  it.todo("calls onToggle when clicked");
  it.todo("renders children when isExpanded is true");
  it.todo("applies opacity 0.3 and scale 0.95 when isOtherExpanded");
  it.todo("shows close button when expanded");
});
