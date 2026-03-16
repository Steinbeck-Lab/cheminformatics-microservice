import { describe, it, expect, vi } from "vitest";
import { render, screen, fireEvent, act } from "@testing-library/react";

vi.mock("@/services/depictService", () => ({
  get2DDepictionUrl: (smiles: string) =>
    `https://api.test/depict/2D?smiles=${encodeURIComponent(smiles)}&width=240&height=240`,
}));

import { SMILESPreview } from "@/components/common/SMILESPreview";

/**
 * Helper: wraps test body with fake timers, restoring real timers before
 * the function returns so that the global afterEach cleanup() (from setup.ts)
 * runs with real timers and does not hang on React effect teardowns.
 */
function withFakeTimers(fn: () => void) {
  vi.useFakeTimers();
  try {
    fn();
  } finally {
    vi.useRealTimers();
  }
}

describe("useDebounce", () => {
  it("returns initial value immediately, then debounces subsequent changes", () =>
    withFakeTimers(() => {
      const { rerender } = render(<SMILESPreview smiles="" />);
      expect(screen.queryByTestId("smiles-preview")).not.toBeInTheDocument();

      rerender(<SMILESPreview smiles="CCO" />);

      // Before debounce fires, preview not visible (debouncedSmiles still "")
      expect(screen.queryByTestId("smiles-preview")).not.toBeInTheDocument();

      act(() => {
        vi.advanceTimersByTime(500);
      });

      // After debounce, preview appears
      expect(screen.getByTestId("smiles-preview")).toBeInTheDocument();
    }));
});

describe("SMILESPreview", () => {
  it("renders nothing when smiles is empty", () =>
    withFakeTimers(() => {
      render(<SMILESPreview smiles="" />);
      expect(screen.queryByTestId("smiles-preview")).not.toBeInTheDocument();
    }));

  it("renders nothing when smiles is whitespace only", () =>
    withFakeTimers(() => {
      render(<SMILESPreview smiles="   " />);

      act(() => {
        vi.advanceTimersByTime(500);
      });

      expect(screen.queryByTestId("smiles-preview")).not.toBeInTheDocument();
    }));

  it("shows a glass skeleton while image is loading", () =>
    withFakeTimers(() => {
      const { rerender } = render(<SMILESPreview smiles="" />);
      rerender(<SMILESPreview smiles="CCO" />);

      act(() => {
        vi.advanceTimersByTime(500);
      });

      expect(screen.getByTestId("glass-skeleton")).toBeInTheDocument();

      const img = screen.getByAltText("Structure preview");
      expect(img).toBeInTheDocument();
      expect(img.className).toContain("opacity-0");
    }));

  it("shows the structure image when loaded successfully", () =>
    withFakeTimers(() => {
      const { rerender } = render(<SMILESPreview smiles="" />);
      rerender(<SMILESPreview smiles="CCO" />);

      act(() => {
        vi.advanceTimersByTime(500);
      });

      const img = screen.getByAltText("Structure preview");

      act(() => {
        fireEvent.load(img);
      });

      expect(img.className).toContain("opacity-100");
      expect(screen.queryByTestId("glass-skeleton")).not.toBeInTheDocument();
    }));

  it("hides silently when image fails to load (onError)", () =>
    withFakeTimers(() => {
      const { rerender } = render(<SMILESPreview smiles="" />);
      rerender(<SMILESPreview smiles="invalid_smiles_xyz" />);

      act(() => {
        vi.advanceTimersByTime(500);
      });

      const img = screen.getByAltText("Structure preview");

      act(() => {
        fireEvent.error(img);
      });

      expect(screen.queryByTestId("smiles-preview")).not.toBeInTheDocument();
    }));

  it("debounces -- does not update image URL on every keystroke", () =>
    withFakeTimers(() => {
      const { rerender } = render(<SMILESPreview smiles="" />);

      rerender(<SMILESPreview smiles="C" />);
      rerender(<SMILESPreview smiles="CC" />);
      rerender(<SMILESPreview smiles="CCO" />);

      act(() => {
        vi.advanceTimersByTime(300);
      });

      expect(screen.queryByTestId("smiles-preview")).not.toBeInTheDocument();

      act(() => {
        vi.advanceTimersByTime(200);
      });

      expect(screen.getByTestId("smiles-preview")).toBeInTheDocument();
      const img = screen.getByAltText("Structure preview");
      expect(img.getAttribute("src")).toContain("CCO");
    }));
});
