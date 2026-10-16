import { describe, it, expect, vi, afterEach } from "vitest";
import { act, render } from "@testing-library/react";
import { silkRibbons } from "@/lib/silk";
import { SilkBackground } from "@/components/common/SilkBackground";

describe("silkRibbons", () => {
  it("returns three closed ribbons with six warp lines each", () => {
    const ribbons = silkRibbons(1280, 800);
    expect(ribbons.map((r) => r.id)).toEqual([1, 2, 3]);
    for (const r of ribbons) {
      expect(r.fill).toMatch(/^M-?\d/);
      expect(r.fill.endsWith("Z")).toBe(true);
      expect(r.fill.split("L")).toHaveLength(57 * 2); // 57 points per edge
      expect(r.warp.match(/M/g)).toHaveLength(6);
      expect(r.fill + r.warp).not.toContain("NaN");
    }
  });

  it("starts the grey ribbon off the left edge", () => {
    // first control point: x = -0.15 * W, y = 0.02 * H, pushed up by half the ribbon width
    const [first] = silkRibbons(1000, 1000);
    const x = parseFloat(first.fill.slice(1));
    expect(x).toBeLessThan(-100);
  });
});

describe("silkRibbons over time", () => {
  it("is the same shape at t = 0 and ripples later", () => {
    const still = silkRibbons(1280, 800);
    expect(silkRibbons(1280, 800, 0)).toEqual(still);
    const later = silkRibbons(1280, 800, 3.7);
    expect(later[0].fill).not.toEqual(still[0].fill);
    expect(later[0].fill + later[0].warp).not.toContain("NaN");
  });

  it("loops: the shape repeats after one full period", () => {
    const a = silkRibbons(1000, 700, 1.3);
    const b = silkRibbons(1000, 700, 1.3 + 2 * Math.PI);
    expect(b[1].fill).toEqual(a[1].fill);
  });
});

describe("SilkBackground", () => {
  afterEach(() => {
    vi.useRealTimers();
    window.innerWidth = 1024; // jsdom defaults
    window.innerHeight = 768;
  });

  it("draws three hidden ribbon layers at the viewport size", () => {
    const { container } = render(<SilkBackground />);
    const root = container.firstElementChild!;
    expect(root).toHaveAttribute("aria-hidden", "true");
    const svgs = root.querySelectorAll("svg");
    expect(svgs).toHaveLength(3);
    expect(svgs[0]!.getAttribute("viewBox")).toBe(`0 0 ${window.innerWidth} ${window.innerHeight}`);
  });

  it("redraws after the window is resized", () => {
    vi.useFakeTimers();
    const { container } = render(<SilkBackground />);
    act(() => {
      window.innerWidth = 900;
      window.innerHeight = 1400;
      window.dispatchEvent(new Event("resize"));
      vi.advanceTimersByTime(250);
    });
    expect(container.querySelector("svg")!.getAttribute("viewBox")).toBe("0 0 900 1400");
  });

  it("never draws smaller than 320 × 480", () => {
    vi.useFakeTimers();
    const { container } = render(<SilkBackground />);
    act(() => {
      window.innerWidth = 200;
      window.innerHeight = 300;
      window.dispatchEvent(new Event("resize"));
      vi.advanceTimersByTime(250);
    });
    expect(container.querySelector("svg")!.getAttribute("viewBox")).toBe("0 0 320 480");
  });
});

describe("SilkBackground live mode", () => {
  it("marks the silk as live only when asked", () => {
    const { container, rerender } = render(<SilkBackground />);
    expect(container.firstElementChild).not.toHaveAttribute("data-live");
    rerender(<SilkBackground live />);
    expect(container.firstElementChild).toHaveAttribute("data-live");
  });
});
