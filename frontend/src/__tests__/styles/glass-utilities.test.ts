import { describe, it, expect } from "vitest";

describe("Glass utility classes", () => {
  describe("glass-bold", () => {
    it.todo("applies backdrop-filter blur(24px) in light mode");
    it.todo("applies backdrop-filter blur(24px) in dark mode");
    it.todo("applies background-color with white/10 transparency");
    it.todo("includes -webkit-backdrop-filter for Safari");
    it.todo("includes border with white/30 transparency");
  });

  describe("glass-mobile", () => {
    it.todo("applies reduced backdrop-filter blur(12px)");
    it.todo("applies background-color with white/10 transparency");
  });

  describe("glass-refraction-hover", () => {
    it.todo("only activates on hover-capable devices");
    it.todo("applies SVG filter url on hover");
    it.todo("applies subtle scale(1.005) on hover");
  });
});
