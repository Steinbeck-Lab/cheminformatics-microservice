import { readFileSync } from "node:fs";
import { resolve } from "node:path";
import { describe, it, expect } from "vitest";

const css = readFileSync(resolve(__dirname, "../../styles/porcelain.css"), "utf8");

describe("porcelain.css class names", () => {
  // .pc-about once meant both "the About page" and "a tool's about paragraph";
  // the second rule broke the page's centring. One class, one top-level rule block.
  // Two-space indent = top level of a @layer (rules inside media queries are indented further).
  it.each(["pc-about", "pc-wrap"])(".%s is defined by one top-level rule only", (name) => {
    const blocks = css.match(new RegExp(`^  \\.${name}\\s*\\{`, "gm")) ?? [];
    expect(blocks).toHaveLength(1);
  });
});
