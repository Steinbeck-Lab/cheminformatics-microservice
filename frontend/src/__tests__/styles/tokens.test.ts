import { readFileSync } from "node:fs";
import { resolve } from "node:path";
import { describe, it, expect } from "vitest";

const read = (file: string) => readFileSync(resolve(__dirname, "../../styles", file), "utf8");
// tailwind.css holds the base tokens; porcelain.css adds the second tier (--link, --el-h, ...)
const sheets = [read("tailwind.css"), read("porcelain.css")];

// Hex tokens of the first ":root" or ".dark" block of each sheet, merged.
function tokens(selector: string): Record<string, string> {
  const esc = selector.replace(".", "\\.");
  const out: Record<string, string> = {};
  for (const css of sheets) {
    const m = css.match(new RegExp(`(?:^|\\n)\\s*${esc}\\s*\\{([^}]*)\\}`));
    if (!m) throw new Error(`no ${selector} block`);
    for (const x of (m[1] ?? "").matchAll(/--([\w-]+):\s*(#[0-9a-fA-F]{6})\b/g)) {
      if (x[1] && x[2]) out[x[1]] = x[2];
    }
  }
  return out;
}

// One sRGB channel (hex digits at i, i + 1) as linear light.
function channel(hex: string, i: number): number {
  const c = parseInt(hex.slice(i, i + 2), 16) / 255;
  return c <= 0.04045 ? c / 12.92 : ((c + 0.055) / 1.055) ** 2.4;
}

function luminance(hex: string): number {
  return 0.2126 * channel(hex, 1) + 0.7152 * channel(hex, 3) + 0.0722 * channel(hex, 5);
}

function contrast(a: string, b: string): number {
  const [la, lb] = [luminance(a), luminance(b)];
  return (Math.max(la, lb) + 0.05) / (Math.min(la, lb) + 0.05);
}

// A token's value; a missing token fails the test.
function tok(t: Record<string, string>, name: string): string {
  const v = t[name];
  if (v === undefined) throw new Error(`--${name} missing`);
  return v;
}

// [text token, surface token, minimum ratio]
const PAIRS: [string, string, number][] = [
  ["foreground", "background", 4.5],
  ["muted-foreground", "background", 4.5],
  ["muted-foreground", "glass-solid", 4.5],
  ["head", "background", 4.5],
  ["primary-foreground", "primary", 4.5],
  ["secondary-foreground", "secondary", 4.5],
  ["saffron-ink", "glass-solid", 4.5],
  ["madder-ink", "glass-solid", 4.5],
  ["danger-ink", "danger", 4.5],
  ["destructive-foreground", "destructive", 4.5],
  ["madder-ink", "secondary", 4.5],
  ["ring", "background", 3],
  ["link", "background", 4.5],
  ["link", "glass-solid", 4.5],
  // a hydrogen atom in the home 3D molecule: a pale sphere that must still part from the page
  ["el-h", "background", 1.75],
];

describe.each([":root", ".dark"])("Porcelain tokens in %s", (selector) => {
  const t = tokens(selector);

  it.each(PAIRS)("%s on %s meets %d:1", (fg, bg, min) => {
    expect(contrast(tok(t, fg), tok(t, bg))).toBeGreaterThanOrEqual(min);
  });

  it("uses the brand saffron fill", () => {
    expect(t.saffron?.toLowerCase()).toBe("#f2b31a");
  });

  it("keeps the depiction stage white", () => {
    expect(t.stage?.toLowerCase()).toBe("#ffffff");
  });
});
