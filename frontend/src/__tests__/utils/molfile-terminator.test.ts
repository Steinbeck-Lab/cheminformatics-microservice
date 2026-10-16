import { readFileSync, readdirSync, statSync } from "node:fs";
import { join, resolve } from "node:path";
import { describe, it, expect } from "vitest";

// The V2000 terminator is "M  END" with TWO spaces. A one-space check rejects every real molfile.
const SRC = resolve(__dirname, "../..");
const files = (dir: string): string[] =>
  readdirSync(dir).flatMap((n) => {
    const p = join(dir, n);
    if (statSync(p).isDirectory()) return n === "__tests__" ? [] : files(p);
    return /\.(ts|tsx)$/.test(n) ? [p] : [];
  });

describe("molfile terminator", () => {
  it('no source checks for "M END" with one space', () => {
    // a one-space marker inside a string: "M END", "M END\n" or the end of an example block
    const bad = files(SRC).filter((f) => /M END(\\n|"|`)/.test(readFileSync(f, "utf8")));
    expect(bad.map((f) => f.replace(SRC, ""))).toEqual([]);
  });

  it("example molfiles keep their fixed-width columns", () => {
    // a squeezed counts line looks like "4 3 0 0 0 0 0 0 0 0999 V2000"
    const bad = files(SRC).filter((f) =>
      /\d \d \d \d \d \d \d \d \d 0999 V2000/.test(readFileSync(f, "utf8"))
    );
    expect(bad.map((f) => f.replace(SRC, ""))).toEqual([]);
  });
});
