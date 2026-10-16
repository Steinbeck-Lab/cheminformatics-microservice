import { describe, it, expect } from "vitest";
import { parseMolblock, assembly } from "@/lib/molblock";

// ethanol-like V2000 block: C C O, two bonds, one 3D conformer
const MOL = `
     RDKit          3D

  3  2  0  0  0  0  0  0  0  0999 V2000
   -1.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.5000    0.2000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.1000    1.4000   -0.3000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  2  0
M  END
`;

describe("parseMolblock", () => {
  it("reads atoms and bonds, centred on their mean", () => {
    const m = parseMolblock(MOL);
    expect(m.atoms.map((a) => a.el)).toEqual(["C", "C", "O"]);
    expect(m.bonds).toEqual([
      { a: 0, b: 1, order: 1 },
      { a: 1, b: 2, order: 2 },
    ]);
    const meanOf = (k: 0 | 1 | 2) => m.atoms.reduce((s, a) => s + a.pos[k], 0) / m.atoms.length;
    [meanOf(0), meanOf(1), meanOf(2)].forEach((v) => expect(Math.abs(v)).toBeLessThan(1e-9));
  });

  it("accepts the API's JSON-quoted answer", () => {
    expect(parseMolblock(JSON.stringify(MOL)).atoms).toHaveLength(3);
  });

  it("returns nothing for text that is not a molblock", () => {
    expect(parseMolblock("Internal Server Error")).toEqual({ atoms: [], bonds: [] });
  });
});

describe("assembly", () => {
  it("brings each atom in turn from 0 to 1, then each bond", () => {
    const n = { atoms: 10, bonds: 9 };
    expect(assembly(0, n).atom(0)).toBe(0);
    expect(assembly(0.4, n).atom(0)).toBeGreaterThan(assembly(0.4, n).atom(9));
    const done = assembly(99, n);
    expect(done.atom(9)).toBe(1);
    expect(done.bond(8)).toBe(1);
    expect(done.finished).toBe(true);
    expect(assembly(0.2, n).bond(0)).toBe(0); // bonds wait for the atoms
  });
});
