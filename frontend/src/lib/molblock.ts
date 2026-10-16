type Vec3 = [number, number, number];

/** A molecule for drawing in 3D: atoms centred on their mean, bonds by atom index. */
export interface Mol3 {
  atoms: { el: string; pos: Vec3 }[];
  bonds: { a: number; b: number; order: number }[];
}

/**
 * Read a V2000 molblock (the API's /convert/mol3D answer, possibly JSON-quoted).
 * Returns no atoms for anything it cannot read.
 */
export function parseMolblock(text: string): Mol3 {
  let src = text;
  if (src.trimStart().startsWith('"')) {
    try {
      src = JSON.parse(src);
    } catch {
      return { atoms: [], bonds: [] };
    }
  }
  const lines = src.split(/\r?\n/);
  const countsAt = lines.findIndex((l) => /V2000\s*$/.test(l));
  if (countsAt < 0) return { atoms: [], bonds: [] };
  const nAtoms = parseInt(lines[countsAt]!.slice(0, 3), 10);
  const nBonds = parseInt(lines[countsAt]!.slice(3, 6), 10);
  const atoms: Mol3["atoms"] = [];
  for (let i = 0; i < nAtoms; i++) {
    const l = lines[countsAt + 1 + i] ?? "";
    const pos: Vec3 = [
      parseFloat(l.slice(0, 10)),
      parseFloat(l.slice(10, 20)),
      parseFloat(l.slice(20, 30)),
    ];
    const el = l.slice(31, 34).trim();
    if (!el || pos.some(Number.isNaN)) return { atoms: [], bonds: [] };
    atoms.push({ el, pos });
  }
  const bonds: Mol3["bonds"] = [];
  for (let j = 0; j < nBonds; j++) {
    const l = lines[countsAt + 1 + nAtoms + j] ?? "";
    const a = parseInt(l.slice(0, 3), 10) - 1;
    const b = parseInt(l.slice(3, 6), 10) - 1;
    const order = parseInt(l.slice(6, 9), 10);
    if (atoms[a] && atoms[b]) bonds.push({ a, b, order });
  }
  // centre on the mean so the molecule turns about its middle
  const meanOf = (k: 0 | 1 | 2) => atoms.reduce((s, x) => s + x.pos[k], 0) / (atoms.length || 1);
  const mean: Vec3 = [meanOf(0), meanOf(1), meanOf(2)];
  for (const x of atoms) x.pos = [x.pos[0] - mean[0], x.pos[1] - mean[1], x.pos[2] - mean[2]];
  return { atoms, bonds };
}

const ATOM_TIME = 0.75; // seconds for one atom to fly in
const BOND_TIME = 0.35; // seconds for one bond to grow
const expoOut = (x: number) => (x >= 1 ? 1 : 1 - Math.pow(2, -10 * x));
const clamp01 = (x: number) => Math.min(1, Math.max(0, x));

/**
 * The assembly at time t (seconds): atoms fly in one after another, then the bonds grow.
 * atom(i) and bond(j) are progress from 0 to 1, eased.
 */
export function assembly(t: number, n: { atoms: number; bonds: number }) {
  const atomStagger = Math.min(0.07, 1.1 / Math.max(1, n.atoms));
  const bondStagger = Math.min(0.05, 0.7 / Math.max(1, n.bonds));
  const atomsDone = Math.max(0, n.atoms - 1) * atomStagger + ATOM_TIME;
  const end = atomsDone + Math.max(0, n.bonds - 1) * bondStagger + BOND_TIME;
  return {
    atom: (i: number) => expoOut(clamp01((t - i * atomStagger) / ATOM_TIME)),
    bond: (j: number) => expoOut(clamp01((t - atomsDone - j * bondStagger) / BOND_TIME)),
    finished: t >= end,
  };
}
