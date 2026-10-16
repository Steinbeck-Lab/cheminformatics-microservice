/**
 * "Silk crossings": three wide ribbons on cubic Bézier centrelines. The width swells and
 * narrows like draped cloth; six fine warp lines run along each ribbon.
 * Grey (1) and madder (2) are wide, saffron (3) is the narrowest.
 */
export interface SilkRibbon {
  id: 1 | 2 | 3;
  /** Closed outline of the ribbon. */
  fill: string;
  /** Six open warp lines inside it. */
  warp: string;
}

/** Cubic Bézier control points: x0 y0 x1 y1 x2 y2 x3 y3. */
type Curve = [number, number, number, number, number, number, number, number];

const STEPS = 56;
const WARP_OFFSETS = [-0.96, -0.55, -0.12, 0.34, 0.78, 0.96];
const f = (v: number) => v.toFixed(1);

function ribbon(p: Curve, width: number, phase: number, t = 0) {
  const up: string[] = [];
  const down: string[] = [];
  const warps: string[][] = WARP_OFFSETS.map(() => []);
  for (let k = 0; k <= STEPS; k++) {
    const u = k / STEPS;
    const mt = 1 - u;
    const x = mt ** 3 * p[0] + 3 * mt * mt * u * p[2] + 3 * mt * u * u * p[4] + u ** 3 * p[6];
    const y = mt ** 3 * p[1] + 3 * mt * mt * u * p[3] + 3 * mt * u * u * p[5] + u ** 3 * p[7];
    const dx = 3 * mt * mt * (p[2] - p[0]) + 6 * mt * u * (p[4] - p[2]) + 3 * u * u * (p[6] - p[4]);
    const dy = 3 * mt * mt * (p[3] - p[1]) + 6 * mt * u * (p[5] - p[3]) + 3 * u * u * (p[7] - p[5]);
    const len = Math.hypot(dx, dy) || 1;
    const nx = -dy / len;
    const ny = dx / len;
    const h = (width * (0.7 + 0.4 * Math.sin(Math.PI * u * 1.3 + phase + t))) / 2;
    up.push(`${f(x + nx * h)} ${f(y + ny * h)}`);
    down.push(`${f(x - nx * h)} ${f(y - ny * h)}`);
    WARP_OFFSETS.forEach((o, i) => warps[i]!.push(`${f(x + nx * h * o)} ${f(y + ny * h * o)}`));
  }
  return {
    fill: `M${up.join("L")}L${down.reverse().join("L")}Z`,
    warp: warps.map((w) => `M${w.join("L")}`).join(""),
  };
}

/**
 * @param t time in radians. Every motion is a whole-number multiple of t, so the shapes repeat
 *   exactly every 2π (a seamless loop), and t = 0 is the still drawing.
 */
export function silkRibbons(W: number, H: number, t = 0): [SilkRibbon, SilkRibbon, SilkRibbon] {
  const S = Math.sqrt(W * H);
  // the inner control points sway; zero at t = 0
  const sway = (amp: number, k: number) => amp * (Math.sin(t + k) - Math.sin(k));
  const curve = (c: Curve, k: number): Curve => [
    c[0],
    c[1],
    c[2] + sway(0.03 * W, k),
    c[3] + sway(0.06 * H, k + 1),
    c[4] + sway(0.03 * W, k + 2),
    c[5] + sway(0.06 * H, k + 3),
    c[6],
    c[7],
  ];
  return [
    {
      id: 1,
      ...ribbon(
        curve([-0.15 * W, 0.02 * H, 0.32 * W, 0.06 * H, 0.58 * W, 0.98 * H, 1.15 * W, 0.86 * H], 0),
        0.3 * S,
        0.2,
        t
      ),
    },
    {
      id: 2,
      ...ribbon(
        curve([-0.15 * W, 1.02 * H, 0.36 * W, 0.78 * H, 0.62 * W, 0.06 * H, 1.15 * W, 0.14 * H], 2),
        0.27 * S,
        1.4,
        -t
      ),
    },
    {
      id: 3,
      ...ribbon(
        curve([-0.15 * W, 0.6 * H, 0.3 * W, 0.3 * H, 0.7 * W, 0.8 * H, 1.15 * W, 0.36 * H], 4),
        0.105 * S,
        2.6,
        t
      ),
    },
  ];
}
