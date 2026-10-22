import { useEffect, useMemo, useRef, useState } from "react";
import { useMediaQuery } from "@/hooks/useMediaQuery";
import { silkRibbons } from "@/lib/silk";

const viewport = () => ({
  w: Math.max(320, Math.round(window.innerWidth)),
  h: Math.max(480, Math.round(window.innerHeight)),
});

const RIBBON_IDS = [1, 2, 3] as const;
const FPS = 30;
const SPEED = 0.16; // radians per second: one full ripple cycle takes about 40 s
const SHEEN_SPEED = 60; // px per second for the band of light

/**
 * Fixed page background: three translucent silk ribbons that deepen where they cross
 * (multiply in light, screen in dark), each lit by a soft sheen. Colours come from the
 * --silk-* tokens. `live` (home, about): the ribbons ripple and the sheen glides across them,
 * redrawn about 30 times a second; paused while the tab is hidden and off under reduced motion.
 */
export function SilkBackground({ live = false }: { live?: boolean }) {
  const [size, setSize] = useState(viewport);
  const fills = useRef<(SVGPathElement | null)[]>([]);
  const warps = useRef<(SVGPathElement | null)[]>([]);
  const sheens = useRef<(SVGLinearGradientElement | null)[]>([]);
  const reduce = useMediaQuery("(prefers-reduced-motion: reduce)"); // live, unlike motion's useReducedMotion

  useEffect(() => {
    let timer: ReturnType<typeof setTimeout>;
    const onResize = () => {
      clearTimeout(timer);
      timer = setTimeout(() => setSize(viewport()), 220);
    };
    window.addEventListener("resize", onResize);
    return () => {
      clearTimeout(timer);
      window.removeEventListener("resize", onResize);
    };
  }, []);

  const ribbons = useMemo(() => silkRibbons(size.w, size.h), [size]);

  useEffect(() => {
    if (!live || reduce) return;
    if (typeof requestAnimationFrame !== "function") return;
    let raf = 0;
    let last = 0;
    const start = performance.now();
    const frame = (now: number) => {
      raf = requestAnimationFrame(frame);
      if (document.hidden || now - last < 1000 / FPS) return;
      last = now;
      const s = (now - start) / 1000;
      silkRibbons(size.w, size.h, s * SPEED).forEach((r, i) => {
        fills.current[i]?.setAttribute("d", r.fill);
        warps.current[i]?.setAttribute("d", r.warp);
      });
      // gradients repeat, so shifting by any amount keeps the sheen seamless
      sheens.current.forEach((g, i) =>
        g?.setAttribute(
          "gradientTransform",
          `translate(${(s * SHEEN_SPEED * (i === 1 ? -1 : 1)) % (size.w * 2)} 0)`
        )
      );
    };
    raf = requestAnimationFrame(frame);
    return () => cancelAnimationFrame(raf);
  }, [live, size, reduce]);

  return (
    <div className="silk" aria-hidden="true" data-live={live || undefined}>
      {ribbons.map((r, i) => (
        <svg
          key={r.id}
          className={`silk-g${r.id}`}
          width={size.w}
          height={size.h}
          viewBox={`0 0 ${size.w} ${size.h}`}
          focusable="false"
        >
          <defs>
            <linearGradient
              id={`silk-sheen-${RIBBON_IDS[i]}`}
              ref={(el) => {
                sheens.current[i] = el;
              }}
              gradientUnits="userSpaceOnUse"
              x1={0}
              y1={0}
              x2={size.w * 0.9}
              y2={size.h * 0.35}
              spreadMethod="repeat"
            >
              <stop offset="0" className={`silk-stop silk-stop-${r.id}`} />
              <stop offset="0.42" className={`silk-stop silk-stop-${r.id}`} />
              <stop offset="0.5" className="silk-stop silk-stop-hi" />
              <stop offset="0.58" className={`silk-stop silk-stop-${r.id}`} />
              <stop offset="1" className={`silk-stop silk-stop-${r.id}`} />
            </linearGradient>
          </defs>
          <path
            ref={(el) => {
              fills.current[i] = el;
            }}
            className={`silk-r${r.id}`}
            fill={`url(#silk-sheen-${RIBBON_IDS[i]})`}
            d={r.fill}
          />
          <path
            ref={(el) => {
              warps.current[i] = el;
            }}
            className="silk-warp"
            d={r.warp}
          />
        </svg>
      ))}
    </div>
  );
}
