/**
 * Hidden SVG filter definition for liquid glass refraction effect.
 *
 * This component renders a zero-size SVG containing the feTurbulence +
 * feDisplacementMap filter used by the .glass-refraction-hover CSS class.
 * It should be rendered once in the app layout (e.g., App.tsx or a layout wrapper).
 *
 * The filter creates a subtle organic distortion simulating glass refraction.
 * It is applied on hover via CSS and restricted to desktop devices only.
 */
export function LiquidGlassFilter() {
  return (
    <svg width="0" height="0" style={{ position: "absolute" }} aria-hidden="true">
      <defs>
        <filter id="liquid-glass-refraction" x="-10%" y="-10%" width="120%" height="120%">
          <feTurbulence
            type="fractalNoise"
            baseFrequency="0.015"
            numOctaves={3}
            seed={2}
            result="turbulence"
          />
          <feDisplacementMap
            in="SourceGraphic"
            in2="turbulence"
            scale={6}
            xChannelSelector="R"
            yChannelSelector="G"
          />
        </filter>
      </defs>
    </svg>
  );
}
