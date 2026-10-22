import { GlassSkeleton } from "./GlassSkeleton";

/**
 * Shown while a section's code loads: the same tabs, head and bento as a tool page, so the real
 * page lands in place. It fades in only after a short delay (fast loads never flash it) and a
 * silk sheen glides across the bars while the viewer allows motion.
 */
export function RouteLoadingFallback() {
  return (
    <div
      data-testid="route-loading-fallback"
      role="status"
      aria-label="Loading page"
      className="pc-wrap pc-skel"
    >
      <div className="pc-tabs glass pc-skel-tabs">
        <GlassSkeleton className="w-24" />
        <GlassSkeleton className="w-28" />
        <GlassSkeleton className="w-20" />
      </div>
      <div className="glass my-3 rounded-cell px-4 py-[18px] sm:px-6">
        <GlassSkeleton className="h-8 w-64" />
        <GlassSkeleton className="mt-3 w-96 max-w-full" />
      </div>
      <div className="pc-bento">
        <div className="pc-cell glass pc-t-input grid content-start gap-3">
          <GlassSkeleton className="w-32" />
          <GlassSkeleton variant="image" className="h-14" />
          <GlassSkeleton className="w-3/4" />
          <GlassSkeleton className="h-10 w-40 rounded-full" />
        </div>
        <div className="pc-cell glass pc-t-side grid content-start gap-3">
          <GlassSkeleton className="w-24" />
          <GlassSkeleton />
          <GlassSkeleton className="w-5/6" />
          <GlassSkeleton className="w-2/3" />
        </div>
        <div className="pc-cell glass pc-t-wide">
          <GlassSkeleton variant="image" className="h-20" />
        </div>
      </div>
    </div>
  );
}
