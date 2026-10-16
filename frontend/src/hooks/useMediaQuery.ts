import { useState, useEffect } from "react";

/**
 * Hook for responsive breakpoint detection using window.matchMedia.
 *
 * Starts from the current match (client-only app, no SSR), so the first render is already
 * right, then listens for changes.
 *
 * @param query - CSS media query string (e.g., '(max-width: 767px)')
 * @returns Whether the media query currently matches
 */
export function useMediaQuery(query: string): boolean {
  const [matches, setMatches] = useState(() => window.matchMedia(query).matches);

  useEffect(() => {
    const mql = window.matchMedia(query);
    setMatches(mql.matches);
    const handler = (e: MediaQueryListEvent) => setMatches(e.matches);
    mql.addEventListener("change", handler);
    return () => mql.removeEventListener("change", handler);
  }, [query]);

  return matches;
}
