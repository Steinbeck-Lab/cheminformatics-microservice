import { useState, useEffect } from "react";

/**
 * Hook for responsive breakpoint detection using window.matchMedia.
 *
 * Initializes to false and updates after mount to avoid SSR mismatch.
 * Listens for changes and updates the boolean state accordingly.
 *
 * @param query - CSS media query string (e.g., '(max-width: 767px)')
 * @returns Whether the media query currently matches
 */
export function useMediaQuery(query: string): boolean {
  const [matches, setMatches] = useState(false);

  useEffect(() => {
    const mql = window.matchMedia(query);
    setMatches(mql.matches);
    const handler = (e: MediaQueryListEvent) => setMatches(e.matches);
    mql.addEventListener("change", handler);
    return () => mql.removeEventListener("change", handler);
  }, [query]);

  return matches;
}
