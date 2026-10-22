import { useEffect } from "react";

/**
 * Scroll the bento cell whose title has this id into view whenever one of `triggers` turns on
 * (pass the run's loading flag, its result and its error). Each run flips the loading flag, so
 * every run scrolls, not only the first. Cells already in the upper part of the screen stay put.
 */
export function useRevealOnResult(id: string, ...triggers: unknown[]) {
  const key = triggers.map((t) => (t ? 1 : 0)).join("");
  useEffect(() => {
    if (!key.includes("1")) return;
    const cell = document.getElementById(id)?.closest("section");
    if (!cell || cell.getBoundingClientRect().top < window.innerHeight * 0.6) return;
    const reduce = window.matchMedia?.("(prefers-reduced-motion: reduce)").matches;
    cell.scrollIntoView?.({ behavior: reduce ? "auto" : "smooth", block: "start" });
  }, [id, key]);
}
