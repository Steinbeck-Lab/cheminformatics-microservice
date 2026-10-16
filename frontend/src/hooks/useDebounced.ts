import { useCallback, useEffect, useLayoutEffect, useRef, useState } from "react";

/**
 * Follows `value`, but only once it has stopped changing for `ms` milliseconds. Controls stay
 * live on the raw value while anything expensive (a server render, say) reads the settled one.
 *
 * The second item commits at once: with no argument it takes the latest value (on pointer
 * release, say); with one it sets both sides together (a reset that should not lag).
 */
export function useDebounced<T>(value: T, ms: number): [T, (...next: [] | [T]) => void] {
  const [settled, setSettled] = useState(value);
  const latest = useRef(value);

  useLayoutEffect(() => {
    latest.current = value;
  }, [value]);

  useEffect(() => {
    const timer = setTimeout(() => setSettled(() => value), ms);
    return () => clearTimeout(timer);
  }, [value, ms]);

  const commit = useCallback((...next: [] | [T]) => {
    const v = next.length === 1 ? next[0] : latest.current;
    setSettled(() => v);
  }, []);

  return [settled, commit];
}
