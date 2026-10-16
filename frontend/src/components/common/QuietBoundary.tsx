import { Component, type ReactNode } from "react";

/**
 * Catches a failure in decorative content (a lazy chunk that will not load, a WebGL crash) and
 * renders the fallback instead, so the rest of the page keeps working.
 */
export class QuietBoundary extends Component<
  { children: ReactNode; fallback?: ReactNode },
  { failed: boolean }
> {
  state = { failed: false };

  static getDerivedStateFromError() {
    return { failed: true };
  }

  render() {
    return this.state.failed ? (this.props.fallback ?? null) : this.props.children;
  }
}
