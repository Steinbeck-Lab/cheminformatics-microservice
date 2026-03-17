import { render, screen } from "@testing-library/react";
import { describe, expect, it, vi } from "vitest";
import { MemoryRouter, Route, Routes } from "react-router-dom";
import { AnimatedOutlet } from "@/components/common/AnimatedOutlet";

// Mock motion/react to render plain divs with data attributes for inspection
vi.mock("motion/react", () => ({
  AnimatePresence: ({ children }: { children: React.ReactNode }) => <>{children}</>,
  motion: {
    div: ({
      children,
      initial,
      animate,
      exit,
      ...rest
    }: {
      children?: React.ReactNode;
      initial?: Record<string, unknown>;
      animate?: Record<string, unknown>;
      exit?: Record<string, unknown>;
      [key: string]: unknown;
    }) => (
      <div
        data-testid="motion-div"
        data-initial={JSON.stringify(initial)}
        data-animate={JSON.stringify(animate)}
        data-exit={JSON.stringify(exit)}
        {...(rest.style ? { style: rest.style as React.CSSProperties } : {})}
      >
        {children}
      </div>
    ),
  },
}));

function renderWithRouter(initialPath: string) {
  return render(
    <MemoryRouter initialEntries={[initialPath]}>
      <Routes>
        <Route
          path="/"
          element={
            <div>
              <AnimatedOutlet />
            </div>
          }
        >
          <Route index element={<div data-testid="home-page">Home Page</div>} />
          <Route path="about" element={<div data-testid="about-page">About Page</div>} />
        </Route>
      </Routes>
    </MemoryRouter>
  );
}

describe("AnimatedOutlet", () => {
  it("renders without crashing", () => {
    renderWithRouter("/");
    expect(screen.getByTestId("home-page")).toBeInTheDocument();
  });

  it("wraps outlet content in a motion.div with correct initial/animate props", () => {
    renderWithRouter("/");
    const motionDiv = screen.getByTestId("motion-div");
    const initial = JSON.parse(motionDiv.getAttribute("data-initial") || "{}");
    const animate = JSON.parse(motionDiv.getAttribute("data-animate") || "{}");

    expect(initial).toEqual({ opacity: 0 });
    expect(animate).toEqual({ opacity: 1 });
  });

  it("sets exit animation with opacity fade", () => {
    renderWithRouter("/");
    const motionDiv = screen.getByTestId("motion-div");
    const exit = JSON.parse(motionDiv.getAttribute("data-exit") || "{}");

    expect(exit.opacity).toBe(0);
  });

  it("renders the correct child route for a nested path", () => {
    renderWithRouter("/about");
    expect(screen.getByTestId("about-page")).toBeInTheDocument();
    expect(screen.getByText("About Page")).toBeInTheDocument();
  });
});
