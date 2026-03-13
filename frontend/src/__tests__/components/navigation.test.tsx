import { render, screen } from "@testing-library/react";
import { describe, it, expect } from "vitest";
import {
  createMemoryRouter,
  RouterProvider,
  createRoutesFromElements,
  Route,
} from "react-router-dom";
import { AppProvider } from "../../context/AppContext";
import Navigation from "../../components/common/Navigation";

/**
 * Helpers to render Navigation within a memory router at a given path.
 */
function renderNavigation({
  initialPath = "/chem",
  isMobile = false,
}: { initialPath?: string; isMobile?: boolean } = {}) {
  const router = createMemoryRouter(
    createRoutesFromElements(
      <Route path="*" element={<Navigation isMobile={isMobile} closeMenu={() => {}} />} />
    ),
    { initialEntries: [initialPath] }
  );

  return render(
    <AppProvider>
      <RouterProvider router={router} />
    </AppProvider>
  );
}

describe("Navigation wayfinding", () => {
  it("active desktop nav link has distinct visual indicator", () => {
    renderNavigation({ initialPath: "/chem", isMobile: false });

    // Find the Chemical Analysis link -- the NavLink wraps a motion.div with active styling
    const chemLink = screen.getByText("Chemical Analysis");
    // Walk up from the text to find the element with active styling (text-white)
    const el: HTMLElement | null = chemLink.closest("a");
    expect(el).toBeTruthy();
    // Inside the NavLink <a>, there should be a div with "text-white" for the active state
    const activeDiv = el!.querySelector("[class*='text-white']");
    expect(activeDiv).toBeTruthy();
  });

  it("mobile nav shows current page highlighted with active styling", () => {
    renderNavigation({ initialPath: "/chem", isMobile: true });

    // In mobile mode, the full label is shown
    const chemLink = screen.getByText("Chemical Analysis");
    // The NavLink should have the active class applied
    const navLink = chemLink.closest("a");
    expect(navLink).toBeTruthy();
    expect(navLink!.className).toMatch(/bg-slate-900/);
    expect(navLink!.className).toMatch(/font-semibold/);
    expect(navLink!.className).toMatch(/border-l-4/);
  });

  it("mobile nav items have minimum 44px tap target", () => {
    renderNavigation({ initialPath: "/chem", isMobile: true });

    const links = screen.getAllByRole("link");
    for (const link of links) {
      expect(link.className).toContain("min-h-[44px]");
    }
  });
});
