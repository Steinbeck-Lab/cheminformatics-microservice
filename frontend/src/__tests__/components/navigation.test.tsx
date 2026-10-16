import { render, screen, within, fireEvent, act } from "@testing-library/react";
import { describe, it, expect } from "vitest";
import {
  createMemoryRouter,
  RouterProvider,
  createRoutesFromElements,
  Route,
} from "react-router-dom";
import { AppProvider } from "../../context/AppContext";
import Navigation from "../../components/common/Navigation";
import Header from "../../components/common/Header";

function renderAt(element: React.ReactElement, initialPath = "/chem") {
  const router = createMemoryRouter(
    createRoutesFromElements(<Route path="*" element={element} />),
    { initialEntries: [initialPath] }
  );
  return render(
    <AppProvider>
      <RouterProvider router={router} />
    </AppProvider>
  );
}

const SECTIONS = ["Home", "Analysis", "Convert", "Depict", "Tools", "OCSR", "About"];

describe("Navigation", () => {
  it("desktop bar lists the seven sections and marks the current one", () => {
    renderAt(<Navigation />, "/chem/descriptors");
    const nav = screen.getByRole("navigation", { name: "Sections" });
    const links = within(nav).getAllByRole("link");
    expect(links.map((l) => l.textContent)).toEqual(SECTIONS);
    expect(within(nav).getByRole("link", { name: "Analysis" })).toHaveAttribute(
      "aria-current",
      "page"
    );
    expect(within(nav).getByRole("link", { name: "Depict" })).not.toHaveAttribute("aria-current");
  });

  it("Home is current only on the home page", () => {
    renderAt(<Navigation />, "/home");
    expect(screen.getByRole("link", { name: "Home" })).toHaveAttribute("aria-current", "page");
  });

  it("menu variant has its own label and calls closeMenu on click", () => {
    let closed = 0;
    renderAt(<Navigation isMobile closeMenu={() => closed++} />, "/convert");
    const nav = screen.getByRole("navigation", { name: "Sections (menu)" });
    const convert = within(nav).getByRole("link", { name: "Convert" });
    expect(convert).toHaveAttribute("aria-current", "page");
    fireEvent.click(convert);
    expect(closed).toBe(1);
  });
});

describe("Header", () => {
  it("theme switch toggles dark mode", () => {
    renderAt(<Header />, "/depict");
    const sw = screen.getAllByRole("switch", { name: "Dark mode" })[0]!;
    const before = document.documentElement.classList.contains("dark");
    expect(sw).toHaveAttribute("aria-checked", String(before));
    fireEvent.click(sw);
    expect(document.documentElement.classList.contains("dark")).toBe(!before);
    expect(sw).toHaveAttribute("aria-checked", String(!before));
  });

  it("turns to clear glass once the page scrolls", () => {
    const { container } = renderAt(<Header />, "/depict");
    const header = container.querySelector("header")!;
    expect(header).not.toHaveClass("scrolled");
    act(() => {
      window.scrollY = 200;
      window.dispatchEvent(new Event("scroll"));
    });
    expect(header).toHaveClass("scrolled");
    act(() => {
      window.scrollY = 0;
      window.dispatchEvent(new Event("scroll"));
    });
    expect(header).not.toHaveClass("scrolled");
  });

  it("has a home link with the logo and a menu button for small screens", () => {
    renderAt(<Header />, "/depict");
    expect(
      screen.getByRole("link", { name: /Cheminformatics Microservice, home/ })
    ).toHaveAttribute("href", "/");
    expect(screen.getByRole("button", { name: "Open menu" })).toBeInTheDocument();
  });
});
