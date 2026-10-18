import { describe, it, expect, vi } from "vitest";
import { fireEvent, render, screen, within } from "@testing-library/react";
import { MemoryRouter } from "react-router-dom";
import { AppProvider } from "@/context/AppContext";
import HomePage from "@/pages/HomePage";

vi.mock("@/services/convertService", () => ({
  generate3DCoordinates: vi.fn(async () => "MOLBLOCK"),
}));

function renderHome() {
  return render(
    <AppProvider>
      <MemoryRouter>
        <HomePage />
      </MemoryRouter>
    </AppProvider>
  );
}

describe("Home page", () => {
  it("says what it is in a few words and offers a way in", () => {
    renderHome();
    const h1 = screen.getByRole("heading", { level: 1 });
    expect(h1.textContent!.trim().split(/\s+/).length).toBeLessThanOrEqual(4);
    const lead = document.querySelector(".pc-hero-lead")!;
    expect(lead.textContent!.trim().split(/\s+/).length).toBeLessThanOrEqual(6);
    expect(screen.getByRole("link", { name: /Start with Depict/ })).toHaveAttribute(
      "href",
      "/depict/2ddepiction"
    );
    expect(screen.getByRole("link", { name: /^API docs/ })).toHaveAttribute(
      "href",
      "https://api.naturalproducts.net/latest/docs"
    );
  });

  it("shows the molecule in 3D inside the lens and switches it", () => {
    renderHome();
    expect(screen.getByRole("img", { name: "Caffeine in 3D" })).toBeInTheDocument();
    fireEvent.click(screen.getByRole("button", { name: "Morphine" }));
    expect(screen.getByRole("img", { name: "Morphine in 3D" })).toBeInTheDocument();
    expect(screen.getByRole("button", { name: "Morphine" })).toHaveAttribute(
      "aria-pressed",
      "true"
    );
  });

  it("circles the lens with 2D drawings that open in Depict", () => {
    renderHome();
    const ring = screen.getByRole("list", { name: "2D drawings" });
    expect(within(ring).getAllByRole("link")).toHaveLength(8);
    for (const a of within(ring).getAllByRole("link")) {
      expect(a.getAttribute("href")).toMatch(/^\/depict\/2ddepiction\?smiles=/);
    }
  });

  it("has no moving row of structures under the hero", () => {
    renderHome();
    expect(screen.queryByRole("list", { name: "Structures drawn live by the API" })).toBeNull();
  });

  it("links the five sections by name only", () => {
    renderHome();
    const index = screen.getByRole("navigation", { name: "Sections" });
    const links = within(index).getAllByRole("link");
    expect(links.map((a) => a.getAttribute("href"))).toEqual([
      "/chem",
      "/convert",
      "/depict",
      "/tools",
      "/ocsr",
    ]);
    expect(links.map((a) => a.textContent)).toEqual([
      "Analysis",
      "Convert",
      "Depict",
      "Tools",
      "OCSR",
    ]);
  });

  it("ends with the sections, no extra proof strip", () => {
    renderHome();
    expect(screen.queryByRole("link", { name: /How to cite/ })).toBeNull();
    expect(screen.queryByText(/Peer reviewed/)).toBeNull();
  });

  it("closes with the open-source band and its credits", () => {
    renderHome();
    expect(
      screen.getByRole("heading", {
        level: 2,
        name: "Designed for Chemists. Fueled by Open-Source.",
      })
    ).toBeInTheDocument();
    expect(screen.getByRole("link", { name: /^Steinbeck Lab/ })).toHaveAttribute(
      "href",
      "https://cheminf.uni-jena.de"
    );
    expect(
      screen.getByRole("link", { name: /^Friedrich Schiller University Jena/ })
    ).toHaveAttribute("href", "https://www.uni-jena.de/en");
  });

  it("keeps the words short: under 120 words in the page body", () => {
    const { container } = renderHome();
    const words = (container.textContent ?? "").split(/\s+/).filter(Boolean).length;
    expect(words).toBeLessThan(120);
  });
});
