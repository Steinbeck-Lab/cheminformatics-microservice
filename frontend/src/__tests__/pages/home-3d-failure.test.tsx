import { describe, it, expect, vi } from "vitest";
import { act, render, screen } from "@testing-library/react";
import { MemoryRouter } from "react-router-dom";
import { AppProvider } from "@/context/AppContext";
import HomePage from "@/pages/HomePage";

// the 3D chunk fails to load (flaky network, or a deploy replaced the chunk)
vi.mock("@/components/common/Molecule3D", () => {
  throw new Error("Failed to fetch dynamically imported module");
});

describe("Home when the 3D code cannot load", () => {
  it("keeps the page and the lens; only the 3D model is missing", async () => {
    vi.spyOn(console, "error").mockImplementation(() => {});
    render(
      <AppProvider>
        <MemoryRouter>
          <HomePage />
        </MemoryRouter>
      </AppProvider>
    );
    await act(async () => {
      await new Promise((r) => setTimeout(r, 50));
    });
    expect(screen.getByRole("heading", { level: 1 })).toBeInTheDocument();
    expect(screen.getByRole("img", { name: "Caffeine in 3D" })).toBeInTheDocument();
  });
});
