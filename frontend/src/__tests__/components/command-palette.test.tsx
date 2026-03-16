import React from "react";
import { render, screen, fireEvent, waitFor } from "@testing-library/react";
import { describe, it, expect, vi, beforeEach } from "vitest";
import {
  createMemoryRouter,
  RouterProvider,
  createRoutesFromElements,
  Route,
} from "react-router-dom";
import { AppProvider } from "../../context/AppContext";
import { CommandPalette } from "../../components/common/CommandPalette";

// Track navigate calls
const navigateMock = vi.fn();

vi.mock("react-router-dom", async () => {
  const actual = await vi.importActual("react-router-dom");
  return {
    ...actual,
    useNavigate: () => navigateMock,
  };
});

function renderPalette() {
  const router = createMemoryRouter(
    createRoutesFromElements(<Route path="*" element={<CommandPalette />} />),
    { initialEntries: ["/home"] }
  );

  return render(
    <AppProvider>
      <RouterProvider router={router} />
    </AppProvider>
  );
}

/** Open the command palette and wait for it to render */
async function openPalette() {
  fireEvent.keyDown(document, { key: "k", metaKey: true });
  await waitFor(() => {
    expect(screen.getByPlaceholderText("Search pages, tools, molecules...")).toBeInTheDocument();
  });
}

describe("CommandPalette", () => {
  beforeEach(() => {
    navigateMock.mockClear();
  });

  it("opens when Cmd+K is pressed", async () => {
    renderPalette();
    // Palette should not be visible initially
    expect(screen.queryByPlaceholderText("Search pages, tools, molecules...")).toBeNull();
    await openPalette();
  });

  it("opens when Ctrl+K is pressed", async () => {
    renderPalette();
    fireEvent.keyDown(document, { key: "k", ctrlKey: true });
    await waitFor(() => {
      expect(screen.getByPlaceholderText("Search pages, tools, molecules...")).toBeInTheDocument();
    });
  });

  it("shows Pages, Tools, and Example Molecules group headings", async () => {
    renderPalette();
    await openPalette();

    // Group headings have cmdk-group-heading attribute
    const headings = document.querySelectorAll("[cmdk-group-heading]");
    const headingTexts = Array.from(headings).map((h) => h.textContent);
    expect(headingTexts).toContain("Pages");
    expect(headingTexts).toContain("Tools");
    expect(headingTexts).toContain("Example Molecules");
  });

  it("shows page entries", async () => {
    renderPalette();
    await openPalette();
    expect(screen.getByText("Home")).toBeInTheDocument();
    // "Chemical Analysis" appears as a page entry and as section labels on tools,
    // so use getAllByText and verify at least one exists
    expect(screen.getAllByText("Chemical Analysis").length).toBeGreaterThanOrEqual(1);
    expect(screen.getByText("OCSR")).toBeInTheDocument();
  });

  it("shows tool entries", async () => {
    renderPalette();
    await openPalette();
    expect(screen.getByText("Descriptors")).toBeInTheDocument();
    expect(screen.getByText("Structure Explorer")).toBeInTheDocument();
  });

  it("shows example molecule entries", async () => {
    renderPalette();
    await openPalette();
    expect(screen.getByText("Caffeine")).toBeInTheDocument();
    expect(screen.getByText("Aspirin")).toBeInTheDocument();
  });
});
