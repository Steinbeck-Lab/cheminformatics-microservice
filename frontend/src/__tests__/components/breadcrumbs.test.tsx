import React from "react";
import { render, screen } from "@testing-library/react";
import { describe, it, expect } from "vitest";
import {
  createMemoryRouter,
  RouterProvider,
  createRoutesFromElements,
  Route,
} from "react-router-dom";
import { AppProvider } from "../../context/AppContext";
import { Breadcrumbs } from "../../components/common/Breadcrumbs";

function renderBreadcrumbs(path: string) {
  // Use route patterns that match what App.tsx defines so useParams works correctly
  const router = createMemoryRouter(
    createRoutesFromElements(
      <>
        <Route path="/chem" element={<Breadcrumbs />} />
        <Route path="/chem/:toolId" element={<Breadcrumbs />} />
        <Route path="/convert" element={<Breadcrumbs />} />
        <Route path="/convert/:convertId" element={<Breadcrumbs />} />
        <Route path="/depict" element={<Breadcrumbs />} />
        <Route path="/depict/:depictId" element={<Breadcrumbs />} />
        <Route path="/tools" element={<Breadcrumbs />} />
        <Route path="/tools/:toolId" element={<Breadcrumbs />} />
        <Route path="/home" element={<Breadcrumbs />} />
        <Route path="/about" element={<Breadcrumbs />} />
        <Route path="/ocsr" element={<Breadcrumbs />} />
      </>
    ),
    { initialEntries: [path] }
  );

  return render(
    <AppProvider>
      <RouterProvider router={router} />
    </AppProvider>
  );
}

describe("Breadcrumbs", () => {
  it("renders section name on /chem route", () => {
    renderBreadcrumbs("/chem");
    const nav = screen.getByTestId("breadcrumbs");
    expect(nav).toBeInTheDocument();
    expect(screen.getByText("Chemical Analysis")).toBeInTheDocument();
  });

  it("renders section + tool on /chem/descriptors", () => {
    renderBreadcrumbs("/chem/descriptors");
    expect(screen.getByText("Chemical Analysis")).toBeInTheDocument();
    expect(screen.getByText("Descriptors")).toBeInTheDocument();
  });

  it("does NOT render on /home", () => {
    renderBreadcrumbs("/home");
    expect(screen.queryByTestId("breadcrumbs")).toBeNull();
  });

  it("does NOT render on /about", () => {
    renderBreadcrumbs("/about");
    expect(screen.queryByTestId("breadcrumbs")).toBeNull();
  });

  it("does NOT render on /ocsr", () => {
    renderBreadcrumbs("/ocsr");
    expect(screen.queryByTestId("breadcrumbs")).toBeNull();
  });

  it("has hidden md:block class for desktop-only display", () => {
    renderBreadcrumbs("/chem");
    const nav = screen.getByTestId("breadcrumbs");
    expect(nav.className).toContain("hidden");
    expect(nav.className).toContain("md:block");
  });

  it("renders section + tool on /depict/structureexplorer", () => {
    renderBreadcrumbs("/depict/structureexplorer");
    expect(screen.getByText("Depiction")).toBeInTheDocument();
    expect(screen.getByText("Structure Explorer")).toBeInTheDocument();
  });

  it("renders section + tool on /tools/sugardetection", () => {
    renderBreadcrumbs("/tools/sugardetection");
    expect(screen.getByText("Tools")).toBeInTheDocument();
    expect(screen.getByText("Sugar Detection")).toBeInTheDocument();
  });
});
