import { fireEvent, render, screen, within } from "@testing-library/react";
import { describe, it, expect } from "vitest";
import { createMemoryRouter, RouterProvider } from "react-router-dom";
import { SectionPage } from "@/components/common/SectionPage";

const tools = [
  { id: "one", name: "Tool One", description: "First tool", component: () => <p>one body</p> },
  { id: "two", name: "Tool Two", description: "Second tool", component: () => <p>two body</p> },
];

function renderAt(path: string) {
  const page = <SectionPage label="Demo tools" basePath="/demo" tools={tools} />;
  const router = createMemoryRouter(
    [
      { path: "/demo", element: page },
      { path: "/demo/:toolId", element: page },
    ],
    { initialEntries: [path] }
  );
  render(<RouterProvider router={router} />);
  return router;
}

describe("SectionPage", () => {
  it("shows a tab per tool and marks the open one", () => {
    renderAt("/demo/two");
    const tabs = within(screen.getByRole("navigation", { name: "Demo tools" })).getAllByRole(
      "link"
    );
    expect(tabs.map((t) => t.textContent)).toEqual(["Tool One", "Tool Two"]);
    expect(tabs[1]).toHaveAttribute("aria-current", "page");
    expect(tabs[1]).toHaveAttribute("href", "/demo/two");
  });

  it("titles the page with the open tool and renders it", () => {
    renderAt("/demo/two");
    expect(screen.getByRole("heading", { level: 1, name: "Tool Two" })).toBeInTheDocument();
    expect(screen.getByText("Second tool")).toBeInTheDocument();
    expect(screen.getByText("two body")).toBeInTheDocument();
  });

  it.each(["/demo", "/demo/nope"])("sends %s to the first tool", async (path) => {
    const router = renderAt(path);
    expect(await screen.findByText("one body")).toBeInTheDocument();
    expect(router.state.location.pathname).toBe("/demo/one");
  });
});

describe("SectionPage with a sidebar", () => {
  const grouped = [
    { ...tools[0]!, category: "Structure", icon: () => <svg /> },
    { ...tools[1]!, category: "Analysis", icon: () => <svg /> },
  ];

  function renderSidebar(path = "/demo/two") {
    const page = (
      <SectionPage label="Demo tools" basePath="/demo" tools={grouped} layout="sidebar" />
    );
    const router = createMemoryRouter([{ path: "/demo/:toolId", element: page }], {
      initialEntries: [path],
    });
    render(<RouterProvider router={router} />);
  }

  it("groups tools under category headings by default", () => {
    renderSidebar();
    const nav = screen.getByRole("navigation", { name: "Demo tools" });
    expect(within(nav).getByRole("heading", { name: "Structure" })).toBeInTheDocument();
    expect(within(nav).getByRole("link", { name: "Tool Two" })).toHaveAttribute(
      "aria-current",
      "page"
    );
  });

  it("switches to one plain list without headings", () => {
    renderSidebar();
    fireEvent.click(screen.getByRole("radio", { name: "All tools" }));
    expect(screen.getByRole("radio", { name: "All tools" })).toHaveAttribute(
      "aria-checked",
      "true"
    );
    const nav = screen.getByRole("navigation", { name: "Demo tools" });
    expect(within(nav).queryByRole("heading")).toBeNull();
    expect(within(nav).getAllByRole("link")).toHaveLength(2);
  });

  it("has a show/hide button for the list on small screens", () => {
    renderSidebar();
    const btn = screen.getByRole("button", { name: "Show tools" });
    expect(btn).toHaveAttribute("aria-expanded", "false");
    fireEvent.click(btn);
    expect(screen.getByRole("button", { name: "Hide tools" })).toHaveAttribute(
      "aria-expanded",
      "true"
    );
  });
});
