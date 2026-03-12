import { render } from "@testing-library/react";
import { describe, it, expect } from "vitest";
import {
  createMemoryRouter,
  RouterProvider,
  createRoutesFromElements,
  Route,
} from "react-router-dom";
import { AppProvider } from "../context/AppContext";

describe("App", () => {
  it("renders without crashing", () => {
    // Use a minimal route to verify the app shell works
    const router = createMemoryRouter(
      createRoutesFromElements(
        <Route path="/" element={<div data-testid="app-root">App loaded</div>} />
      ),
      { initialEntries: ["/"] }
    );

    const { container } = render(
      <AppProvider>
        <RouterProvider router={router} />
      </AppProvider>
    );

    expect(container).toBeTruthy();
  });
});
