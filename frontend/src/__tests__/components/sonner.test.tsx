import { describe, it, expect } from "vitest";
import { render } from "@testing-library/react";
import { Toaster } from "@/components/ui/sonner";
import { AppProvider } from "@/context/AppContext";

describe("Sonner Toaster", () => {
  it("renders without crashing", () => {
    const { container } = render(
      <AppProvider>
        <Toaster />
      </AppProvider>
    );
    expect(container).toBeTruthy();
  });
});
