import { describe, it, expect, beforeEach, vi } from "vitest";
import { render, act } from "@testing-library/react";
import { AppProvider, useAppContext } from "@/context/AppContext";

// Helper component to access context values
function ThemeConsumer({
  onRender,
}: {
  onRender: (ctx: ReturnType<typeof useAppContext>) => void;
}) {
  const ctx = useAppContext();
  onRender(ctx);
  return null;
}

describe("Theme context", () => {
  beforeEach(() => {
    localStorage.clear();
    document.documentElement.classList.remove("dark");
    vi.restoreAllMocks();
  });

  it.todo("THEME-01: toggleDarkMode switches isDarkMode state");

  it("THEME-02: persists theme preference to localStorage", () => {
    let contextValue: ReturnType<typeof useAppContext> | undefined;

    render(
      <AppProvider>
        <ThemeConsumer
          onRender={(ctx) => {
            contextValue = ctx;
          }}
        />
      </AppProvider>
    );

    // Toggle dark mode on
    act(() => {
      contextValue!.toggleDarkMode();
    });

    expect(localStorage.getItem("darkMode")).toBe("true");

    // Toggle dark mode off
    act(() => {
      contextValue!.toggleDarkMode();
    });

    expect(localStorage.getItem("darkMode")).toBe("false");
  });

  it("THEME-03: detects system preference when no localStorage value exists", () => {
    // Mock system preference as dark
    const matchMediaMock = vi.fn().mockImplementation((query: string) => ({
      matches: query === "(prefers-color-scheme: dark)",
      media: query,
      onchange: null,
      addListener: vi.fn(),
      removeListener: vi.fn(),
      addEventListener: vi.fn(),
      removeEventListener: vi.fn(),
      dispatchEvent: vi.fn(),
    }));
    vi.stubGlobal("matchMedia", matchMediaMock);

    let contextValue: ReturnType<typeof useAppContext> | undefined;

    render(
      <AppProvider>
        <ThemeConsumer
          onRender={(ctx) => {
            contextValue = ctx;
          }}
        />
      </AppProvider>
    );

    // Should detect system dark preference
    expect(contextValue!.isDarkMode).toBe(true);
    expect(matchMediaMock).toHaveBeenCalledWith("(prefers-color-scheme: dark)");
  });

  it("defaults to system preference for first-time visitors", () => {
    // Mock system preference as light
    const matchMediaMock = vi.fn().mockImplementation((query: string) => ({
      matches: false, // system prefers light
      media: query,
      onchange: null,
      addListener: vi.fn(),
      removeListener: vi.fn(),
      addEventListener: vi.fn(),
      removeEventListener: vi.fn(),
      dispatchEvent: vi.fn(),
    }));
    vi.stubGlobal("matchMedia", matchMediaMock);

    let contextValue: ReturnType<typeof useAppContext> | undefined;

    render(
      <AppProvider>
        <ThemeConsumer
          onRender={(ctx) => {
            contextValue = ctx;
          }}
        />
      </AppProvider>
    );

    // Should default to light mode when system prefers light
    expect(contextValue!.isDarkMode).toBe(false);
  });
});
