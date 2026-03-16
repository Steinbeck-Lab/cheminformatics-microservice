// Description: This file contains the context provider for the application, managing global state such as theme, recent molecules, API settings, and loading/error states.
import { createContext, useContext, useState, useEffect, type ReactNode } from "react";
import type { RecentMolecule, ApiConfig, AppContextValue } from "../types/molecule";

// Create the context
const AppContext = createContext<AppContextValue | undefined>(undefined);

// Custom hook for using the context
export const useAppContext = (): AppContextValue => {
  const context = useContext(AppContext);
  if (!context) {
    throw new Error("useAppContext must be used within an AppProvider");
  }
  return context;
};

// Provider component
export const AppProvider = ({ children }: { children: ReactNode }) => {
  // Theme state (dark/light) -- default false; FOUC script handles initial render,
  // and the mount useEffect below will set the correct value from localStorage or system preference.
  const [isDarkMode, setIsDarkMode] = useState(false);

  // History of recent molecules
  const [recentMolecules, setRecentMolecules] = useState<RecentMolecule[]>([]);

  // API settings
  const [apiConfig, setApiConfig] = useState<ApiConfig>({
    baseUrl: import.meta.env.VITE_API_URL || "https://dev.api.naturalproducts.net/latest",
    timeout: 30000,
  });

  // Global loading state
  const [isLoading, setIsLoading] = useState(false);

  // Global error state
  const [globalError, setGlobalError] = useState<string | null>(null);

  // Toggle dark/light mode
  const toggleDarkMode = () => {
    setIsDarkMode((prev) => !prev);
  };

  // Add a molecule to recent history
  const addRecentMolecule = (molecule: Omit<RecentMolecule, "_savedAt">) => {
    if (!molecule.smiles) return;

    setRecentMolecules((prev) => {
      // Remove duplicate if exists
      const filtered = prev.filter((m) => m.smiles !== molecule.smiles);

      // Add timestamp and place at the beginning, limit to 10 items
      return [{ ...molecule, _savedAt: Date.now() }, ...filtered].slice(0, 10);
    });
  };

  // Clear recent molecules history
  const clearRecentMolecules = () => {
    setRecentMolecules([]);
  };

  // Update API configuration
  const updateApiConfig = (newConfig: Partial<ApiConfig>) => {
    setApiConfig((prev) => ({ ...prev, ...newConfig }));
  };

  // Initialize theme from localStorage on mount, with system preference fallback
  useEffect(() => {
    const savedDarkMode = localStorage.getItem("darkMode");
    if (savedDarkMode !== null) {
      setIsDarkMode(savedDarkMode === "true");
    } else {
      // No saved preference -- detect system preference
      const prefersDark = window.matchMedia("(prefers-color-scheme: dark)").matches;
      setIsDarkMode(prefersDark);
    }

    // Load recent molecules from localStorage, filtering out entries older than 24h
    const savedMolecules = localStorage.getItem("recentMolecules");
    if (savedMolecules) {
      try {
        const parsed = JSON.parse(savedMolecules) as RecentMolecule[];
        const cutoff = Date.now() - 24 * 60 * 60 * 1000;
        const fresh = Array.isArray(parsed)
          ? parsed.filter((m) => m._savedAt && m._savedAt > cutoff)
          : [];
        setRecentMolecules(fresh);
      } catch (error) {
        console.error("Failed to parse saved molecules", error);
      }
    }
  }, []);

  // Save theme preference to localStorage when it changes
  useEffect(() => {
    localStorage.setItem("darkMode", String(isDarkMode));

    // Update document class for global styling
    if (isDarkMode) {
      document.documentElement.classList.add("dark");
    } else {
      document.documentElement.classList.remove("dark");
    }
  }, [isDarkMode]);

  // Save recent molecules to localStorage when they change
  useEffect(() => {
    localStorage.setItem("recentMolecules", JSON.stringify(recentMolecules));
  }, [recentMolecules]);

  // Context value with all state and functions
  const contextValue: AppContextValue = {
    isDarkMode,
    toggleDarkMode,
    recentMolecules,
    addRecentMolecule,
    clearRecentMolecules,
    apiConfig,
    updateApiConfig,
    isLoading,
    setIsLoading,
    globalError,
    setGlobalError,
  };

  return <AppContext.Provider value={contextValue}>{children}</AppContext.Provider>;
};

export default AppContext;
