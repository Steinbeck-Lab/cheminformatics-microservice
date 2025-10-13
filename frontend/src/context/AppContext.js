// Description: This file contains the context provider for the application, managing global state such as theme, recent molecules, API settings, and loading/error states.
import React, { createContext, useContext, useState, useEffect } from "react";

// Create the context
const AppContext = createContext();

// Custom hook for using the context
export const useAppContext = () => useContext(AppContext);

// Provider component
export const AppProvider = ({ children }) => {
  // Theme state (dark/light)
  const [isDarkMode, setIsDarkMode] = useState(true);

  // History of recent molecules
  const [recentMolecules, setRecentMolecules] = useState([]);

  // API settings
  const [apiConfig, setApiConfig] = useState({
    baseUrl: process.env.REACT_APP_API_URL || "https://dev.api.naturalproducts.net/latest",
    timeout: 30000,
  });

  // Global loading state
  const [isLoading, setIsLoading] = useState(false);

  // Global error state
  const [globalError, setGlobalError] = useState(null);

  // Toggle dark/light mode
  const toggleDarkMode = () => {
    setIsDarkMode((prev) => !prev);
  };

  // Add a molecule to recent history
  const addRecentMolecule = (molecule) => {
    if (!molecule.smiles) return;

    setRecentMolecules((prev) => {
      // Remove duplicate if exists
      const filtered = prev.filter((m) => m.smiles !== molecule.smiles);

      // Add to the beginning and limit to 10 items
      return [molecule, ...filtered].slice(0, 10);
    });
  };

  // Clear recent molecules history
  const clearRecentMolecules = () => {
    setRecentMolecules([]);
  };

  // Update API configuration
  const updateApiConfig = (newConfig) => {
    setApiConfig((prev) => ({ ...prev, ...newConfig }));
  };

  // Initialize theme from localStorage on mount
  useEffect(() => {
    const savedDarkMode = localStorage.getItem("darkMode");
    if (savedDarkMode !== null) {
      setIsDarkMode(savedDarkMode === "true");
    }

    // Load recent molecules from localStorage
    const savedMolecules = localStorage.getItem("recentMolecules");
    if (savedMolecules) {
      try {
        setRecentMolecules(JSON.parse(savedMolecules));
      } catch (error) {
        console.error("Failed to parse saved molecules", error);
      }
    }
  }, []);

  // Save theme preference to localStorage when it changes
  useEffect(() => {
    localStorage.setItem("darkMode", isDarkMode);

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
  const contextValue = {
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
