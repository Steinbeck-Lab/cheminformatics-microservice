/**
 * ComparisonContext -- manages molecule comparison state (max 2 molecules).
 *
 * Provides a persistent comparison tray (like a shopping cart drawer) that
 * lets users accumulate molecules and compare them side-by-side.
 *
 * Uses React Context (not zustand) since the state is limited to max 2 items.
 */
import React, { createContext, useContext, useState, useCallback } from "react";

/** A molecule added to the comparison tray. */
export interface ComparisonMolecule {
  smiles: string;
  title: string;
  imageUrl?: string;
  descriptors?: Record<string, string | number>;
  sourceToolId?: string;
}

/** Context value exposed by ComparisonProvider. */
export interface ComparisonContextValue {
  molecules: ComparisonMolecule[];
  isOpen: boolean;
  isComparing: boolean;
  canAdd: boolean;
  addMolecule: (mol: ComparisonMolecule) => void;
  removeMolecule: (smiles: string) => void;
  clearAll: () => void;
  openTray: () => void;
  closeTray: () => void;
  startCompare: () => void;
  stopCompare: () => void;
}

const ComparisonContext = createContext<ComparisonContextValue | undefined>(undefined);

/** Hook to access comparison state. Must be used within ComparisonProvider. */
export const useComparison = (): ComparisonContextValue => {
  const context = useContext(ComparisonContext);
  if (!context) {
    throw new Error("useComparison must be used within a ComparisonProvider");
  }
  return context;
};

const MAX_MOLECULES = 2;

/** Provider component for molecule comparison state. */
export const ComparisonProvider = ({ children }: { children: React.ReactNode }) => {
  const [molecules, setMolecules] = useState<ComparisonMolecule[]>([]);
  const [isOpen, setIsOpen] = useState(false);
  const [isComparing, setIsComparing] = useState(false);

  const canAdd = molecules.length < MAX_MOLECULES;

  const addMolecule = useCallback((mol: ComparisonMolecule) => {
    setMolecules((prev) => {
      // Reject if already at max
      if (prev.length >= MAX_MOLECULES) return prev;
      // Reject duplicate SMILES
      if (prev.some((m) => m.smiles === mol.smiles)) return prev;
      return [...prev, mol];
    });
    // Auto-open tray on first add
    setIsOpen(true);
  }, []);

  const removeMolecule = useCallback((smiles: string) => {
    setMolecules((prev) => {
      const next = prev.filter((m) => m.smiles !== smiles);
      // Close tray if last molecule removed
      if (next.length === 0) {
        setIsOpen(false);
        setIsComparing(false);
      }
      // Stop comparing if below 2
      if (next.length < MAX_MOLECULES) {
        setIsComparing(false);
      }
      return next;
    });
  }, []);

  const clearAll = useCallback(() => {
    setMolecules([]);
    setIsOpen(false);
    setIsComparing(false);
  }, []);

  const openTray = useCallback(() => setIsOpen(true), []);
  const closeTray = useCallback(() => setIsOpen(false), []);

  const startCompare = useCallback(() => {
    setIsComparing(true);
  }, []);

  const stopCompare = useCallback(() => {
    setIsComparing(false);
  }, []);

  const contextValue: ComparisonContextValue = {
    molecules,
    isOpen,
    isComparing,
    canAdd,
    addMolecule,
    removeMolecule,
    clearAll,
    openTray,
    closeTray,
    startCompare,
    stopCompare,
  };

  return <ComparisonContext.Provider value={contextValue}>{children}</ComparisonContext.Provider>;
};

export default ComparisonContext;
