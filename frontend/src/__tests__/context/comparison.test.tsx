import { describe, it, expect } from "vitest";
import { renderHook, act } from "@testing-library/react";
import React from "react";
import { ComparisonProvider, useComparison } from "../../context/ComparisonContext";

const wrapper = ({ children }: { children: React.ReactNode }) => (
  <ComparisonProvider>{children}</ComparisonProvider>
);

const molA = { smiles: "CCO", title: "Ethanol" };
const molB = { smiles: "CC(=O)O", title: "Acetic acid" };
const molC = { smiles: "C1=CC=CC=C1", title: "Benzene" };

describe("ComparisonContext", () => {
  it("starts with empty molecules array", () => {
    const { result } = renderHook(() => useComparison(), { wrapper });
    expect(result.current.molecules).toEqual([]);
    expect(result.current.canAdd).toBe(true);
    expect(result.current.isOpen).toBe(false);
    expect(result.current.isComparing).toBe(false);
  });

  it("addMolecule adds a molecule", () => {
    const { result } = renderHook(() => useComparison(), { wrapper });
    act(() => result.current.addMolecule(molA));
    expect(result.current.molecules).toHaveLength(1);
    expect(result.current.molecules[0].smiles).toBe("CCO");
  });

  it("addMolecule rejects duplicate SMILES", () => {
    const { result } = renderHook(() => useComparison(), { wrapper });
    act(() => result.current.addMolecule(molA));
    act(() => result.current.addMolecule({ ...molA, title: "Duplicate" }));
    expect(result.current.molecules).toHaveLength(1);
  });

  it("addMolecule rejects when already at 2 molecules", () => {
    const { result } = renderHook(() => useComparison(), { wrapper });
    act(() => result.current.addMolecule(molA));
    act(() => result.current.addMolecule(molB));
    act(() => result.current.addMolecule(molC));
    expect(result.current.molecules).toHaveLength(2);
  });

  it("canAdd returns false when 2 molecules present", () => {
    const { result } = renderHook(() => useComparison(), { wrapper });
    act(() => result.current.addMolecule(molA));
    act(() => result.current.addMolecule(molB));
    expect(result.current.canAdd).toBe(false);
  });

  it("removeMolecule removes by SMILES", () => {
    const { result } = renderHook(() => useComparison(), { wrapper });
    act(() => result.current.addMolecule(molA));
    act(() => result.current.addMolecule(molB));
    act(() => result.current.removeMolecule("CCO"));
    expect(result.current.molecules).toHaveLength(1);
    expect(result.current.molecules[0].smiles).toBe("CC(=O)O");
  });

  it("clearAll resets to empty state", () => {
    const { result } = renderHook(() => useComparison(), { wrapper });
    act(() => result.current.addMolecule(molA));
    act(() => result.current.addMolecule(molB));
    act(() => result.current.clearAll());
    expect(result.current.molecules).toEqual([]);
    expect(result.current.isOpen).toBe(false);
    expect(result.current.isComparing).toBe(false);
  });

  it("auto-opens tray on first molecule add", () => {
    const { result } = renderHook(() => useComparison(), { wrapper });
    expect(result.current.isOpen).toBe(false);
    act(() => result.current.addMolecule(molA));
    expect(result.current.isOpen).toBe(true);
  });

  it("closes tray when last molecule removed", () => {
    const { result } = renderHook(() => useComparison(), { wrapper });
    act(() => result.current.addMolecule(molA));
    expect(result.current.isOpen).toBe(true);
    act(() => result.current.removeMolecule("CCO"));
    expect(result.current.isOpen).toBe(false);
  });
});

describe("ComparisonTray", () => {
  it.todo("renders with glass-bold class");
  it.todo("shows molecule titles when molecules present");
  it.todo("disables Compare button when fewer than 2 molecules");
  it.todo("enables Compare button when 2 molecules present");
});
