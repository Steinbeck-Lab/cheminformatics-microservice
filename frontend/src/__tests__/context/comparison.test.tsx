import { describe, it, expect } from "vitest";

describe("ComparisonContext", () => {
  it.todo("starts with empty molecules array");
  it.todo("addMolecule adds a molecule");
  it.todo("addMolecule rejects duplicate SMILES");
  it.todo("addMolecule rejects when already at 2 molecules");
  it.todo("canAdd returns false when 2 molecules present");
  it.todo("removeMolecule removes by SMILES");
  it.todo("clearAll resets to empty state");
  it.todo("auto-opens tray on first molecule add");
  it.todo("closes tray when last molecule removed");
});

describe("ComparisonTray", () => {
  it.todo("renders with glass-bold class");
  it.todo("shows molecule titles when molecules present");
  it.todo("disables Compare button when fewer than 2 molecules");
  it.todo("enables Compare button when 2 molecules present");
});
