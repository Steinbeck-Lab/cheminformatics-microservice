import { describe, it, expect, vi } from "vitest";
import { fireEvent, render, screen, within } from "@testing-library/react";
import { AppProvider } from "@/context/AppContext";
import SMILESInput from "@/components/common/SMILESInput";

function renderInput(value = "", onChange = vi.fn()) {
  render(
    <AppProvider>
      <SMILESInput value={value} onChange={onChange} required />
    </AppProvider>
  );
  return onChange;
}

describe("SMILESInput", () => {
  it("is a labelled monospace field", () => {
    const onChange = renderInput("CCO");
    const field = screen.getByLabelText(/SMILES/);
    expect(field).toHaveValue("CCO");
    expect(field).toHaveClass("pc-well");
    fireEvent.change(field, { target: { value: "CCN" } });
    expect(onChange).toHaveBeenCalledWith("CCN");
  });

  it("shows the examples as chips that fill the field", () => {
    const onChange = renderInput();
    const examples = screen.getByRole("group", { name: "Examples" });
    fireEvent.click(within(examples).getByRole("button", { name: "Caffeine" }));
    expect(onChange).toHaveBeenCalledWith("CN1C=NC2=C1C(=O)N(C(=O)N2C)C");
  });

  it("marks the example that is in the field", () => {
    renderInput("CC(=O)OC1=CC=CC=C1C(=O)O");
    expect(screen.getByRole("button", { name: "Aspirin" })).toHaveAttribute("aria-pressed", "true");
  });

  it("clears the field", () => {
    const onChange = renderInput("CCO");
    fireEvent.click(screen.getByRole("button", { name: "Clear input" }));
    expect(onChange).toHaveBeenCalledWith("");
  });
});
