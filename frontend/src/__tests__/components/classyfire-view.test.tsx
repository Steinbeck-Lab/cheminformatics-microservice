import { afterEach, describe, it, expect, vi } from "vitest";
import { fireEvent, render, screen, within } from "@testing-library/react";
import { AppProvider } from "@/context/AppContext";
import ClassyfireView from "@/components/chem/ClassyfireView";

// The view talks to the ClassyFire endpoints with fetch; the service module is mocked so
// nothing else reaches the network.
vi.mock("@/services/chemService", () => ({}));

const RESULT = {
  classification_status: "Done",
  entities: [
    {
      smiles: "CCO",
      kingdom: { name: "Organic compounds", description: "Carbon-based compounds" },
      superclass: { name: "Organic oxygen compounds" },
      class_: { name: "Organooxygen compounds" },
      direct_parent: { name: "Primary alcohols" },
      molecular_framework: "Aliphatic acyclic compounds",
      substituents: ["Primary alcohol"],
      ancestors: ["Organic compounds", "Alcohols"],
    },
  ],
};

const json = (body: unknown, status = 200) =>
  new Response(JSON.stringify(body), { status, statusText: status === 200 ? "OK" : "Error" });

const renderView = () =>
  render(
    <AppProvider>
      <ClassyfireView />
    </AppProvider>
  );

describe("ClassyfireView (bento)", () => {
  afterEach(() => vi.unstubAllGlobals());

  it("lays out input, about and results panels", () => {
    renderView();
    const input = screen.getByRole("region", { name: "Structure input" });
    expect(within(input).getByRole("radiogroup", { name: "Mode" })).toBeInTheDocument();
    expect(screen.getByRole("region", { name: "About ClassyFire" })).toBeInTheDocument();
    expect(screen.getByRole("region", { name: "Results" })).toHaveTextContent(/Enter a SMILES/);
  });

  it("submits, polls and shows the classification", async () => {
    vi.stubGlobal(
      "fetch",
      vi.fn(async (url: string) => (url.includes("/classify") ? json({ id: "42" }) : json(RESULT)))
    );
    renderView();
    fireEvent.change(screen.getByLabelText(/SMILES/), { target: { value: "CCO" } });
    fireEvent.click(screen.getByRole("button", { name: /Classify molecule/i }));
    const results = await screen.findByRole("region", { name: "Classification" });
    expect(within(results).getByText("Primary alcohols")).toBeInTheDocument();
    expect(within(results).getByText("Primary alcohol")).toBeInTheDocument();
  });

  it("explains calmly when ClassyFire does not answer", async () => {
    vi.stubGlobal(
      "fetch",
      vi.fn(async () => json({}, 503))
    );
    renderView();
    fireEvent.change(screen.getByLabelText(/SMILES/), { target: { value: "CCO" } });
    fireEvent.click(screen.getByRole("button", { name: /Classify molecule/i }));
    const results = screen.getByRole("region", { name: "Results" });
    expect(await within(results).findByText(/ClassyFire did not answer/)).toBeInTheDocument();
  });
});
