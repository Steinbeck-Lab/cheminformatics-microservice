import { describe, it, expect } from "vitest";
import { render, screen, within } from "@testing-library/react";
import { MemoryRouter } from "react-router-dom";
import AboutPage from "@/pages/AboutPage";

function renderAbout() {
  return render(
    <MemoryRouter>
      <AboutPage />
    </MemoryRouter>
  );
}

describe("About page", () => {
  it("names the service and states it is free, open source and needs no login", () => {
    renderAbout();
    expect(
      screen.getByRole("heading", { level: 1, name: /Cheminformatics Microservice/ })
    ).toBeInTheDocument();
    const facts = screen.getByRole("list", { name: "Terms of use" });
    expect(facts).toHaveTextContent(/MIT/);
    expect(facts).toHaveTextContent(/No login/);
  });

  it("shows the stack from browser to funders, top to bottom", () => {
    renderAbout();
    const layers = within(screen.getByRole("list", { name: "The service, layer by layer" }))
      .getAllByRole("heading", { level: 2 })
      .map((h) => h.textContent);
    expect(layers).toEqual([
      "In your browser",
      "One REST API",
      "Open toolkits",
      "A public instance",
      "People and funding",
    ]);
  });

  it("links the browser layer to every section", () => {
    renderAbout();
    for (const [name, href] of [
      ["Analysis", "/chem"],
      ["Convert", "/convert"],
      ["Depict", "/depict"],
      ["Tools", "/tools"],
      ["OCSR", "/ocsr"],
    ]) {
      expect(screen.getByRole("link", { name: new RegExp(`^${name}`) })).toHaveAttribute(
        "href",
        href
      );
    }
  });

  it("lists the papers and the software DOI to cite, each copyable", () => {
    renderAbout();
    const refs = screen.getByRole("region", { name: "How to cite" });
    for (const doi of [
      "https://doi.org/10.1186/s13321-025-01094-1",
      "https://doi.org/10.1186/s13321-023-00762-4",
      "https://doi.org/10.5281/zenodo.7745987",
    ]) {
      expect(within(refs).getByRole("link", { name: doi })).toHaveAttribute("href", doi);
    }
    expect(within(refs).getAllByRole("button", { name: /^Copy citation/ })).toHaveLength(3);
  });

  it("credits the lab, the university and the DFG projects", () => {
    renderAbout();
    expect(screen.getAllByRole("link", { name: /Steinbeck/ }).length).toBeGreaterThan(0);
    expect(
      screen.getAllByRole("link", { name: /Friedrich Schiller University Jena/ }).length
    ).toBeGreaterThan(0);
    expect(screen.getByText(/441958208/)).toBeInTheDocument();
    expect(screen.getByText(/239748522/)).toBeInTheDocument();
  });

  it("does not offer IUPAC naming (STOUT is unavailable)", () => {
    renderAbout();
    expect(document.body).not.toHaveTextContent(/IUPAC name/i);
    expect(document.body).not.toHaveTextContent(/STOUT/);
  });
});
