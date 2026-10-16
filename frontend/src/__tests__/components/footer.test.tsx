import { render, screen, within } from "@testing-library/react";
import { MemoryRouter } from "react-router-dom";
import { describe, it, expect } from "vitest";
import { AppProvider } from "../../context/AppContext";
import Footer from "../../components/common/Footer";

function renderFooter() {
  return render(
    <AppProvider>
      <MemoryRouter>
        <Footer />
      </MemoryRouter>
    </AppProvider>
  );
}

describe("Footer", () => {
  it("spells the name as committed: Cheminformatics Microservice", () => {
    renderFooter();
    expect(screen.getByRole("contentinfo")).not.toHaveTextContent(/Microservices/);
  });

  it("links Terms and Privacy", () => {
    renderFooter();
    expect(screen.getByRole("link", { name: "Terms of Service" })).toHaveAttribute(
      "href",
      "/terms"
    );
    expect(screen.getByRole("link", { name: "Privacy Policy" })).toHaveAttribute(
      "href",
      "/privacy"
    );
  });

  it("has four resource cards that open in a new tab", () => {
    renderFooter();
    const links = within(screen.getByRole("navigation", { name: "Resources" })).getAllByRole(
      "link"
    );
    expect(links).toHaveLength(4);
    for (const a of links) {
      expect(a).toHaveAttribute("target", "_blank");
      expect(a).toHaveAttribute("rel", expect.stringContaining("noopener"));
      expect(a).toHaveTextContent("(opens in a new tab)");
    }
  });

  it("credits the lab, the university and the DFG, with no citation row", () => {
    renderFooter();
    const footer = screen.getByRole("contentinfo");
    expect(within(footer).getByRole("link", { name: "Steinbeck Lab" })).toBeInTheDocument();
    expect(
      within(footer).getByRole("link", { name: "Friedrich Schiller University Jena" })
    ).toBeInTheDocument();
    expect(footer).toHaveTextContent("Deutsche Forschungsgemeinschaft (DFG)");
    expect(footer).toHaveTextContent(`© ${new Date().getFullYear()}`);
    expect(footer).not.toHaveTextContent("J Cheminform");
  });
});
