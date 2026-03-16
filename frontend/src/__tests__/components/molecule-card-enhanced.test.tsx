import { type PropsWithChildren } from "react";
import { describe, it, expect, vi } from "vitest";
import { render, screen } from "@testing-library/react";
import userEvent from "@testing-library/user-event";

// Mock AppContext
vi.mock("../../context/AppContext", () => ({
  useAppContext: () => ({
    addRecentMolecule: vi.fn(),
  }),
}));

// Mock motion/react to simplify animation testing
vi.mock("motion/react", () => ({
  motion: {
    div: ({ children, ...props }: PropsWithChildren<Record<string, unknown>>) => (
      <div data-testid="motion-div" {...props}>
        {children}
      </div>
    ),
  },
  AnimatePresence: ({ children }: PropsWithChildren) => <>{children}</>,
}));

// Lazy import after mocks
import MoleculeCard from "@/components/common/MoleculeCard";
import { ComparisonProvider } from "@/context/ComparisonContext";

function Wrapper({ children }: PropsWithChildren) {
  return <ComparisonProvider>{children}</ComparisonProvider>;
}

describe("MoleculeCard enhanced", () => {
  const defaultProps = {
    smiles: "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
    title: "Caffeine",
  };

  it("renders with glass-bold class", () => {
    const { container } = render(<MoleculeCard {...defaultProps} />, {
      wrapper: Wrapper,
    });
    const card = container.querySelector("[role='figure'], [role='button']");
    expect(card?.className).toContain("glass-bold");
  });

  it("renders with clay-interactive class", () => {
    const { container } = render(<MoleculeCard {...defaultProps} />, {
      wrapper: Wrapper,
    });
    const card = container.querySelector("[role='figure'], [role='button']");
    expect(card?.className).toContain("clay-interactive");
  });

  it("has hover translate and shadow classes", () => {
    const { container } = render(<MoleculeCard {...defaultProps} />, {
      wrapper: Wrapper,
    });
    const card = container.querySelector("[role='figure'], [role='button']");
    expect(card?.className).toContain("hover:-translate-y-1");
    expect(card?.className).toContain("hover:shadow-lg");
    expect(card?.className).toContain("hover:shadow-primary/20");
  });

  it("renders expand/collapse toggle when isExpandable is true", () => {
    render(<MoleculeCard {...defaultProps} isExpandable={true} />, {
      wrapper: Wrapper,
    });
    expect(screen.getByLabelText("Expand details")).toBeInTheDocument();
  });

  it("does not render expand toggle when isExpandable is false", () => {
    render(<MoleculeCard {...defaultProps} isExpandable={false} />, {
      wrapper: Wrapper,
    });
    expect(screen.queryByLabelText("Expand details")).not.toBeInTheDocument();
  });

  it("shows additional details when expanded", async () => {
    const user = userEvent.setup();
    render(<MoleculeCard {...defaultProps} isExpandable={true} />, {
      wrapper: Wrapper,
    });

    const toggle = screen.getByLabelText("Expand details");
    await user.click(toggle);

    expect(screen.getByLabelText("Collapse details")).toBeInTheDocument();
  });

  it("hides additional details when collapsed", () => {
    render(<MoleculeCard {...defaultProps} isExpandable={true} />, {
      wrapper: Wrapper,
    });
    expect(screen.queryByLabelText("Collapse details")).not.toBeInTheDocument();
  });
});
