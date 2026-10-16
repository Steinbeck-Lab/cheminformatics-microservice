import { describe, it, expect, vi } from "vitest";
import { render, screen } from "@testing-library/react";
import userEvent from "@testing-library/user-event";
import { Button } from "@/components/ui/button";

describe("Button", () => {
  it("default variant is indigo clay", () => {
    render(<Button>Depict</Button>);
    const b = screen.getByRole("button", { name: "Depict" });
    expect(b).toHaveClass("clay", "clay-primary");
  });

  it.each(["secondary", "outline"] as const)("%s variant is madder-outlined clay", (variant) => {
    render(<Button variant={variant}>Paste</Button>);
    const b = screen.getByRole("button", { name: "Paste" });
    expect(b).toHaveClass("clay", "clay-outline");
    expect(b).not.toHaveClass("clay-primary");
  });

  it("destructive variant is madder clay", () => {
    render(<Button variant="destructive">Clear</Button>);
    expect(screen.getByRole("button", { name: "Clear" })).toHaveClass("clay", "clay-danger");
  });

  it.each(["ghost", "link"] as const)("%s variant is not clay", (variant) => {
    render(<Button variant={variant}>More</Button>);
    expect(screen.getByRole("button", { name: "More" })).not.toHaveClass("clay");
  });

  it("keeps the touch target on every size", () => {
    render(<Button size="icon-sm" aria-label="Copy" />);
    expect(screen.getByRole("button", { name: "Copy" })).toHaveClass("touch-target");
  });

  it("forwards onClick", async () => {
    const onClick = vi.fn();
    render(<Button onClick={onClick}>Go</Button>);
    await userEvent.click(screen.getByRole("button", { name: "Go" }));
    expect(onClick).toHaveBeenCalledOnce();
  });
});

describe("Button disabled look", () => {
  it.each(["ghost", "link"] as const)("%s variant dims when disabled", (variant) => {
    render(
      <Button variant={variant} disabled>
        Off
      </Button>
    );
    expect(screen.getByRole("button", { name: "Off" })).toHaveClass("disabled:opacity-50");
  });
});
