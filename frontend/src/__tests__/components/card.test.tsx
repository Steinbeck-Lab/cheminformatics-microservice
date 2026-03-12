import { describe, it, expect } from "vitest";
import { render, screen } from "@testing-library/react";
import { Card, CardHeader, CardContent, CardTitle } from "@/components/ui/card";

describe("shadcn/ui Card", () => {
  it("renders Card with CardHeader and CardContent", () => {
    render(
      <Card>
        <CardHeader>
          <CardTitle>Test Title</CardTitle>
        </CardHeader>
        <CardContent>Test Content</CardContent>
      </Card>
    );
    expect(screen.getByText("Test Title")).toBeInTheDocument();
    expect(screen.getByText("Test Content")).toBeInTheDocument();
  });

  it.todo("renders CardTitle and CardDescription");
  it.todo("applies custom className via cn()");
});
