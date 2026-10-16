import { readFileSync } from "node:fs";
import { resolve } from "node:path";
import { describe, it, expect } from "vitest";

const css = readFileSync(resolve(__dirname, "../../styles/tailwind.css"), "utf8");

// Body of one `@utility name { … }` block (brace-balanced).
function utility(name: string): string {
  const start = css.indexOf(`@utility ${name} {`);
  if (start < 0) throw new Error(`no @utility ${name}`);
  let depth = 0;
  for (let i = css.indexOf("{", start); i < css.length; i++) {
    if (css[i] === "{") depth++;
    if (css[i] === "}" && --depth === 0) return css.slice(start, i + 1);
  }
  throw new Error(`unclosed @utility ${name}`);
}

// Pages still pass their own hover:bg-*, bg-gradient and text colours to <Button>.
// Clay must not out-rank them: zero-specificity states, and no `background` shorthand
// (it would wipe a caller's background-image).
describe.each(["clay", "clay-primary", "clay-outline", "clay-danger"])("@utility %s", (name) => {
  const body = utility(name);

  it("never uses the background shorthand", () => {
    expect(body).not.toMatch(/\bbackground\s*:/);
  });

  it("keeps hover, active and disabled states at zero specificity", () => {
    expect(body).not.toMatch(/&:(hover|active|disabled)/);
  });
});
