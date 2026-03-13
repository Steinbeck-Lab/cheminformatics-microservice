# Phase 4: Component System + Dark Mode - Gap Closure Research

**Researched:** 2026-03-13
**Domain:** shadcn/ui Button asChild + motion/react animation for OCSRPage and HomePage gap closure
**Confidence:** HIGH

## Summary

Phase 4 is verified complete except for two pages (OCSRPage and HomePage) that still use hand-rolled button-like elements and card-like containers instead of shadcn/ui components. The gaps are styling inconsistencies, not functional failures. The verification flagged three specific issues: (1) OCSRPage `motion.a` external link buttons (lines 103-124), (2) OCSRPage hand-rolled card divs (lines 133-168), (3) HomePage primary CTA button (`<Link>` with inline gradient, line 228-237) and secondary external link buttons (lines 239-260), and (4) HomePage feature cards using `feature-card-enhanced` class with TiltCard wrapper.

This research focuses narrowly on the patterns needed to close these gaps: wrapping shadcn Button with motion for external links, using `asChild` with anchor tags, and applying shadcn Card to info-card and feature-card containers.

**Primary recommendation:** Use `motion.div` wrappers around shadcn `<Button asChild>` + `<a>` for external links (preserving motion animation variants on the wrapper div, not on Button itself). For cards, use shadcn `<Card>` + `<CardContent>` with className overrides to match existing visual styles. Do NOT use `motion.create(Button)` -- it has known issues with Radix Slot forwarding.

<user_constraints>
## User Constraints (from CONTEXT.md)

### Locked Decisions
- Install shadcn/ui components: Button, Card, Dialog, Sheet, Input, Select, Textarea, Tabs
- Migration approach: component-by-component across ALL 9 pages
- Wrap shadcn/ui components inside motion.div wrappers to preserve existing page transitions, stagger animations, and spring physics
- The shadcn component handles styling and accessibility; motion handles animation

### Claude's Discretion
- Any additional shadcn/ui components beyond the required set if clearly needed during migration

### Deferred Ideas (OUT OF SCOPE)
- Bento grid layouts for pages -- Phase 5 or Phase 6
- Toast notifications (Sonner) -- Phase 5
- Skeleton loading states -- Phase 5
</user_constraints>

<phase_requirements>
## Phase Requirements

| ID | Description | Research Support |
|----|-------------|-----------------|
| COMP-02 | Core buttons replaced with shadcn/ui Button component | OCSRPage motion.a buttons and HomePage CTA/secondary buttons must use shadcn Button |
| COMP-03 | Cards and content containers replaced with shadcn/ui Card | OCSRPage info cards and HomePage feature cards must use shadcn Card |
| COMP-07 | Consistent component usage across all 9 pages | OCSRPage and HomePage are the two remaining non-consistent pages |
</phase_requirements>

## Standard Stack

### Core (Already Installed)
| Library | Version | Purpose | Relevant API |
|---------|---------|---------|--------------|
| motion | ^12.36.0 | Animation | `motion.div` wrapper, `variants`, `whileHover` |
| radix-ui | ^1.4.3 | Primitive | `Slot.Root` for `asChild` prop forwarding |
| class-variance-authority | ^0.7.1 | Variant styling | `buttonVariants` for className generation |
| @/components/ui/button | local | shadcn Button | `asChild`, `variant`, `size`, `className` |
| @/components/ui/card | local | shadcn Card | `Card`, `CardContent`, `CardHeader`, `CardTitle` |

### No New Dependencies Needed
All required libraries are already installed. This is purely a component migration within existing pages.

## Architecture Patterns

### Pattern 1: External Link Button with Motion Animation (OCSRPage)

**What:** Replace `motion.a` styled-as-button with `motion.div` wrapper around `Button asChild` + `<a>`
**When to use:** External links that need to look like buttons AND have motion animation variants
**Why this pattern:** The `asChild` prop delegates rendering to the child `<a>` tag via Radix Slot, preserving the semantic anchor element for accessibility and external navigation. The `motion.div` wrapper handles animation independently from the Button component.

**Current code (OCSRPage lines 103-113):**
```tsx
<motion.a
  href="https://decimer.ai"
  target="_blank"
  rel="noopener noreferrer"
  className="flex items-center justify-center px-6 py-4 bg-linear-to-r from-blue-600 to-blue-700 hover:from-blue-500 hover:to-blue-600 text-white font-semibold rounded-lg shadow-lg w-full md:w-auto"
  variants={buttonVariant}
  whileHover="hover"
>
  <span className="mr-2 text-lg">Visit DECIMER</span>
  <ExternalLink className="h-5 w-5" />
</motion.a>
```

**Target pattern:**
```tsx
<motion.div variants={buttonVariant} whileHover="hover">
  <Button
    asChild
    size="lg"
    className="px-6 py-4 h-auto bg-linear-to-r from-blue-600 to-blue-700 hover:from-blue-500 hover:to-blue-600 text-white font-semibold rounded-lg shadow-lg w-full md:w-auto text-lg"
  >
    <a href="https://decimer.ai" target="_blank" rel="noopener noreferrer">
      Visit DECIMER
      <ExternalLink className="h-5 w-5" />
    </a>
  </Button>
</motion.div>
```

**Key details:**
- `motion.div` carries the animation variants and `whileHover` -- NOT the Button or anchor
- `Button asChild` renders the `<a>` as the root element (via Radix Slot), merging Button's className with the anchor
- The gradient `className` override on Button preserves the existing visual appearance
- `size="lg"` + `h-auto` + explicit padding overrides shadcn's default height constraint
- The `<a>` gets proper `target="_blank"` and `rel="noopener noreferrer"` for external links

**WARNING: Do NOT use `motion.create(Button)` or `motion(Button)`** -- multiple community reports confirm that combining `motion.create()` with Radix Slot's `asChild` causes unstable animations. The wrapper pattern is the established safe approach used in this project (see Header.tsx lines 203-214 for SheetTrigger asChild + motion.div wrapper precedent).

### Pattern 2: Internal Link Button with Motion Animation (HomePage)

**What:** Replace styled `<Link>` with `Button asChild` wrapping React Router `<Link>`
**When to use:** Internal navigation links styled as primary CTA buttons

**Current code (HomePage lines 228-237):**
```tsx
<Link
  to="/depict"
  className="group primary-button w-full sm:w-auto relative inline-flex items-center justify-center px-7 py-3 sm:px-9 sm:py-4 text-base sm:text-lg font-semibold text-white bg-linear-to-r from-sky-600 to-cyan-500 dark:from-blue-600 dark:to-cyan-500 rounded-xl shadow-lg hover:shadow-xl hover:brightness-110 ..."
>
  ...
</Link>
```

**Target pattern:**
```tsx
<Button
  asChild
  size="lg"
  className="group w-full sm:w-auto relative px-7 py-3 sm:px-9 sm:py-4 h-auto text-base sm:text-lg font-semibold text-white bg-linear-to-r from-sky-600 to-cyan-500 dark:from-blue-600 dark:to-cyan-500 rounded-xl shadow-lg hover:shadow-xl hover:brightness-110 dark:hover:shadow-cyan-500/30 transition-all duration-300 ease-out transform hover:scale-[1.06] focus:outline-hidden focus:ring-4 ring-offset-2 ring-offset-slate-100 dark:ring-offset-gray-950 focus:ring-cyan-500/50 dark:focus:ring-cyan-400/50"
>
  <Link to="/depict">
    <span className="absolute inset-0 rounded-xl bg-linear-to-r from-white/10 to-transparent opacity-0 group-hover:opacity-100 dark:from-white/5 transition-opacity duration-300 blur-xs"></span>
    <span className="relative z-10 flex items-center">
      Get Started <ArrowRight className="ml-2.5 h-5 w-5 transition-transform duration-300 group-hover:translate-x-1.5" />
    </span>
  </Link>
</Button>
```

**Key details:**
- The CTA button is already inside a `motion.div` (heroTextVariants) so no additional motion wrapper needed
- `Button asChild` renders the React Router `<Link>` as the root element
- The extensive className override preserves the exact gradient + hover + focus styling
- `h-auto` overrides shadcn's fixed height to allow padding-based sizing

### Pattern 3: Secondary External Link Buttons (HomePage)

**What:** Replace styled `<a>` tags with `Button asChild variant="outline"` wrapping `<a>`

**Current code (HomePage lines 239-260):** Two `<a>` tags for "Guides" and "API Docs" with outline/glass styling.

**Target pattern:**
```tsx
<Button
  asChild
  variant="outline"
  size="lg"
  className="group w-full sm:w-auto relative px-7 py-3 sm:px-9 sm:py-4 h-auto text-base sm:text-lg font-semibold text-slate-700 dark:text-slate-200 bg-white/90 dark:bg-slate-800/80 backdrop-blur-md border-slate-300 dark:border-slate-700 hover:bg-white dark:hover:bg-slate-700/90 hover:border-slate-400 dark:hover:border-slate-500 rounded-xl shadow-md hover:shadow-lg dark:hover:shadow-slate-700/40 transition-all duration-300 ease-out transform hover:scale-[1.06] focus:outline-hidden focus:ring-4 focus:ring-slate-400/50 dark:focus:ring-slate-600/50 ring-offset-2 ring-offset-slate-100 dark:ring-offset-gray-950"
>
  <a href="https://docs.api.naturalproducts.net/introduction.html" target="_blank" rel="noopener noreferrer">
    <span className="absolute inset-0 rounded-xl bg-linear-to-t from-black/5 to-transparent opacity-0 group-hover:opacity-50 dark:from-white/5 dark:group-hover:opacity-100 transition-opacity duration-300"></span>
    <span className="relative z-10 flex items-center">
      <Code className="mr-2 h-5 w-5 text-green-600 dark:text-green-400" /> Guides
    </span>
  </a>
</Button>
```

**Key details:**
- `variant="outline"` provides the base outline styling; className override handles the specific glass/backdrop-blur aesthetic
- These are already inside the `motion.div` (heroTextVariants custom={3}) so no additional wrapper needed
- Same `asChild` + `<a>` pattern as Pattern 1 but without needing a separate motion wrapper

### Pattern 4: Info Cards with shadcn Card (OCSRPage)

**What:** Replace hand-rolled card divs with shadcn `Card` + `CardContent`
**When to use:** Static information cards that need consistent border/bg theming

**Current code (OCSRPage lines 133-148):**
```tsx
<motion.div
  className="bg-white dark:bg-slate-800/80 rounded-lg p-5 border border-slate-200 dark:border-slate-700/50 shadow-md"
  variants={headerItemVariants}
>
  <h3 className="text-lg font-semibold text-blue-600 dark:text-blue-400 mb-2">
    What is DECIMER?
  </h3>
  <p className="text-muted-foreground mb-3">...</p>
</motion.div>
```

**Target pattern:**
```tsx
<motion.div variants={headerItemVariants}>
  <Card className="bg-white dark:bg-slate-800/80 rounded-lg border-slate-200 dark:border-slate-700/50 shadow-md">
    <CardContent className="p-5">
      <h3 className="text-lg font-semibold text-blue-600 dark:text-blue-400 mb-2">
        What is DECIMER?
      </h3>
      <p className="text-muted-foreground mb-3">...</p>
    </CardContent>
  </Card>
</motion.div>
```

**Key details:**
- `motion.div` wraps `Card` (not the other way around) -- motion handles animation, Card handles structure/theming
- Card's default gap-6 and py-6 must be overridden: use `CardContent className="p-5"` to match current spacing
- Card's default `bg-card` is overridden by explicit `bg-white dark:bg-slate-800/80` className to preserve the current visual
- Card's default `rounded-xl` is overridden to `rounded-lg` to match existing style
- The same approach applies to the sidebar info box (line 210) and MARCUS card (line 149)

### Pattern 5: Feature Cards with shadcn Card (HomePage)

**What:** Replace `feature-card-enhanced` + TiltCard wrapper with shadcn Card inside TiltCard
**When to use:** Interactive feature cards with 3D tilt effect and hover border reveal

**Current code (HomePage lines 285-317):**
```tsx
<TiltCard tiltIntensity={7}>
  <Link
    to={feature.link}
    className="feature-card-enhanced group relative flex flex-col h-full min-h-[270px] ..."
  >
    ...
  </Link>
</TiltCard>
```

**Target pattern:**
```tsx
<TiltCard tiltIntensity={7}>
  <Card className="group relative flex flex-col h-full min-h-[270px] sm:min-h-[290px] bg-white/80 dark:bg-slate-800/75 backdrop-blur-xl border-transparent rounded-3xl shadow-lg dark:shadow-xl hover:shadow-xl dark:hover:shadow-blue-900/30 overflow-hidden transform-style-3d p-0 gap-0 transition-all duration-400 ease-out">
    <Link to={feature.link} className="flex flex-col h-full p-4 sm:p-6" style={{ transform: "translateZ(0)", willChange: "border-color, backdrop-filter" }}>
      <div className="absolute inset-0 rounded-3xl bg-gradient-radial from-(--primary-accent-faint) via-transparent to-transparent opacity-0 group-hover:opacity-60 dark:group-hover:opacity-100 transition-opacity duration-500 pointer-events-none"></div>
      {/* ... icon, title, description, arrow ... */}
    </Link>
  </Card>
</TiltCard>
```

**Key details:**
- Card wraps the `<Link>` (Card is the visual container, Link is the interaction layer)
- `p-0 gap-0` on Card resets default padding/gap since Link provides its own `p-4 sm:p-6`
- The `feature-card-enhanced` CSS class can remain in the `<style>` block for the `::before` pseudo-element border reveal effect, applied to Card via additional className
- TiltCard wraps Card just as it wraps the current div -- no change to the animation layer
- The `border-transparent` + hover color change is critical for the border reveal effect
- Consider keeping `feature-card-enhanced` class on Card for the `::before` pseudo-element hover effect, OR migrate the effect to Tailwind utilities if possible

### Anti-Patterns to Avoid

- **Do NOT use `motion.create(Button)`:** Conflicts with Radix Slot's `asChild` prop forwarding. Animations become unstable. Use `motion.div` wrapper pattern instead.
- **Do NOT remove motion animation variants:** The `buttonVariant` (whileHover scale) and `headerItemVariants` (spring entrance) are part of the page's animation system. Preserve them on wrapper divs.
- **Do NOT restructure the TiltCard component:** TiltCard manages its own `motion.div` with perspective transforms. Keep it as-is; just change what goes inside it.
- **Do NOT use Card without className overrides:** shadcn Card's defaults (`bg-card`, `rounded-xl`, `gap-6`, `py-6`) differ from the existing visual. Always override to match.

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Button-styled links | Custom `<a>` with inline Tailwind button classes | `Button asChild` + `<a>` | Consistent variant system, focus ring, disabled states, accessibility attributes |
| Card containers | `<div>` with bg/border/shadow classes | `Card` + `CardContent` with className overrides | data-slot attributes, consistent theming, semantic structure |
| Button-styled React Router links | Custom `<Link>` with inline classes | `Button asChild` + `<Link>` | Preserves React Router navigation + shadcn Button semantics |

## Common Pitfalls

### Pitfall 1: Radix Slot Merging with className
**What goes wrong:** When `asChild` merges Button's className with the child's className, Tailwind specificity conflicts can occur (e.g., Button's `rounded-md` vs your `rounded-xl`).
**Why it happens:** Radix Slot concatenates classNames; last-wins in CSS but not in Tailwind's atomic classes.
**How to avoid:** Put ALL styling overrides on Button's `className` prop, not on the child element. The child should only have structural props (href, to, target, rel).
**Warning signs:** Border radius, padding, or background colors not matching expected values.

### Pitfall 2: shadcn Button Default Height Constraint
**What goes wrong:** Button has a fixed `h-9` (default), `h-8` (sm), or `h-10` (lg) that clips content styled with padding-based sizing.
**Why it happens:** The existing buttons use `py-3`/`py-4` padding for height; shadcn uses fixed height values.
**How to avoid:** Always add `h-auto` in className override when using custom padding instead of shadcn's fixed heights.
**Warning signs:** Text or icons cut off, button appears compressed vertically.

### Pitfall 3: Card Default Padding/Gap Conflicts
**What goes wrong:** Card has default `py-6` and `gap-6` that add unwanted spacing.
**Why it happens:** shadcn Card is designed for structured content (CardHeader + CardContent + CardFooter); direct children without these sub-components get unexpected spacing.
**How to avoid:** When using Card as a thin wrapper, add `p-0 gap-0` to Card and handle all padding in CardContent or the direct child.
**Warning signs:** Extra whitespace above/below card content.

### Pitfall 4: Feature Card Pseudo-element Border Effect
**What goes wrong:** The `feature-card-enhanced::before` CSS pseudo-element creates a border reveal effect on hover. If Card replaces the element that has this class, the effect breaks.
**Why it happens:** CSS pseudo-elements are tied to the element they're applied to. Card renders a `<div>` just like the current implementation, but class application must be explicit.
**How to avoid:** Apply `feature-card-enhanced` as an additional className on Card. The `::before` pseudo-element will still work since Card renders a `<div>`.
**Warning signs:** Missing gradient border effect on hover.

### Pitfall 5: motion.div Wrapper Disrupting Layout
**What goes wrong:** Adding `motion.div` as a wrapper around Card changes the flexbox/grid child behavior (the wrapper becomes the grid item, not Card).
**Why it happens:** An additional DOM element is inserted between the grid container and the Card.
**How to avoid:** For the OCSRPage info cards, the `motion.div` already exists as the grid item (it's the element with `variants={headerItemVariants}`). Just replace the inner div with Card -- no new wrapper needed.
**Warning signs:** Cards lose their grid sizing or alignment.

## Code Examples

### Complete OCSRPage DECIMER Button Migration
```tsx
// Before (lines 103-113)
<motion.a
  href="https://decimer.ai"
  target="_blank"
  rel="noopener noreferrer"
  className="flex items-center justify-center px-6 py-4 bg-linear-to-r from-blue-600 to-blue-700 hover:from-blue-500 hover:to-blue-600 text-white font-semibold rounded-lg shadow-lg w-full md:w-auto"
  variants={buttonVariant}
  whileHover="hover"
>
  <span className="mr-2 text-lg">Visit DECIMER</span>
  <ExternalLink className="h-5 w-5" />
</motion.a>

// After
<motion.div variants={buttonVariant} whileHover="hover">
  <Button
    asChild
    size="lg"
    className="h-auto px-6 py-4 bg-linear-to-r from-blue-600 to-blue-700 hover:from-blue-500 hover:to-blue-600 text-white font-semibold rounded-lg shadow-lg w-full md:w-auto text-lg"
  >
    <a href="https://decimer.ai" target="_blank" rel="noopener noreferrer">
      Visit DECIMER
      <ExternalLink className="h-5 w-5" />
    </a>
  </Button>
</motion.div>
```

### Complete OCSRPage Info Card Migration
```tsx
// Before (lines 133-148)
<motion.div
  className="bg-white dark:bg-slate-800/80 rounded-lg p-5 border border-slate-200 dark:border-slate-700/50 shadow-md"
  variants={headerItemVariants}
>
  <h3 className="text-lg font-semibold text-blue-600 dark:text-blue-400 mb-2">
    What is DECIMER?
  </h3>
  <p className="text-muted-foreground mb-3">...</p>
  <p className="text-muted-foreground">...</p>
</motion.div>

// After
<motion.div variants={headerItemVariants}>
  <Card className="bg-white dark:bg-slate-800/80 rounded-lg border-slate-200 dark:border-slate-700/50 shadow-md py-0 gap-0">
    <CardContent className="p-5">
      <h3 className="text-lg font-semibold text-blue-600 dark:text-blue-400 mb-2">
        What is DECIMER?
      </h3>
      <p className="text-muted-foreground mb-3">...</p>
      <p className="text-muted-foreground">...</p>
    </CardContent>
  </Card>
</motion.div>
```

### Import Statement for Both Pages
```tsx
// OCSRPage.tsx - add these imports
import { Button } from "@/components/ui/button";
import { Card, CardContent } from "@/components/ui/card";

// HomePage.tsx - add these imports
import { Button } from "@/components/ui/button";
import { Card } from "@/components/ui/card";
```

## Scope Assessment

### OCSRPage.tsx Changes
| Element | Lines | Change | Complexity |
|---------|-------|--------|------------|
| DECIMER motion.a button | 103-113 | motion.div wrapper + Button asChild + a | Low |
| MARCUS motion.a button | 114-124 | motion.div wrapper + Button asChild + a | Low |
| DECIMER info card div | 133-148 | Card + CardContent inside existing motion.div | Low |
| MARCUS info card div | 149-167 | Card + CardContent inside existing motion.div | Low |
| Sidebar info box div | 210-274 | Card + CardContent (same pattern) | Low |
| Add imports | top | Button, Card, CardContent | Trivial |

### HomePage.tsx Changes
| Element | Lines | Change | Complexity |
|---------|-------|--------|------------|
| Primary CTA Link | 228-237 | Button asChild + Link (already in motion container) | Low |
| Guides external a | 239-249 | Button asChild variant="outline" + a | Low |
| API Docs external a | 250-260 | Button asChild variant="outline" + a | Low |
| Feature cards (x5) | 285-317 | Card wrapping Link inside TiltCard | Medium |
| Recent molecule cards (x3) | 354-366 | Card wrapping MoleculeCard inside TiltCard | Medium |
| Add imports | top | Button, Card | Trivial |

### Pages NOT Flagged (No Changes Needed)
- ChemPage, DepictPage, ConvertPage, ToolsPage: Already use Button imports
- PrivacyPage, TermsPage: Static content pages, no buttons/cards
- AboutPage: Imports Button already; has `motion.a` instances but was NOT flagged in verification

### Note on AboutPage and Footer
AboutPage (lines 613, 658, 895) and Footer (line 350) also contain `motion.a` styled-as-button patterns identical to OCSRPage. These were NOT flagged in the verification report. The planner should decide whether to include them in scope for completeness (addressing the root pattern everywhere) or strictly limit to the two flagged pages. Recommendation: include them for true consistency per COMP-07 "all 9 pages."

## Validation Architecture

### Test Framework
| Property | Value |
|----------|-------|
| Framework | vitest 3.0.7 |
| Config file | frontend/vite.config.mts (test section) |
| Quick run command | `cd frontend && npx vitest run --reporter=verbose` |
| Full suite command | `cd frontend && npx vitest run --reporter=verbose` |

### Phase Requirements -> Test Map
| Req ID | Behavior | Test Type | Automated Command | File Exists? |
|--------|----------|-----------|-------------------|-------------|
| COMP-02 | OCSRPage imports Button from @/components/ui/button | unit (file content) | `cd frontend && npx vitest run src/__tests__/icon-migration.test.ts -t "OCSRPage" --reporter=verbose` | New test needed |
| COMP-02 | HomePage imports Button from @/components/ui/button | unit (file content) | `cd frontend && npx vitest run src/__tests__/icon-migration.test.ts -t "HomePage" --reporter=verbose` | New test needed |
| COMP-03 | OCSRPage imports Card from @/components/ui/card | unit (file content) | Same test file | New test needed |
| COMP-03 | HomePage imports Card from @/components/ui/card | unit (file content) | Same test file | New test needed |
| COMP-07 | No hand-rolled button-like elements in OCSRPage | unit (file content check: no motion.a with button styling) | Same test file | New test needed |
| COMP-07 | No feature-card-enhanced without Card wrapper in HomePage | manual-only | Visual verification | N/A |

### Sampling Rate
- **Per task commit:** `cd frontend && npx vitest run --reporter=verbose`
- **Per wave merge:** Same (single wave for gap closure)
- **Phase gate:** Full suite green + visual inspection of OCSRPage and HomePage in light/dark mode

### Wave 0 Gaps
- [ ] Add component-migration tests to existing `icon-migration.test.ts` (or create new `component-migration.test.ts`): verify OCSRPage imports Button and Card from @/components/ui/, verify HomePage imports Button and Card from @/components/ui/, verify no `motion.a` styled-as-button patterns remain in flagged pages

## Open Questions

1. **Should AboutPage and Footer motion.a patterns be included in scope?**
   - What we know: Verification only flagged OCSRPage and HomePage. But AboutPage has 3 instances of `motion.a` styled-as-button and Footer has 1 instance.
   - What's unclear: Whether the verifier intentionally excluded these (because AboutPage already imports Button for other elements) or missed them.
   - Recommendation: Include them in the gap closure plan for true COMP-07 consistency. The additional effort is minimal (same pattern, already understood).

2. **Should feature-card-enhanced CSS be removed or kept alongside Card?**
   - What we know: The `::before` pseudo-element on `.feature-card-enhanced` creates a gradient border hover effect. Card renders a `<div>` that can receive this class.
   - What's unclear: Whether the border reveal effect is worth preserving or should be simplified.
   - Recommendation: Keep the class on Card for visual continuity. Removing it is a design decision, not a component consistency issue.

## Sources

### Primary (HIGH confidence)
- Direct code inspection of `frontend/src/pages/OCSRPage.tsx` (303 lines)
- Direct code inspection of `frontend/src/pages/HomePage.tsx` (479 lines)
- Direct code inspection of `frontend/src/components/ui/button.tsx` (62 lines) -- confirmed `asChild` + Radix Slot pattern
- Direct code inspection of `frontend/src/components/ui/card.tsx` (75 lines) -- confirmed Card/CardContent structure
- Direct code inspection of `frontend/src/components/common/Header.tsx` (239 lines) -- confirmed motion.div + asChild pattern precedent at lines 203-214
- `frontend/package.json` -- confirmed motion ^12.36.0, radix-ui ^1.4.3

### Secondary (MEDIUM confidence)
- [shadcn/ui Button documentation](https://ui.shadcn.com/docs/components/radix/button) -- asChild usage with Link
- [shadcn/ui + Framer Motion discussion #1636](https://github.com/shadcn-ui/ui/discussions/1636) -- wrapper pattern recommended, motion.create + asChild unstable
- [Motion component docs](https://motion.dev/docs/react-motion-component) -- motion.create API reference

### Tertiary (LOW confidence)
- [shadcn/ui Issue #480](https://github.com/shadcn-ui/ui/issues/480) -- external link + asChild pattern discussions

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH -- all libraries already installed and in use across the project
- Architecture patterns: HIGH -- based on direct code inspection of working patterns in Header.tsx and other migrated pages
- Pitfalls: HIGH -- based on community reports and understanding of Radix Slot + className merging mechanics
- Card migration approach: HIGH -- Card component is simple (just a styled div), className overrides are straightforward

**Research date:** 2026-03-13
**Valid until:** 2026-04-13 (stable -- no dependency changes expected)
