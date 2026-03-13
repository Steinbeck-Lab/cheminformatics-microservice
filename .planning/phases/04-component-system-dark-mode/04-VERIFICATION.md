---
phase: 04-component-system-dark-mode
verified: 2026-03-13T08:48:22Z
status: gaps_found
score: 10/12 must-haves verified
re_verification: false
gaps:
  - truth: "All 9 pages use shadcn/ui components consistently -- no mix of hand-rolled and shadcn/ui buttons, cards, or form inputs"
    status: partial
    reason: "OCSRPage uses hand-rolled motion.a anchor elements styled as buttons and hand-rolled card divs instead of shadcn Button and Card components. HomePage uses a hand-rolled primary button (inline Tailwind gradient) and hand-rolled TiltCard / feature-card divs instead of shadcn Button and Card."
    artifacts:
      - path: "frontend/src/pages/OCSRPage.tsx"
        issue: "motion.a links styled as buttons (lines 103-124) and card divs (lines 133-168) are not shadcn Button or Card components. The page has no @/components/ui import."
      - path: "frontend/src/pages/HomePage.tsx"
        issue: "Primary button at line 230 is a styled <a> tag with inline gradient classes. Feature card containers are hand-rolled divs with className 'feature-card-enhanced'. No @/components/ui import in this page file."
    missing:
      - "OCSRPage: Replace motion.a external link buttons with shadcn Button variant='default' or variant='outline' wrapping an anchor, or use asChild with motion"
      - "OCSRPage: Replace hand-rolled info card divs (lines 133-168) with shadcn Card/CardContent"
      - "HomePage: Replace primary CTA button with shadcn Button (gradient via className override is acceptable)"
      - "HomePage: Replace TiltCard/feature-card-enhanced divs with shadcn Card where visually appropriate"
  - truth: "All buttons across all 9 pages use shadcn/ui Button component with appropriate variant"
    status: partial
    reason: "Same root cause as above: OCSRPage and HomePage still render hand-rolled button-like elements not backed by shadcn Button."
    artifacts:
      - path: "frontend/src/pages/OCSRPage.tsx"
        issue: "Two motion.a elements styled as buttons (DECIMER and MARCUS links at lines 103-124)"
      - path: "frontend/src/pages/HomePage.tsx"
        issue: "Primary CTA button at line 230 is a raw <a> element with inline Tailwind classes, not shadcn Button"
    missing:
      - "Migrate CTA and external link buttons to shadcn Button"
human_verification:
  - test: "Visual light/dark toggle on all 9 pages"
    expected: "Switching theme changes the entire UI between light and dark. No page shows incorrect colors or invisible text in either mode."
    why_human: "Visual correctness of OKLCH color palette, contrast ratios, and element visibility cannot be verified programmatically."
  - test: "Mobile Sheet navigation (375px viewport)"
    expected: "Hamburger button opens shadcn Sheet slide-out panel. Navigation links inside Sheet are clickable and close the Sheet. Mobile theme toggle inside Sheet works."
    why_human: "Sheet interaction, slide animation, and link closure behavior require browser rendering."
  - test: "FOUC prevention on hard refresh with dark mode active"
    expected: "Hard refresh (Cmd+Shift+R) with dark mode active shows no flash of light theme before dark mode applies."
    why_human: "Flash of Unstyled Content is a rendering timing issue observable only in the browser."
  - test: "Tab animations on ChemPage, DepictPage, ConvertPage, ToolsPage"
    expected: "Tab indicator animates smoothly when switching tabs. No broken motion wrappers."
    why_human: "Animation quality and motion wrapper correctness require visual inspection."
---

# Phase 4: Component System + Dark Mode Verification Report

**Phase Goal:** All UI elements use shadcn/ui components with a unified CSS variable theming system and functional dark mode toggle
**Verified:** 2026-03-13T08:48:22Z
**Status:** gaps_found
**Re-verification:** No — initial verification

## Goal Achievement

### Observable Truths (Success Criteria)

| # | Truth | Status | Evidence |
|---|-------|--------|----------|
| SC-1 | shadcn/ui installed with Button, Card, Dialog, Input, Select, Textarea, Sheet in components/ui/ | VERIFIED | `frontend/src/components/ui/` contains all 8 required files (button.tsx, card.tsx, dialog.tsx, input.tsx, select.tsx, sheet.tsx, tabs.tsx, textarea.tsx) |
| SC-2 | All 9 pages use shadcn/ui components consistently — no mix of hand-rolled and shadcn/ui elements | PARTIAL | 7 pages clean; OCSRPage uses hand-rolled motion.a buttons and card divs; HomePage uses hand-rolled gradient CTA button and feature card divs |
| SC-3 | Dark mode toggle visible in header/nav using lucide Sun/Moon icons | VERIFIED | Header.tsx line 197: `<ThemeToggle>` rendered in desktop nav; line 230: in Sheet for mobile; uses `import { Menu, Moon, Sun } from "lucide-react"` |
| SC-4 | Theme preference persists in localStorage; system preference auto-detected on first visit | VERIFIED | AppContext.tsx lines 68-75: reads localStorage on mount, falls back to `window.matchMedia("(prefers-color-scheme: dark)")`. FOUC script in index.html lines 57-68 prevents flash. Tests THEME-02 and THEME-03 pass. |
| SC-5 | Icons throughout the app use lucide-react instead of react-icons/@fortawesome | VERIFIED | 47 source files import from lucide-react; react-icons and @fortawesome removed from package.json; icon-migration.test.ts has 6 passing real assertions |

**Score:** 4/5 success criteria fully verified (SC-2 partial)

---

## Required Artifacts (from PLAN frontmatter)

### Plan 04-01 Artifacts

| Artifact | Status | Details |
|----------|--------|---------|
| `frontend/components.json` | VERIFIED | Exists; `rsc: false`, `style: new-york`, correct aliases, `iconLibrary: lucide` |
| `frontend/src/lib/utils.ts` | VERIFIED | Exports `cn()` using clsx + tailwind-merge; 3 passing unit tests |
| `frontend/src/styles/tailwind.css` | VERIFIED | Contains `--background: oklch(...)` and full shadcn/ui variable set for light/dark |
| `frontend/index.html` | VERIFIED | FOUC prevention blocking script at lines 57-68 reads `localStorage.getItem('darkMode')` |
| `frontend/src/context/AppContext.tsx` | VERIFIED | Line 73: `window.matchMedia("(prefers-color-scheme: dark)").matches` |
| `frontend/src/__tests__/utils.test.ts` | VERIFIED | 3 real assertions (not stubs), all passing |
| `frontend/src/__tests__/components/button.test.tsx` | VERIFIED | 1 real render test + 3 todos |
| `frontend/src/__tests__/components/card.test.tsx` | VERIFIED | 1 real render test + 2 todos |
| `frontend/src/__tests__/components/dialog.test.tsx` | VERIFIED | Exists with 3 todos (dialog has no render test — COMP-04 is stub-only) |
| `frontend/src/__tests__/components/input.test.tsx` | VERIFIED | 1 real Input render test + 4 todos |
| `frontend/src/__tests__/icon-migration.test.ts` | VERIFIED | 6 real assertions, all passing |
| `frontend/src/__tests__/context/theme.test.tsx` | VERIFIED | THEME-02 and THEME-03 have real assertions; THEME-01 is still `it.todo` |

### Plan 04-02 Artifacts

| Artifact | Status | Details |
|----------|--------|---------|
| `frontend/package.json` | VERIFIED | `lucide-react: ^0.577.0` present; react-icons and @fortawesome absent |
| `frontend/src/components/common/Header.tsx` | VERIFIED | Imports `{ Menu, Moon, Sun } from "lucide-react"` |

### Plan 04-03 Artifacts

| Artifact | Status | Details |
|----------|--------|---------|
| `frontend/src/pages/HomePage.tsx` | PARTIAL | Imports lucide-react icons but does NOT import `@/components/ui/button` or `@/components/ui/card` |
| `frontend/src/components/common/SMILESInput.tsx` | VERIFIED | Imports both `@/components/ui/input` and `@/components/ui/button` |

### Plan 04-04 Artifacts

| Artifact | Status | Details |
|----------|--------|---------|
| `frontend/src/components/common/Navigation.tsx` | NOT_PRIMARY | Plan listed Navigation.tsx but SUMMARY shows only Header.tsx was modified; Sheet is implemented in Header.tsx |
| `frontend/src/components/common/Header.tsx` | VERIFIED | Uses `Sheet, SheetTrigger, SheetContent, SheetHeader, SheetTitle` from `@/components/ui/sheet` (line 8) |

---

## Key Link Verification

| From | To | Via | Status | Details |
|------|----|-----|--------|---------|
| `frontend/index.html` | localStorage | blocking script reads `darkMode` key | VERIFIED | `localStorage.getItem('darkMode')` at line 59 |
| `frontend/src/context/AppContext.tsx` | window.matchMedia | system preference detection | VERIFIED | `window.matchMedia("(prefers-color-scheme: dark)").matches` at line 73 |
| `frontend/vite.config.mts` | `frontend/src` | @ path alias | VERIFIED | `resolve.alias: { "@": path.resolve(__dirname, "./src") }` at lines 15-16; also in test alias at lines 42-43 |
| all component/page TSX files | `@/components/ui/button` | shadcn Button imports | PARTIAL | 38 files import Button; OCSRPage and HomePage page-level files do not |
| all component/page TSX files | `@/components/ui/card` | shadcn Card imports | PARTIAL | Only 1 component file imports Card directly; OCSRPage and HomePage use hand-rolled card divs |
| `frontend/src/components/common/Header.tsx` | `@/components/ui/sheet` | Sheet for mobile navigation | VERIFIED | `Sheet, SheetTrigger, SheetContent` imported and used |

---

## Requirements Coverage

| Requirement | Source Plan | Description | Status | Evidence |
|-------------|-------------|-------------|--------|---------|
| COMP-01 | 04-01 | shadcn/ui installed and configured for Vite + Tailwind v4 | SATISFIED | components.json, src/components/ui/ (8 components), cn() utility, path aliases |
| COMP-02 | 04-03 | Core buttons replaced with shadcn/ui Button | PARTIAL | 38 files use shadcn Button; OCSRPage and HomePage page files use hand-rolled button-like elements |
| COMP-03 | 04-03 | Cards and content containers replaced with shadcn/ui Card | PARTIAL | Most components use Card; OCSRPage and HomePage page-level containers are hand-rolled divs |
| COMP-04 | 04-04 | Modals and dialogs replaced with shadcn/ui Dialog/Sheet | SATISFIED | Sheet implemented in Header for mobile nav; no dialog patterns existed in codebase (confirmed by plan) |
| COMP-05 | 04-03 | Form inputs replaced with shadcn/ui Input/Select/Textarea | SATISFIED | 24 files use shadcn Input; remaining native inputs in InChIView/RInChIView are checkbox/radio types (explicitly excluded per plan) |
| COMP-06 | 04-02 | Icons migrated from react-icons to lucide-react | SATISFIED | 47 files use lucide-react; react-icons/fortawesome removed from package.json; 6 automated tests pass |
| COMP-07 | 04-04 | Consistent component usage across all 9 pages | PARTIAL | 7 pages consistent; OCSRPage and HomePage page-level still have hand-rolled elements |
| THEME-01 | 04-01, 04-04 | Dark mode toggle available in UI (header/nav) | SATISFIED | ThemeToggle component in Header.tsx rendered in both desktop nav (line 197) and mobile Sheet (line 230) |
| THEME-02 | 04-01 | System preference auto-detection (prefers-color-scheme) | SATISFIED | AppContext.tsx line 73; FOUC script in index.html also detects system preference; THEME-03 test passes |
| THEME-03 | 04-01 | Theme preference persisted in localStorage | SATISFIED | AppContext.tsx lines 94-95; THEME-02 test passes |
| THEME-04 | 04-01 | CSS variable-based theming using shadcn/ui approach | SATISFIED | tailwind.css uses full shadcn/ui OKLCH variable set (--background, --foreground, --primary, etc.) with @theme inline block |
| THEME-05 | 04-04 | All pages and components render correctly in both light and dark modes | NEEDS HUMAN | All dark: variant classes in place; visual correctness requires browser verification |

---

## Anti-Patterns Found

| File | Line | Pattern | Severity | Impact |
|------|------|---------|----------|--------|
| `frontend/src/pages/HomePage.tsx` | 230 | Hand-rolled `<a>` styled as primary button with inline gradient classes instead of shadcn Button | Warning | Inconsistency with COMP-02/COMP-07; no functional regression |
| `frontend/src/pages/OCSRPage.tsx` | 103-124 | `motion.a` elements styled as buttons (no shadcn Button) | Warning | Inconsistency with COMP-02/COMP-07; no functional regression |
| `frontend/src/pages/OCSRPage.tsx` | 133-168 | Hand-rolled card divs with `bg-white dark:bg-slate-800/80` instead of shadcn Card | Warning | Uses isDarkMode-style dark: ternary in className; inconsistency with COMP-03/COMP-07 |
| `frontend/src/__tests__/context/theme.test.tsx` | 27 | `it.todo("THEME-01: toggleDarkMode switches isDarkMode state")` | Info | THEME-01 has no automated assertion; toggleDarkMode itself is tested indirectly by THEME-02 |
| `frontend/src/__tests__/components/dialog.test.tsx` | all | All 3 tests are `it.todo` — no real Dialog assertions | Info | Dialog component installed but untested |
| `frontend/src/components/common/Header.tsx` | 87,90 | Two `isDarkMode ? "text-..." : "text-..."` className ternaries in ThemeToggle for Sun/Moon icon colors | Info | Legitimate non-CSS use (icon-specific color that cannot use dark: variant without structural change); not a blocking issue |

---

## Human Verification Required

### 1. Full visual inspection: light and dark modes across all 9 pages

**Test:** Start `npm run dev`, visit all 9 pages, toggle between light and dark mode on each.
**Expected:** Colors are correct (sky/slate palette), text is readable, cards have visible boundaries, buttons are properly styled, no invisible elements.
**Why human:** Visual correctness of OKLCH color values and contrast ratios cannot be verified programmatically.

### 2. Mobile Sheet navigation

**Test:** Open DevTools, set viewport to 375px, click the hamburger icon in the header.
**Expected:** shadcn/ui Sheet slides in from the right. Navigation links are visible and clickable. Clicking a link closes the Sheet. Mobile theme toggle inside the Sheet switches modes.
**Why human:** Sheet slide animation, focus trap, and link-close behavior require browser rendering.

### 3. FOUC prevention verification

**Test:** Enable dark mode, then hard-refresh the page (Cmd+Shift+R on Mac / Ctrl+Shift+R on Windows).
**Expected:** Page renders directly in dark mode — no visible flash of light theme before dark mode is applied.
**Why human:** Flash of Unstyled Content is a rendering timing issue only observable in a real browser.

### 4. Tab animations on 4 tabbed pages

**Test:** Visit ChemPage, DepictPage, ConvertPage, ToolsPage. Switch between tabs using both desktop sidebar and mobile dropdown.
**Expected:** Tab selection works correctly. On desktop, the animated tab indicator moves smoothly. No console errors.
**Why human:** LayoutGroup animation quality and motion wrapper integrity require visual inspection.

---

## Gaps Summary

Two gaps share the same root cause: **OCSRPage and HomePage page-level files were not migrated to shadcn/ui components during Plan 04-03.**

The component files rendering the actual tool content (OCRView.tsx etc.) were migrated correctly. But the outer page wrapper components for OCSRPage and HomePage still contain hand-rolled button-like elements and card-like containers. This breaks success criterion SC-2 ("no mix of hand-rolled and shadcn/ui buttons, cards") and leaves COMP-02, COMP-03, and COMP-07 partially unsatisfied.

**Scope:** The gaps are limited to 2 of 9 pages at the page wrapper level. All 35+ component files are correctly migrated. The gaps are styling inconsistencies, not functional failures — the dark mode toggle, theme persistence, system preference detection, icon migration, CSS variables, and all form inputs are fully working.

**Priority:** Both gaps are in page-level wrappers, not in tool functionality. They are Warning-severity (not Blocker) since they do not prevent any feature from working.

---

_Verified: 2026-03-13T08:48:22Z_
_Verifier: Claude (gsd-verifier)_
