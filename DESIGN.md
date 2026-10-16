---
name: Cheminformatics Microservice
description: Porcelain workbench under silk light. Glass bento cells over three translucent silk ribbons; the molecule on its white stage is the only solid thing.
colors:
  porcelain: "#f6f4f1"
  porcelain-dark: "#171517"
  ink: "#221e26"
  ink-2: "#4f4856"
  head: "#4a3d5e"
  glass-solid: "#fcfbf9"
  indigo: "#2b2f7f"
  indigo-dark: "#aeb1f0"
  madder: "#b23a65"
  madder-dark: "#e06a93"
  madder-ink-dark: "#f4a3c0"
  saffron: "#f2b31a"
  saffron-ink: "#8a5a00"
  saffron-ink-dark: "#ffc93d"
  mauve: "#5b4a72"
  stage: "#ffffff"
typography:
  display:
    fontFamily: "Antonio Variable, Arial Narrow, Inter, sans-serif"
    fontSize: "clamp(42px, 4.5vw, 66px)"
    fontWeight: 300
    lineHeight: 1.1
    letterSpacing: "0.05em"
  title:
    fontFamily: "Antonio Variable, Arial Narrow, Inter, sans-serif"
    fontSize: "30px"
    fontWeight: 300
    lineHeight: 1.1
    letterSpacing: "0.06em"
  label:
    fontFamily: "Antonio Variable, Arial Narrow, Inter, sans-serif"
    fontSize: "12.5px"
    fontWeight: 500
    lineHeight: 1.2
    letterSpacing: "0.24em"
  body:
    fontFamily: "Inter, system-ui, sans-serif"
    fontSize: "15px"
    fontWeight: 400
    lineHeight: 1.55
    letterSpacing: "normal"
  identifier:
    fontFamily: "JetBrains Mono, ui-monospace, monospace"
    fontSize: "13.5px"
    fontWeight: 400
    lineHeight: 1.5
    letterSpacing: "normal"
rounded:
  ctl: "14px"
  stage: "16px"
  cell: "22px"
  pill: "999px"
spacing:
  gutter: "16px"
  cell: "20px"
  page: "24px"
components:
  button-primary:
    backgroundColor: "{colors.indigo}"
    textColor: "#ffffff"
    rounded: "{rounded.ctl}"
    height: "40px"
    padding: "0 16px"
  button-outline:
    backgroundColor: "#ffffff"
    textColor: "{colors.madder}"
    rounded: "{rounded.ctl}"
    height: "40px"
  button-danger:
    backgroundColor: "{colors.madder}"
    textColor: "#ffffff"
    rounded: "{rounded.ctl}"
    height: "40px"
  chip:
    backgroundColor: "#ffffff"
    textColor: "{colors.madder}"
    rounded: "{rounded.pill}"
    height: "34px"
  stage:
    backgroundColor: "{colors.stage}"
    rounded: "{rounded.stage}"
---

# Design System: Cheminformatics Microservice

## Overview

A porcelain workbench under silk light. Every surface is frosted glass floating over three wide translucent silk ribbons (porcelain grey, madder, saffron) that cross the page diagonally and deepen where they overlap. The only solid, opaque thing on a page is the molecule: drawn on a white stage, so depictions read exactly as they will in a paper. Clay is reserved for things you press. The world refuses the blue glass-and-gradient "science dashboard".

The interface serves three groups at once (bench chemists, developers, students), so the system is dense where work happens (tool pages: bento cells, toolkit badges, identifiers) and nearly wordless where people arrive (Home: four-word headline, one 3D molecule). Every result names the toolkit that produced it.

Brand commitments in `PRODUCT.md` are binding: name and hexagon logo, DFG and FSU Jena / Steinbeck Lab credit, citations on About, free and no login.

## Colors

Strategy: **restrained porcelain with three named accents**. Neutrals carry the page; colour is assigned by role, never by decoration.

- **Porcelain** (`#f6f4f1`, dark `#171517`): the ground. **Ink** `#221e26` for text, **Ink-2** `#4f4856` for secondary text, **Head** `#4a3d5e` (mauve-tinted) for headings.
- **Indigo** (`#2b2f7f`, dark `#aeb1f0`): the one primary action per screen, the active nav pill, selected states, focus ring in light mode.
- **Madder** (`#b23a65`): secondary actions (outline buttons, chips), Clear/destructive, section-label dots, panel accents, the "how to cite" lead-ins. Used widely by the owner's request.
- **Saffron** (`#f2b31a`): fills and strokes only, about 10–15% of a screen: the selected dot inside chips and segmented controls, the underline under active tabs and the headline, the tracer line, the lens aura, text selection. **Saffron is never text on white**; saffron-toned text uses **Saffron-ink** `#8a5a00` (dark `#ffc93d`).
- **Mauve** `#5b4a72`: quiet supporting tone (charts, the lens aura).
- **Stage** `#ffffff`: depictions always sit on pure white in both themes.

Dark mode swaps every token under `.dark` on `<html>`; components read tokens only, so nothing hard-codes a palette colour. Text pairs are guarded by `frontend/src/__tests__/styles/tokens.test.ts` (WCAG 4.5:1, focus ring 3:1).

## Typography

Three faces, three jobs:

- **Antonio** (self-hosted, `@fontsource-variable/antonio`): the display layer only: page headlines, tool titles, section and panel labels, the sections strip. Always **uppercase**, light (300) or medium (500), wide tracking (0.05em on headlines, 0.24em on labels). Owner decision: no serif, no italics anywhere.
- **Inter**: all running text, controls, buttons, form labels, notes.
- **JetBrains Mono**: SMILES, InChI, identifiers, endpoints, numeric descriptor values (tabular numerals).

Ramp: hero headline `clamp(42px, 4.5vw, 66px)`; tool and About titles 30–36px; panel labels 12.5px tracked capitals with a hairline rule to the right that ends in a madder dot; body 15px/1.55; notes 12.5–13px.

## Layout

- One container for every page: `.pc-wrap`, max 1312px, 24px side padding (16px on phones). All pages share the same left and right edges.
- **Tool pages are bentos**: a 12-column grid with 16px gaps. Input cell on the left (7 columns), options or "About" cell on the right (5), results full width below. Under 1024px every cell is full width.
- Section frame: glass tab bar (saffron underline on the open tab; tabs wrap on phones) or, for Analysis, a sticky glass sidebar with Categories / All tools; then a glass title cell; then the tool.
- Home: two columns at 1024px+ (5/12 copy, 7/12 lens); the lens group is sized from one variable (`--orbit`) and the chip size, so the ring, lens, chips and picker never overlap at any width.

## Elevation & Depth

Depth comes from glass, not from shadows stacked on shadows:

- **Glass cell**: translucent porcelain (`--glass`, 66% white / 58% dark), `blur(14px) saturate(160%)` (never above 16px), 1px rim, a top-left sheen, and one soft offset drop shadow. Without `backdrop-filter` support it falls back to solid `--glass-solid`.
- **Clay** (buttons only): raised face with an inner highlight, an inner lower shade and a short offset drop; pressed moves down 1px (only with motion allowed) and inverts the shading. Disabled clay is flat.
- **Stage**: white, 1px hairline, no glass.
- The silk sits behind everything (`z-index: -10` inside an isolated layout) and multiplies in light, screens in dark.

## Shapes

Radii: 14px for controls and buttons, 16px for stages and inner boxes, 22px for cells, pills for chips, nav pill and segmented tracks. Molecular forms (the hexagon logo, the lens circle, the ribbons) provide the organic shapes; nothing else is decorative.

## Components

- **Buttons** (`<Button>`): default = indigo clay (one per screen), `outline` = white clay with madder rim and label (secondary), `destructive` = solid madder clay (Clear), `ghost` = flat icon actions. Clay states live inside the Tailwind utilities (`&:where(:hover|:active|:disabled)`) so page classes can still override and disabled stays flat.
- **Cell** (`porcelain-ui.tsx`): a glass `<section>` labelled by its uppercase panel label, optional meta at the right of the label row.
- **Segmented control** (`Seg`): clay track pressed in, raised clay thumb with a saffron dot; arrow keys move the selection.
- **Switch**: flat track, clay knob; indigo when on.
- **SMILES input**: mono well, Paste (outline), Clear (madder), example chips (saffron dot on the active one), recent chips.
- **Stage**: white box, image scales to fit (never cropped), small grey mono tag for size/toolkit.
- **Toolkit compare**: two stages side by side, the active one ringed in indigo.
- **Result cards**: inset card, white stage, mono SMILES line, ghost icon actions (copy, download, details).
- **Home lens**: glass lens with a three.js molecule in site colours (C plum, N indigo, O madder, S saffron, glossy clear-coat) that assembles atom by atom, ringed by eight live 2D drawings; a turning saffron/madder/mauve aura behind; a saffron spark on the outer ring.
- **Silk**: still on tool pages; on Home and About the ribbons ripple in a seamless loop with a gliding sheen (paused when the tab is hidden, off under reduced motion).
- **Header**: floating glass pill that hugs its links; turns to clear glass on scroll. **Footer**: compact, brand block, Terms/Privacy, four link cards, flask divider, one credit line.

## Do's and Don'ts

- **Do** put every surface on glass and every action on clay; nothing else is clay.
- **Do** name the toolkit next to every result.
- **Do** keep depictions on a white stage in both themes.
- **Do** keep page text short; Home and About lead with the molecule, not paragraphs.
- **Do** respect `prefers-reduced-motion`: silk still, orbit still, molecule pre-assembled, no press movement.
- **Do** keep fixed-width chemistry text intact: the molfile terminator is `M  END` with two spaces (guarded by a test).
- **Don't** use saffron as text on a light ground.
- **Don't** use serif faces or italics; Antonio only in uppercase for the display layer.
- **Don't** add gradients on text, coloured side stripes on cards, kicker labels above headings, or bouncy easing (use `cubic-bezier(0.16, 1, 0.3, 1)`).
- **Don't** hard-code palette colours in components; read the tokens so dark mode follows.
- **Don't** add a second primary button to a screen.
