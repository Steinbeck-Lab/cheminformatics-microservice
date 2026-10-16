# Product

<!-- impeccable:product-schema 1 -->

## Platform

web

## Users

Three groups, equally important (confirmed):

- **Bench chemists.** Medicinal and natural-product chemists who do not write code. They paste, draw or photograph a structure and need a fast, correct answer: a depiction, a conversion, a property, a check.
- **Cheminformaticians and developers.** They try a function in the UI to see its input and output, then call the same endpoint from their own code through the public API.
- **Students and teachers.** Courses and workshops that use the tool to learn what SMILES, InChI, depictions, descriptors and filters mean.

## Product Purpose

Cheminformatics Microservice gives free, open web access to the main open cheminformatics toolkits through one UI and one REST API. People convert, depict, analyse and validate chemical structures without installing CDK, RDKit or Open Babel.

The UI has two jobs, equally weighted (confirmed): it is a daily work tool people return to, and it is the front door that shows what the API can do so people start using the API.

Success means a user gets a trustworthy result in seconds, understands which toolkit produced it, and can repeat the same call through the API.

## Positioning

- One place that unifies CDK, RDKit and Open Babel behind the same inputs, so results from different toolkits can be compared side by side.
- Hosts DECIMER for optical chemical structure recognition (image to structure), from the same lab that builds DECIMER.
- Natural-product focus from the Steinbeck Lab: NP-likeness, sugar detection and removal (Sugar Removal Utility), COCONUT preprocessing.
- Every UI function maps to a documented public API endpoint (`api.naturalproducts.net`), and the work is peer-reviewed in the Journal of Cheminformatics.

## Operating Context

- Input is usually a SMILES string pasted from another tool, a structure drawn in the embedded Ketcher editor, a name/CAS/formula looked up via PubChem, a molfile upload, or an image of a structure (OCSR).
- Output is a depiction (2D SVG/PNG, 3D), a converted identifier (SMILES, InChI, InChIKey, RInChI, CXSMILES, molfile, and more), a table of descriptors or filter results, or a classification.
- Users often copy results back into papers, notebooks, spreadsheets or their own code, and download depictions.
- Sections in the current UI: Analysis (`/chem`), Convert (`/convert`), Depict (`/depict`), Tools (`/tools`), OCSR (`/ocsr`), About, plus Terms and Privacy.
- Backend: FastAPI with CDK (via JPype), RDKit, Open Babel and the ChEMBL structure pipeline. Production API at `https://api.naturalproducts.net`, development at `dev.api.naturalproducts.net`. API docs at `/latest/docs`; documentation site built with VitePress (`docs/`).

## Capabilities and Constraints

- Analysis tools: Stereoisomers, HOSE Codes, Standardize, Functional Groups (Ertl), Tautomers, Fix Radicals, Descriptors, NP-likeness, Similarity (Tanimoto), Structure Finder, Check Structure, All Filters, COCONUT Preprocessing, ClassyFire.
- Convert: Format Conversion, 2D Coordinates, 3D Coordinates.
- Depict: Structure Explorer, 2D Depiction (CDK or RDKit, batch), 3D Depiction, Draw a Structure (Ketcher).
- Tools: Sugar Detection and Removal, Structure Generation (Surge), InChI Converter, RInChI Converter.
- OCSR: image to structure through DECIMER.
- STOUT (SMILES to IUPAC name) is unavailable for legal reasons. Do not present it as a feature.
- ClassyFire depends on the external Wishart Lab service, which is often slow or down. The UI must handle its failure gracefully.
- Frontend stack: React 19, Vite, Tailwind CSS v4, shadcn/ui on Radix primitives, motion, three.js for the home-page hero. Ketcher is embedded as a standalone iframe.
- Usage analytics run through Matomo.

## Brand Commitments

All confirmed as binding for every future design:

- Keep the name "Cheminformatics Microservice" and the current hexagon logo (`public/img/logo.svg`, `logo.png`, `logo_small.png`, `logo_small_inverted.png`).
- Show the DFG (Deutsche Forschungsgemeinschaft) funding credit and the Friedrich Schiller University Jena / Steinbeck Lab credit.
- Show the papers so users can cite the tool, on the About page (the owner removed the citation row from the footer):
  - Chandrasekhar, V., Sharma, N., Schaub, J. et al. Cheminformatics Microservice: unifying access to open cheminformatics toolkits. J Cheminform 15, 98 (2023). https://doi.org/10.1186/s13321-023-00762-4
  - Rajan, K., Chandrasekhar, V., Sharma, N. et al. Cheminformatics Microservice V3: a web portal for chemical structure manipulation and analysis. J Cheminform 17, 142 (2025). https://doi.org/10.1186/s13321-025-01094-1
  - Software DOI: https://doi.org/10.5281/zenodo.7745987
- Free, open source (MIT) and no login. The UI must never look like a paid product, and never add sign-up walls, pricing or upsell.
- Visual constraints set by the owner (binding for the redesign):
  - Fonts: Inter for UI text, JetBrains Mono for SMILES and identifiers. Display layer (headlines, page titles, panel labels): Antonio, light condensed capitals with wide tracking; no serif and no italics (owner decision, 2026-10-15).
  - Only buttons use claymorphism (soft, raised, pressable). Nothing else is clay.
  - Surfaces (panels, bento cells, header) use glassmorphism / liquid glass.
  - Tool pages use bento-grid layouts.
  - Atmosphere: "Silk crossings": three wide translucent silk ribbons crossing the page diagonally (porcelain grey, pale mauve, one saffron), deepening where they overlap.
  - Theme: "Porcelain" with a saffron accent, light and dark. Colours come from the "Computational Color Chord Counterpoint" palette: indigo silk #2B2F7F (primary button), saffron silk #F2B31A (secondary accent, about 10–15% of the screen, never as text on white), madder silk #B23A65 (warnings, Clear), mauve #5B4A72, on a porcelain-white ground.
  - The header turns into clear glass when the page is scrolled. The footer has no citation row.
  - Keep today's navigation: the floating header bar with Home plus six sections (Analysis, Convert, Depict, Tools, OCSR, About), tabs inside each section, and the sun/moon theme switch.
  - Keep the original header and footer structure (footer: brand block, Terms/Privacy, four link cards, flask divider, credits), restyled in the new theme rather than copied.
  - The API call stays in the API docs, not on tool pages.

## Evidence on Hand

- The two peer-reviewed papers and the Zenodo DOI above; `CITATION.cff` in the repo root.
- Public API documentation and the VitePress docs in `docs/`.
- Logos in `public/img/`; the DFG funding logo is referenced in the README and About page.
- No user testimonials, usage numbers or customer logos are on hand. Do not invent them.

## Product Principles

1. **The molecule and the result come first.** Every screen exists to get a structure in and a trustworthy result out.
2. **Show the toolkit.** Say which toolkit (CDK, RDKit, Open Babel, DECIMER) produced a result, and make comparing toolkits easy.
3. **Every UI action is an API call people can repeat.** Make the path from "it worked in the UI" to "I can call this from code" short and visible.
4. **Serve the expert and the learner at once.** Fast and dense enough for daily users, clear enough for a student meeting SMILES for the first time.
5. **Open science, plainly.** Free, no login, credited and citable. Never fake commercial polish.

## Accessibility & Inclusion

Target: WCAG 2.2 AA (confirmed). The interface must work fully with a keyboard and a screen reader, respect reduced-motion settings, and keep chemical data (SMILES, identifiers, tables) readable and copyable.
