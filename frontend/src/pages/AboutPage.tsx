import { useState, type ReactNode } from "react";
import { Link } from "react-router-dom";
import {
  ArrowDown,
  BookOpen,
  Camera,
  Check,
  Copy,
  ExternalLink,
  FlaskConical,
  Github,
  LineChart,
  RefreshCw,
  SlidersHorizontal,
} from "lucide-react";
import { Button } from "@/components/ui/button";

const LOGO_LIGHT =
  "https://raw.githubusercontent.com/Steinbeck-Lab/cheminformatics-microservice/main/public/img/logo_small.png";
const LOGO_DARK =
  "https://raw.githubusercontent.com/Steinbeck-Lab/cheminformatics-microservice/main/public/img/logo_small_inverted.png";
const TRACER = "CN1C=NC2=C1C(=O)N(C(=O)N2C)C";

const REFERENCES = [
  {
    id: "v3",
    kind: "Paper, web portal (V3)",
    text: "Rajan, K., Chandrasekhar, V., Sharma, N. et al. Cheminformatics Microservice V3: a web portal for chemical structure manipulation and analysis. J Cheminform 17, 142 (2025).",
    doi: "https://doi.org/10.1186/s13321-025-01094-1",
  },
  {
    id: "v1",
    kind: "Paper, original API",
    text: "Chandrasekhar, V., Sharma, N., Schaub, J. et al. Cheminformatics Microservice: unifying access to open cheminformatics toolkits. J Cheminform 15, 98 (2023).",
    doi: "https://doi.org/10.1186/s13321-023-00762-4",
  },
  {
    id: "sw",
    kind: "Software, all versions",
    text: "Venkata, C., Sharma, N., Schaub, J., Steinbeck, C., & Rajan, K. cheminformatics-microservice [Computer software]. Zenodo.",
    doi: "https://doi.org/10.5281/zenodo.7745987",
    note: { label: "This release, v3.5.0", href: "https://doi.org/10.5281/zenodo.16890410" },
  },
];

const SECTIONS = [
  {
    to: "/chem",
    name: "Analysis",
    what: "descriptors, filters, NP-likeness, ClassyFire",
    icon: FlaskConical,
  },
  {
    to: "/convert",
    name: "Convert",
    what: "SMILES, InChI, molfile, 2D and 3D coordinates",
    icon: RefreshCw,
  },
  { to: "/depict", name: "Depict", what: "2D and 3D depictions, Ketcher drawing", icon: LineChart },
  {
    to: "/tools",
    name: "Tools",
    what: "sugar removal, Surge, InChI and RInChI",
    icon: SlidersHorizontal,
  },
  { to: "/ocsr", name: "OCSR", what: "structure images to SMILES with DECIMER", icon: Camera },
];

const ROUTES = ["/chem", "/convert", "/depict", "/tools", "/ocsr"];

const TOOLKITS = [
  {
    name: "CDK",
    href: "https://cdk.github.io/",
    logo: "https://github.com/cdk/cdk.github.io/raw/source/src/img/logo.png",
    what: "Chemistry Development Kit, the Java library for structural cheminformatics.",
  },
  {
    name: "RDKit",
    href: "https://www.rdkit.org/",
    logo: "https://github.com/rdkit/rdkit/raw/71935ecae38c2f469ac0ed7d81dacd6cd6a19d5a/Docs/Images/logo.png",
    what: "Cheminformatics and machine-learning tools in C++ and Python.",
  },
  {
    name: "Open Babel",
    href: "https://openbabel.org/",
    logo: "https://upload.wikimedia.org/wikipedia/commons/c/ce/Open_Babel_computer_icon.png",
    what: "The chemical toolbox to search, convert, analyse and store data.",
  },
  {
    name: "DECIMER",
    href: "https://decimer.ai/",
    logo: "https://github.com/Steinbeck-Lab/cheminf-jena-logos/blob/main/DECIMER/DECIMER-1.png?raw=true",
    what: "Deep learning for chemical image recognition, from the same lab.",
  },
];

const ALSO = [
  {
    name: "ChEMBL structure pipeline",
    href: "https://github.com/chembl/ChEMBL_Structure_Pipeline",
    what: "standardisation",
  },
  {
    name: "Ketcher",
    href: "https://lifescience.opensource.epam.com/ketcher/",
    what: "structure editor (Apache 2.0)",
  },
  {
    name: "PubChem PUG REST",
    href: "https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest",
    what: "name and identifier lookup",
  },
  { name: "InChI", href: "https://www.inchi-trust.org/", what: "InChI and RInChI generation" },
];

const ext = { target: "_blank", rel: "noopener noreferrer" } as const;
const NewTab = () => <span className="sr-only"> (opens in a new tab)</span>;

function CopyCitation({ text, label }: { text: string; label: string }) {
  const [done, setDone] = useState(false);
  return (
    <Button
      variant="outline"
      size="sm"
      aria-label={`Copy citation: ${label}`}
      onClick={async () => {
        await navigator.clipboard?.writeText(text);
        setDone(true);
        setTimeout(() => setDone(false), 1500);
      }}
    >
      {done ? <Check /> : <Copy />}
      {done ? "Copied" : "Copy"}
    </Button>
  );
}

function Layer({
  id,
  name,
  lead,
  children,
}: {
  id: string;
  name: string;
  lead: ReactNode;
  children: ReactNode;
}) {
  return (
    <li className="pc-stratum glass" aria-labelledby={id}>
      <span className="pc-node" aria-hidden="true" />
      <div className="pc-stratum-head">
        <h2 id={id}>{name}</h2>
        <p>{lead}</p>
      </div>
      <div className="pc-stratum-body">{children}</div>
    </li>
  );
}

/** About: the service cut open, from the browser down to the people and funders who hold it up. */
const AboutPage = () => {
  const v3 = REFERENCES[0]!; // a fixed, non-empty list
  return (
    <div className="pc-wrap pc-about">
      <div className="pc-about-top">
        <header className="pc-about-head glass">
          <span className="pc-about-logo">
            <img className="logo-l" src={LOGO_LIGHT} alt="" width={64} height={54} />
            <img className="logo-d" src={LOGO_DARK} alt="" width={64} height={54} />
          </span>
          <h1>About Cheminformatics Microservice</h1>
          <p className="pc-about-lead">
            One interface and one REST API over the main open cheminformatics toolkits. Convert,
            depict, analyse and validate structures without installing CDK, RDKit or Open Babel.
          </p>
          <ul className="pc-facts" aria-label="Terms of use">
            <li>Free to use</li>
            <li>Open source, MIT licence</li>
            <li>No login</li>
          </ul>
        </header>

        <aside className="pc-cite pc-cell-end glass" aria-labelledby="cite-h">
          <h2 id="cite-h">Cite this service</h2>
          <p className="pc-cite-ref">{v3.text}</p>
          <a className="pc-doi" href={v3.doi} {...ext}>
            {v3.doi.replace("https://", "")}
            <NewTab />
          </a>
          <ul className="pc-cite-more" aria-label="Also cite">
            {REFERENCES.slice(1).map((r) => (
              <li key={r.id}>
                <span>{r.kind}</span>
                <a href={r.doi} {...ext}>
                  {r.doi.replace("https://", "")}
                  <NewTab />
                </a>
              </li>
            ))}
          </ul>
          <div className="pc-row">
            <CopyCitation text={`${v3.text} ${v3.doi}`} label="V3 paper" />
            <a className="clay clay-outline sm" href="#how-to-cite">
              <ArrowDown className="pc-i" aria-hidden="true" />
              All references
            </a>
          </div>
        </aside>
      </div>

      <section className="pc-stack" aria-label="How the service is built">
        <div className="pc-tracer" aria-hidden="true">
          <span className="pc-tracer-line" />
          <span className="pc-tracer-dot" />
        </div>
        <p className="pc-tracer-label">
          <span className="pc-mono">{TRACER}</span> caffeine, on its way down
        </p>
        <ol className="pc-strata" aria-label="The service, layer by layer">
          <Layer
            id="l-ui"
            name="In your browser"
            lead="Paste a SMILES, draw in Ketcher, look up a name or upload an image. Every section answers in seconds and says which toolkit produced the result."
          >
            <ul className="pc-section-links">
              {SECTIONS.map(({ to, name, what, icon: Icon }) => (
                <li key={to}>
                  <Link to={to} className="pc-section-link">
                    <Icon className="pc-i" aria-hidden="true" />
                    <span>
                      {name}
                      <small>{what}</small>
                    </span>
                  </Link>
                </li>
              ))}
            </ul>
          </Layer>

          <Layer
            id="l-api"
            name="One REST API"
            lead="Each action in this interface is a call to a documented, public endpoint. What works here can be repeated from a notebook, a pipeline or your own app."
          >
            <div className="pc-routes">
              <code className="pc-host">api.naturalproducts.net/latest</code>
              <ul aria-label="Endpoint families">
                {ROUTES.map((r) => (
                  <li key={r}>
                    <code>{r}</code>
                  </li>
                ))}
              </ul>
              <a
                className="clay clay-primary"
                href="https://api.naturalproducts.net/latest/docs"
                {...ext}
              >
                <BookOpen className="pc-i" aria-hidden="true" />
                Open the API docs
                <NewTab />
              </a>
            </div>
          </Layer>

          <Layer
            id="l-kits"
            name="Open toolkits"
            lead="The same input goes to CDK, RDKit or Open Babel, so you can compare their answers side by side. DECIMER reads structures from images."
          >
            <ul className="pc-kits">
              {TOOLKITS.map((k) => (
                <li key={k.name}>
                  <a className="pc-kit" href={k.href} {...ext}>
                    <span className="pc-kit-logo">
                      <img src={k.logo} alt="" loading="lazy" />
                    </span>
                    <span className="pc-kit-name">{k.name}</span>
                    <span className="pc-kit-what">{k.what}</span>
                    <NewTab />
                  </a>
                </li>
              ))}
            </ul>
            <p className="pc-also">
              Also inside:{" "}
              {ALSO.map((a, i) => (
                <span key={a.name}>
                  <a href={a.href} {...ext}>
                    {a.name}
                    <NewTab />
                  </a>{" "}
                  <span className="pc-also-what">({a.what})</span>
                  {i < ALSO.length - 1 ? ", " : "."}
                </span>
              ))}
            </p>
          </Layer>

          <Layer
            id="l-host"
            name="A public instance"
            lead={
              <>
                Hosted at{" "}
                <a href="https://www.uni-jena.de/en/" {...ext}>
                  Friedrich Schiller University Jena
                  <NewTab />
                </a>{" "}
                in Germany. Free to use, no login. Run your own copy with Docker.
              </>
            }
          >
            <dl className="pc-host-facts">
              <div>
                <dt>Licence</dt>
                <dd>
                  <a href="https://opensource.org/licenses/MIT" {...ext}>
                    MIT
                    <NewTab />
                  </a>
                  ; Ketcher under Apache 2.0, other components under their own licences.
                </dd>
              </div>
              <div>
                <dt>Source code</dt>
                <dd>
                  <a href="https://github.com/Steinbeck-Lab/cheminformatics-microservice" {...ext}>
                    <Github className="pc-i" aria-hidden="true" />
                    Steinbeck-Lab/cheminformatics-microservice
                    <NewTab />
                  </a>
                </dd>
              </div>
              <div>
                <dt>Guides</dt>
                <dd>
                  <a href="https://docs.api.naturalproducts.net" {...ext}>
                    docs.api.naturalproducts.net
                    <NewTab />
                  </a>
                </dd>
              </div>
            </dl>
          </Layer>

          <Layer
            id="l-people"
            name="People and funding"
            lead={
              <>
                Built and maintained by the{" "}
                <a href="https://cheminf.uni-jena.de" {...ext}>
                  Steinbeck group
                  <NewTab />
                </a>{" "}
                at{" "}
                <a href="https://www.uni-jena.de/en/" {...ext}>
                  Friedrich Schiller University Jena
                  <NewTab />
                </a>
                , who also run{" "}
                <a href="https://naturalproducts.net" {...ext}>
                  Natural Products Online
                  <NewTab />
                </a>{" "}
                and{" "}
                <a href="https://chemaudit.naturalproducts.net/" {...ext}>
                  ChemAudit
                  <NewTab />
                </a>
                .
              </>
            }
          >
            <div className="pc-funders">
              <a className="pc-funder-logo" href="https://www.dfg.de/" {...ext}>
                <img
                  src="https://github.com/Steinbeck-Lab/cheminformatics-microservice/blob/main/docs/public/dfg_logo_schriftzug_blau_foerderung_en.gif?raw=true"
                  alt="Funded by the DFG, Deutsche Forschungsgemeinschaft"
                  width={135}
                  height={44}
                  loading="lazy"
                />
                <NewTab />
              </a>
              <p>
                Funded by the Deutsche Forschungsgemeinschaft (DFG) under{" "}
                <a href="https://nfdi4chem.de/" {...ext}>
                  NFDI4Chem
                  <NewTab />
                </a>{" "}
                (project 441958208) and{" "}
                <a href="https://www.chembiosys.de/en/" {...ext}>
                  ChemBioSys
                  <NewTab />
                </a>{" "}
                (project INF, SFB 1127, project 239748522).
              </p>
            </div>
          </Layer>
        </ol>
      </section>

      <section id="how-to-cite" className="pc-refs glass" aria-labelledby="refs-h">
        <div className="pc-refs-head">
          <h2 id="refs-h">How to cite</h2>
          <p>
            Please cite the software and the papers. The DOIs let readers find the exact code and
            data you used.
          </p>
        </div>
        <ol className="pc-ref-list">
          {REFERENCES.map((r) => (
            <li key={r.id}>
              <p className="pc-ref-text">
                <strong>{r.kind}.</strong> {r.text}
              </p>
              <p className="pc-ref-links">
                <a href={r.doi} {...ext}>
                  {r.doi}
                </a>
                {r.note && (
                  <a className="pc-ref-note" href={r.note.href} {...ext}>
                    {r.note.label}
                    <ExternalLink className="pc-i" aria-hidden="true" />
                    <NewTab />
                  </a>
                )}
              </p>
              <CopyCitation text={`${r.text} ${r.doi}`} label={r.kind} />
            </li>
          ))}
        </ol>
      </section>
    </div>
  );
};

export default AboutPage;
