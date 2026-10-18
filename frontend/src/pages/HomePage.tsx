import { lazy, Suspense, useState } from "react";
import { Link } from "react-router-dom";
import {
  BookOpen,
  Camera,
  FlaskConical,
  LineChart,
  Play,
  RefreshCw,
  SlidersHorizontal,
} from "lucide-react";
import { get2DDepictionUrlEnhanced } from "@/services/depictService";
import { QuietBoundary } from "@/components/common/QuietBoundary";
// three.js loads only with the home page
const Molecule3D = lazy(() => import("@/components/common/Molecule3D"));

interface Mol {
  name: string;
  smiles: string;
}

const MOLS: Mol[] = [
  { name: "Caffeine", smiles: "CN1C=NC2=C1C(=O)N(C(=O)N2C)C" },
  { name: "Artemisinin", smiles: "CC1CCC2C(C(=O)OC3C24C1CCC(O3)(OO4)C)C" },
  { name: "Quercetin", smiles: "O=C1C(O)=C(Oc2cc(O)cc(O)c12)c1ccc(O)c(O)c1" },
  { name: "Morphine", smiles: "CN1CCC23C4C1CC5=C2C(=C(C=C5)O)OC3C(C=C4)O" },
  { name: "Penicillin G", smiles: "CC1(C)SC2C(NC(=O)Cc3ccccc3)C(=O)N2C1C(=O)O" },
  { name: "Salicin", smiles: "OCC1OC(Oc2ccccc2CO)C(O)C(O)C1O" },
  { name: "Capsaicin", smiles: "COc1cc(CNC(=O)CCCC/C=C/C(C)C)ccc1O" },
  { name: "Nicotine", smiles: "CN1CCCC1c1cccnc1" },
  { name: "Vanillin", smiles: "O=Cc1ccc(O)c(OC)c1" },
  { name: "Serotonin", smiles: "NCCc1c[nH]c2ccc(O)cc12" },
  { name: "Coumarin", smiles: "O=c1ccc2ccccc2o1" },
  { name: "Menthol", smiles: "CC(C)C1CCC(C)CC1O" },
  { name: "Strychnine", smiles: "O=C1CC2OCC=C3CN4CCC56C4CC3C2C6N1c1ccccc15" },
  { name: "Resveratrol", smiles: "Oc1ccc(/C=C/c2cc(O)cc(O)c2)cc1" },
];

const PICK = MOLS.slice(0, 5);
const RING = MOLS.slice(5, 13);
const RING_ANGLES = [0, 45, 90, 135, 180, 225, 270, 315];

const SECTIONS = [
  { to: "/chem", name: "Analysis", icon: FlaskConical },
  { to: "/convert", name: "Convert", icon: RefreshCw },
  { to: "/depict", name: "Depict", icon: LineChart },
  { to: "/tools", name: "Tools", icon: SlidersHorizontal },
  { to: "/ocsr", name: "OCSR", icon: Camera },
];

const draw = (smiles: string, size = 220) =>
  get2DDepictionUrlEnhanced(smiles, { width: size, height: size, CIP: false });
const depictHref = (smiles: string) => `/depict/2ddepiction?smiles=${encodeURIComponent(smiles)}`;
const ext = { target: "_blank", rel: "noopener noreferrer" } as const;
const NewTab = () => <span className="sr-only"> (opens in a new tab)</span>;

/** Home: the molecule in a glass lens, 2D drawings around it, then a belt, the sections and the proof. */
const HomePage = () => {
  const [mol, setMol] = useState<Mol>(MOLS[0]!); // a fixed, non-empty list;

  return (
    <div className="pc-home">
      <section className="pc-wrap pc-hero pc-hero-lens" aria-labelledby="home-h">
        <div className="pc-hero-copy">
          <h1 id="home-h">
            Open tools. <span>Seamless access.</span>
          </h1>
          <p className="pc-hero-lead">Visualize, analyze, transform. One open API.</p>
          <div className="pc-row pc-hero-cta">
            <Link className="clay clay-primary" to="/depict/2ddepiction">
              <Play className="pc-i" aria-hidden="true" />
              Start with Depict
            </Link>
            <a
              className="clay clay-outline"
              href="https://api.naturalproducts.net/latest/docs"
              {...ext}
            >
              <BookOpen className="pc-i" aria-hidden="true" />
              API docs
              <NewTab />
            </a>
          </div>
        </div>

        <div className="pc-lens-stage">
          <span className="pc-lens-aura" aria-hidden="true" />
          <span className="pc-lens-ring pc-lens-ring-out" aria-hidden="true" />
          <span className="pc-lens-ring" aria-hidden="true" />
          <ul className="pc-orbit" aria-label="2D drawings">
            {RING.map((m, i) => (
              <li key={m.name} style={{ "--a": `${RING_ANGLES[i]}deg` } as React.CSSProperties}>
                <Link to={depictHref(m.smiles)} className="pc-orbit-chip">
                  <img src={draw(m.smiles, 200)} alt="" width={200} height={200} />
                  <span>{m.name}</span>
                </Link>
              </li>
            ))}
          </ul>
          <div className="pc-lens">
            <div className="pc-lens-3d" role="img" aria-label={`${mol.name} in 3D`}>
              <QuietBoundary>
                <Suspense fallback={null}>
                  <Molecule3D smiles={mol.smiles} className="pc-lens-canvas" />
                </Suspense>
              </QuietBoundary>
            </div>
          </div>
          <div className="pc-row pc-hero-pick" role="group" aria-label="Choose a molecule">
            {PICK.map((m) => (
              <button
                key={m.name}
                type="button"
                className="clay clay-outline pc-chip"
                aria-pressed={m.smiles === mol.smiles}
                onClick={() => setMol(m)}
              >
                {m.name}
              </button>
            ))}
          </div>
        </div>
      </section>

      <div className="pc-wrap pc-home-body">
        <nav className="pc-sections-bar glass" aria-label="Sections">
          {SECTIONS.map(({ to, name, icon: Icon }) => (
            <Link key={to} to={to}>
              <Icon className="pc-i" aria-hidden="true" />
              {name}
            </Link>
          ))}
        </nav>

        <section className="pc-band glass" aria-labelledby="band-h">
          <svg className="pc-band-icon" viewBox="0 0 24 24" aria-hidden="true">
            <path d="M14 10l-2 1m0 0l-2-1m2 1v2.5M20 7l-2 1m2-1l-2-1m2 1v2.5M14 4l-2-1-2 1M4 7l2-1M4 7l2 1M4 7v2.5M12 21l-2-1m2 1l2-1m-2 1v-2.5M6 18l-2-1v-2.5M18 18l2-1v-2.5" />
          </svg>
          <div>
            <h2 id="band-h">Designed for Chemists. Fueled by Open-Source.</h2>
            <p>
              Built on the open-source Cheminformatics Microservice API by the{" "}
              <a href="https://cheminf.uni-jena.de" {...ext}>
                Steinbeck Lab
                <NewTab />
              </a>{" "}
              at{" "}
              <a href="https://www.uni-jena.de/en" {...ext}>
                Friedrich Schiller University Jena
                <NewTab />
              </a>
              .
            </p>
          </div>
        </section>
      </div>
    </div>
  );
};

export default HomePage;
