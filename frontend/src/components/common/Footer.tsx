import { Link } from "react-router-dom";
import {
  BookOpen,
  Code,
  ExternalLink,
  FileCheck,
  FlaskConical,
  Github,
  Landmark,
  ShieldCheck,
} from "lucide-react";

const LOGO_LIGHT =
  "https://raw.githubusercontent.com/Steinbeck-Lab/cheminformatics-microservice/main/public/img/logo_small.png";
const LOGO_DARK =
  "https://raw.githubusercontent.com/Steinbeck-Lab/cheminformatics-microservice/main/public/img/logo_small_inverted.png";

const resources = [
  {
    href: "https://api.naturalproducts.net/latest/docs",
    icon: BookOpen,
    title: "API Docs",
    desc: "Explore endpoints",
    tint: "var(--bar1)",
  },
  {
    href: "https://github.com/Steinbeck-Lab/cheminformatics-microservice",
    icon: Github,
    title: "GitHub",
    desc: "View code",
    tint: "var(--bar2)",
  },
  {
    href: "https://docs.api.naturalproducts.net",
    icon: Code,
    title: "Guides",
    desc: "Integration help",
    tint: "var(--bar3)",
  },
  {
    href: "https://cheminf.uni-jena.de",
    icon: Landmark,
    title: "Steinbeck Lab",
    desc: "Research group",
    tint: "var(--bar4)",
  },
];

/** A coffee cup with three wisps of steam (the steam stops under reduced motion). */
const Cup = () => (
  <svg className="pc-cup" viewBox="0 0 24 26" role="img" aria-label="coffee">
    <g className="pc-steam">
      <path d="M8 9.6c-1.1-1.2 1.1-2.2 0-3.6" />
      <path className="w2" d="M11 9.6c-1.1-1.2 1.1-2.2 0-3.6" />
      <path className="w3" d="M14 9.6c-1.1-1.2 1.1-2.2 0-3.6" />
    </g>
    <path d="M4.5 12h12v4.2a4.8 4.8 0 0 1-4.8 4.8H9.3a4.8 4.8 0 0 1-4.8-4.8z" />
    <path d="M16.5 13.4h1.4a2.6 2.6 0 0 1 0 5.2h-1.6" />
  </svg>
);

/** Compact glass footer: brand block and four resource cards, a flask divider, one credit line. */
const Footer = () => (
  <footer className="pc-site-foot glass">
    <div className="pc-wrap pc-foot-in">
      <div className="pc-foot-top">
        <div className="pc-foot-brand">
          <Link className="pc-fb-logo" to="/">
            <span className="pc-fb-tile">
              <img className="logo-l" src={LOGO_LIGHT} alt="" width="36" height="30" />
              <img className="logo-d" src={LOGO_DARK} alt="" width="36" height="30" />
            </span>
            <span className="pc-fb-name">
              <span className="pc-fb-1">Cheminformatics</span>
              <span className="pc-fb-2">Microservice</span>
            </span>
          </Link>
          <p className="pc-fb-desc pc-foot-dev">
            Developed by the{" "}
            <a href="https://cheminf.uni-jena.de" target="_blank" rel="noopener noreferrer">
              Steinbeck Lab
            </a>{" "}
            at{" "}
            <a href="https://www.uni-jena.de/en" target="_blank" rel="noopener noreferrer">
              Friedrich Schiller University Jena
            </a>
            , with funding from the Deutsche Forschungsgemeinschaft (DFG).
          </p>
          <div className="pc-fb-legal">
            <Link className="clay clay-outline sm" to="/terms">
              <FileCheck className="pc-i" aria-hidden="true" />
              Terms of Service
            </Link>
            <Link className="clay clay-outline sm" to="/privacy">
              <ShieldCheck className="pc-i" aria-hidden="true" />
              Privacy Policy
            </Link>
          </div>
        </div>
        <nav aria-label="Resources">
          <ul className="pc-foot-cards">
            {resources.map(({ href, icon: Icon, title, desc, tint }) => (
              <li key={href}>
                <a
                  className="pc-fcard glass"
                  style={{ "--tint": tint } as React.CSSProperties}
                  href={href}
                  target="_blank"
                  rel="noopener noreferrer"
                >
                  <span className="pc-ft">
                    <Icon className="pc-i" aria-hidden="true" />
                    {title}
                  </span>
                  <span className="pc-fs">{desc}</span>
                  <ExternalLink className="pc-i pc-ext" aria-hidden="true" />
                  <span className="sr-only"> (opens in a new tab)</span>
                </a>
              </li>
            ))}
          </ul>
        </nav>
      </div>
      <div className="pc-foot-div" aria-hidden="true">
        <span className="pc-hl" />
        <span className="pc-flask">
          <FlaskConical className="pc-i" />
        </span>
        <span className="pc-hl r" />
      </div>
    </div>
    <p className="pc-foot-base">
      <span className="pc-foot-copy">
        © {new Date().getFullYear()} Cheminformatics Microservice UI. Built with <Cup /> and the
        Cheminformatics Microservice API.
      </span>
    </p>
  </footer>
);

export default Footer;
