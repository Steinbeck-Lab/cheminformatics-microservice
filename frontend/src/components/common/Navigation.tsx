import { useLayoutEffect, useRef } from "react";
import { NavLink, useLocation } from "react-router-dom";
import {
  Camera,
  FlaskConical,
  House,
  LineChart,
  RefreshCw,
  SlidersHorizontal,
  Users,
} from "lucide-react";

const navLinks = [
  { to: "/home", label: "Home", icon: House, exact: true },
  { to: "/chem", label: "Analysis", icon: FlaskConical },
  { to: "/convert", label: "Convert", icon: RefreshCw },
  { to: "/depict", label: "Depict", icon: LineChart },
  { to: "/tools", label: "Tools", icon: SlidersHorizontal },
  { to: "/ocsr", label: "OCSR", icon: Camera },
  { to: "/about", label: "About", icon: Users },
];

/**
 * The seven sections. Desktop: links in the header bar with an indigo clay pill that slides
 * under the current one (saffron glow beneath). Menu: the same links as a list in the phone sheet.
 */
const Navigation = ({
  isMobile = false,
  closeMenu = () => {},
}: {
  isMobile?: boolean;
  closeMenu?: () => void;
}) => {
  const navRef = useRef<HTMLElement>(null);
  const pillRef = useRef<HTMLSpanElement>(null);
  const { pathname } = useLocation();

  // Slide the pill to the current link (desktop only).
  useLayoutEffect(() => {
    const nav = navRef.current;
    const pill = pillRef.current;
    if (isMobile || !nav || !pill) return;
    const place = () => {
      const cur = nav.querySelector<HTMLElement>('a[aria-current="page"]');
      nav.classList.toggle("has-pill", !!cur);
      if (!cur) return;
      pill.style.left = `${cur.offsetLeft}px`;
      pill.style.width = `${cur.offsetWidth}px`;
    };
    place();
    document.fonts?.ready.then(place);
    window.addEventListener("resize", place);
    return () => window.removeEventListener("resize", place);
  }, [isMobile, pathname]);

  const links = navLinks.map(({ to, label, icon: Icon, exact }) => (
    <NavLink key={to} to={to} end={exact} onClick={isMobile ? closeMenu : undefined}>
      <Icon className="pc-i" aria-hidden="true" />
      {label}
    </NavLink>
  ));

  if (isMobile) return <nav aria-label="Sections (menu)">{links}</nav>;

  return (
    <nav ref={navRef} className="pc-sections" aria-label="Sections">
      {links}
      <span ref={pillRef} className="pc-npill" aria-hidden="true" />
    </nav>
  );
};

export default Navigation;
