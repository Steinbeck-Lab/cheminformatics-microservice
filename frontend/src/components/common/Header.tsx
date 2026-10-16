import { useEffect, useState } from "react";
import { Link, useLocation } from "react-router-dom";
import { Menu, Moon, Sun, X } from "lucide-react";
import { useAppContext } from "../../context/AppContext";
import Navigation from "./Navigation";
import { Sheet, SheetClose, SheetContent, SheetTitle, SheetTrigger } from "@/components/ui/sheet";

const LOGO_LIGHT =
  "https://raw.githubusercontent.com/Steinbeck-Lab/cheminformatics-microservice/main/public/img/logo_small.png";
const LOGO_DARK =
  "https://raw.githubusercontent.com/Steinbeck-Lab/cheminformatics-microservice/main/public/img/logo_small_inverted.png";

/** Day/night switch: a pressed-in clay track with a raised clay knob carrying the sun or moon. */
function ThemeSwitch(props: { "aria-label"?: string; "aria-labelledby"?: string }) {
  const { isDarkMode, toggleDarkMode } = useAppContext();
  return (
    <button
      type="button"
      role="switch"
      aria-checked={isDarkMode}
      className="pc-theme"
      onClick={toggleDarkMode}
      {...props}
    >
      <span className="pc-tt-track" aria-hidden="true">
        <i className="pc-st pc-s1" />
        <i className="pc-st pc-s2" />
        <span className="pc-tt-knob">
          <Sun className="pc-i pc-sun" />
          <Moon className="pc-i pc-moon" />
        </span>
      </span>
    </button>
  );
}

/** The floating glass header bar. It turns into clear glass once the page scrolls. */
const Header = () => {
  const [scrolled, setScrolled] = useState(false);
  const [menuOpen, setMenuOpen] = useState(false);
  const location = useLocation();

  useEffect(() => setMenuOpen(false), [location]);

  useEffect(() => {
    const onScroll = () => setScrolled(window.scrollY > 24);
    onScroll();
    window.addEventListener("scroll", onScroll, { passive: true });
    return () => window.removeEventListener("scroll", onScroll);
  }, []);

  return (
    <header className={`pc-top${scrolled ? " scrolled" : ""}`}>
      <div className="pc-wrap">
        <div className="pc-hbar glass">
          <Link to="/" className="pc-brand" aria-label="Cheminformatics Microservice, home">
            <img className="logo-l" src={LOGO_LIGHT} alt="" />
            <img className="logo-d" src={LOGO_DARK} alt="" />
          </Link>
          <span className="pc-vsep" aria-hidden="true" />
          <Navigation />
          <ThemeSwitch aria-label="Dark mode" />

          <Sheet open={menuOpen} onOpenChange={setMenuOpen}>
            <SheetTrigger asChild>
              <button
                type="button"
                className="clay clay-outline icon pc-menu-btn"
                aria-label="Open menu"
              >
                <Menu className="pc-i" aria-hidden="true" />
              </button>
            </SheetTrigger>
            <SheetContent
              side="right"
              showCloseButton={false}
              aria-describedby={undefined}
              className="pc-sheet w-[min(300px,86vw)] gap-1 rounded-l-3xl border-0 border-l border-[var(--rim)] bg-[color-mix(in_srgb,var(--glass-c)_86%,transparent)] px-3 pt-3.5 pb-5 text-foreground backdrop-blur-[16px]"
            >
              <div className="pc-sheet-h">
                <SheetTitle>Navigation</SheetTitle>
                <SheetClose asChild>
                  <button type="button" className="clay clay-outline icon" aria-label="Close menu">
                    <X className="pc-i" aria-hidden="true" />
                  </button>
                </SheetClose>
              </div>
              <Navigation isMobile closeMenu={() => setMenuOpen(false)} />
              <div className="pc-sheet-theme">
                <span id="theme-label-sheet">Dark mode</span>
                <ThemeSwitch aria-labelledby="theme-label-sheet" />
              </div>
            </SheetContent>
          </Sheet>
        </div>
      </div>
    </header>
  );
};

export default Header;
