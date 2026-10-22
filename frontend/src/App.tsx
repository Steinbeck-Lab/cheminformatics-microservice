// src/App.js
import { lazy, useEffect, type ComponentType } from "react";
import {
  createBrowserRouter,
  RouterProvider,
  createRoutesFromElements,
  Navigate,
  Route,
  useLocation,
} from "react-router-dom";
import { AppProvider } from "./context/AppContext";
import { clearReloadFlag, reloadOnce } from "./utils/reloadOnce";
import Header from "./components/common/Header";
import Footer from "./components/common/Footer";
import { AnimatedOutlet } from "./components/common/AnimatedOutlet";
import { Toaster } from "./components/ui/sonner";
import { CommandPalette } from "./components/common/CommandPalette";
import { SilkBackground } from "./components/common/SilkBackground";

// Lazy-loaded pages (route-level code splitting). Once a page's code has arrived, render it
// directly: React.lazy would still suspend for a moment and hold its placeholder for ~300 ms.
type PageModule = { default: ComponentType };
const loaders: (() => Promise<PageModule>)[] = [];
function lazyPage(loader: () => Promise<PageModule>) {
  let Page: ComponentType | null = null;
  const load = () =>
    loader().then((m) => {
      Page = m.default;
      clearReloadFlag();
      return m;
    });
  loaders.push(load);
  // only the real render path reloads; the idle prefetch swallows its own failures
  const Lazy = lazy(() => load().catch(reloadOnce));
  return function LazyPage() {
    return Page ? <Page /> : <Lazy />;
  };
}
const HomePage = lazyPage(() => import("./pages/HomePage"));
const ChemPage = lazyPage(() => import("./pages/ChemPage"));
const ConvertPage = lazyPage(() => import("./pages/ConvertPage"));
const DepictPage = lazyPage(() => import("./pages/DepictPage"));
const ToolsPage = lazyPage(() => import("./pages/ToolsPage"));
const OCSRPage = lazyPage(() => import("./pages/OCSRPage"));
const AboutPage = lazyPage(() => import("./pages/AboutPage"));
const TermsOfService = lazyPage(() => import("./pages/TermsOfService"));
const PrivacyPolicy = lazyPage(() => import("./pages/PrivacyPolicy"));

// Once the first page is up, fetch every section's code while the browser is idle, so switching
// sections later is instant instead of showing the loading placeholder.
const prefetchPages = () => loaders.forEach((l) => void l().catch(() => {}));

// Layout component with header/footer
const Layout = () => {
  const { pathname } = useLocation();
  const liveSilk = pathname === "/home" || pathname === "/about";
  useEffect(() => {
    // a visitor who asked the browser to save data gets each section only when they open it
    if ((navigator as Navigator & { connection?: { saveData?: boolean } }).connection?.saveData)
      return;
    const idle = window.requestIdleCallback ?? ((cb: () => void) => window.setTimeout(cb, 1500));
    idle(prefetchPages);
  }, []);
  return (
    <div className="flex flex-col min-h-screen bg-background text-foreground relative isolate">
      <a
        href="#main-content"
        className="sr-only focus:not-sr-only focus:fixed focus:top-3 focus:left-3 focus:z-100 focus:rounded-lg focus:bg-background focus:px-4 focus:py-2 focus:text-foreground focus:shadow-lg focus:outline-2 focus:outline-ring"
      >
        Skip to content
      </a>
      <SilkBackground live={liveSilk} />
      <Header />
      <main id="main-content" tabIndex={-1} className="grow focus:outline-none">
        <AnimatedOutlet />
      </main>
      <Footer />
      <Toaster />
      <CommandPalette />
    </div>
  );
};

// NotFound component
const NotFoundPage = () => (
  <div className="pc-wrap pc-notfound">
    <section className="glass pc-notfound-card" aria-labelledby="nf-h">
      <h1 id="nf-h" className="pc-display">
        Page not found
      </h1>
      <p>This page does not exist or has moved.</p>
      <a className="clay clay-primary" href="/home">
        Go home
      </a>
    </section>
  </div>
);

// Create router with future flags
const router = createBrowserRouter(
  createRoutesFromElements(
    <Route path="/" element={<Layout />}>
      <Route index element={<Navigate to="/home" replace />} />
      <Route path="/home" element={<HomePage />} />
      <Route path="/chem" element={<ChemPage />} />
      <Route path="/chem/:toolId" element={<ChemPage />} />
      <Route path="/convert" element={<ConvertPage />} />
      <Route path="/convert/:toolId" element={<ConvertPage />} />
      <Route path="/depict" element={<DepictPage />} />
      <Route path="/depict/:toolId" element={<DepictPage />} />
      <Route path="/tools" element={<ToolsPage />} />
      <Route path="/tools/:toolId" element={<ToolsPage />} />
      <Route path="/ocsr" element={<OCSRPage />} />
      <Route path="/about" element={<AboutPage />} />
      <Route path="/terms" element={<TermsOfService />} />
      <Route path="/privacy" element={<PrivacyPolicy />} />
      <Route path="*" element={<NotFoundPage />} />
    </Route>
  ),
  {
    future: {
      v7_fetcherPersist: true,
      v7_normalizeFormMethod: true,
      v7_partialHydration: true,
      v7_relativeSplatPath: true,
      v7_skipActionErrorRevalidation: true,
    },
  }
);

function App() {
  return (
    <AppProvider>
      <RouterProvider router={router} />
    </AppProvider>
  );
}

export default App;
