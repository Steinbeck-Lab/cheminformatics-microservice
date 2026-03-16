// src/App.js
import { lazy, Suspense } from "react";
import {
  createBrowserRouter,
  RouterProvider,
  createRoutesFromElements,
  Route,
  useLocation,
} from "react-router-dom";
import { AppProvider } from "./context/AppContext";
import { ComparisonProvider } from "./context/ComparisonContext";
import Header from "./components/common/Header";
import Footer from "./components/common/Footer";
import { ComparisonTray } from "./components/common/ComparisonTray";
import { ComparisonView } from "./components/common/ComparisonView";
import { RouteLoadingFallback } from "./components/feedback/RouteLoadingFallback";
import { AnimatedOutlet } from "./components/common/AnimatedOutlet";
import { Toaster } from "./components/ui/sonner";
import { CommandPalette } from "./components/common/CommandPalette";
import { GradientMesh } from "./components/common/GradientMesh";
import type { GradientPageKey } from "./config/gradients";

// Lazy-loaded pages (route-level code splitting)
const HomePage = lazy(() => import("./pages/HomePage"));
const ChemPage = lazy(() => import("./pages/ChemPage"));
const ConvertPage = lazy(() => import("./pages/ConvertPage"));
const DepictPage = lazy(() => import("./pages/DepictPage"));
const ToolsPage = lazy(() => import("./pages/ToolsPage"));
const OCSRPage = lazy(() => import("./pages/OCSRPage"));
const AboutPage = lazy(() => import("./pages/AboutPage"));
const TermsOfService = lazy(() => import("./pages/TermsOfService"));
const PrivacyPolicy = lazy(() => import("./pages/PrivacyPolicy"));

// Map route segments to gradient page keys
const ROUTE_TO_GRADIENT: Record<string, GradientPageKey> = {
  "": "depict",
  home: "home",
  chem: "chem",
  convert: "convert",
  depict: "depict",
  tools: "tools",
  ocsr: "ocsr",
  about: "about",
  terms: "terms",
  privacy: "privacy",
};

// Layout component with header/footer
const Layout = () => {
  const location = useLocation();
  const segment = location.pathname.split("/").filter(Boolean)[0] || "";
  const gradientPage = ROUTE_TO_GRADIENT[segment] || "home";

  return (
    <div className="flex flex-col min-h-screen bg-background text-foreground relative">
      <GradientMesh page={gradientPage} />
      <Header />
      <main className="grow pt-20">
        <Suspense fallback={<RouteLoadingFallback />}>
          <AnimatedOutlet />
        </Suspense>
      </main>
      <Footer />
      <ComparisonTray />
      <ComparisonView />
      <Toaster />
      <CommandPalette />
    </div>
  );
};

// NotFound component
const NotFoundPage = () => (
  <div className="max-w-7xl mx-auto px-4 py-16 text-center">
    <h1 className="text-4xl font-bold text-white mb-4">Page Not Found</h1>
    <p className="text-xl text-gray-400 mb-8">
      The page you're looking for doesn't exist or has been moved.
    </p>
    <a
      href="/"
      className="px-6 py-3 bg-blue-600 hover:bg-blue-700 text-white font-medium rounded-lg shadow-lg transition-colors"
    >
      Go Home
    </a>
  </div>
);

// Create router with future flags
const router = createBrowserRouter(
  createRoutesFromElements(
    <Route path="/" element={<Layout />}>
      <Route index element={<DepictPage />} />
      <Route path="/home" element={<HomePage />} />
      <Route path="/chem" element={<ChemPage />} />
      <Route path="/chem/:toolId" element={<ChemPage />} />
      <Route path="/convert" element={<ConvertPage />} />
      <Route path="/convert/:convertId" element={<ConvertPage />} />
      <Route path="/depict" element={<DepictPage />} />
      <Route path="/depict/:depictId" element={<DepictPage />} />
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
      <ComparisonProvider>
        <RouterProvider
          router={router}
          future={{
            v7_startTransition: true,
          }}
        />
      </ComparisonProvider>
    </AppProvider>
  );
}

export default App;
