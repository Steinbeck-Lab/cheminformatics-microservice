// src/App.js
import React from "react";
import {
  createBrowserRouter,
  RouterProvider,
  createRoutesFromElements,
  Route,
  Outlet,
} from "react-router-dom";
import { AppProvider } from "./context/AppContext";
import Header from "./components/common/Header";
import Footer from "./components/common/Footer";
import HomePage from "./pages/HomePage";
import ChemPage from "./pages/ChemPage";
import ConvertPage from "./pages/ConvertPage";
import DepictPage from "./pages/DepictPage";
import ToolsPage from "./pages/ToolsPage";
import OCSRPage from "./pages/OCSRPage";
import AboutPage from "./pages/AboutPage";
import TermsOfService from "./pages/TermsOfService";
import PrivacyPolicy from "./pages/PrivacyPolicy";

// Layout component with header/footer
const Layout = () => (
  <div className="flex flex-col min-h-screen bg-gray-900 text-white">
    <Header />
    <main className="flex-grow">
      <Outlet />
    </main>
    <Footer />
  </div>
);

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
      <Route index element={<HomePage />} />
      <Route path="/chem" element={<ChemPage />} />
      <Route path="/convert" element={<ConvertPage />} />
      <Route path="/depict" element={<DepictPage />} />
      <Route path="/tools" element={<ToolsPage />} />
      <Route path="/ocsr" element={<OCSRPage />} />
      <Route path="/about" element={<AboutPage />} />
      <Route path="/terms" element={<TermsOfService />} />
      <Route path="/privacy" element={<PrivacyPolicy />} />
      <Route path="*" element={<NotFoundPage />} />
    </Route>
  ),
  {
    future: {
      v7_startTransition: true,
      v7_relativeSplatPath: true,
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
