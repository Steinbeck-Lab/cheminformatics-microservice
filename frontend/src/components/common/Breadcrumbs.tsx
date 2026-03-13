/**
 * Breadcrumbs -- contextual section > tool navigation for tool pages.
 *
 * Renders below the fixed header on /chem, /convert, /depict, /tools routes.
 * Hidden on non-tool pages (home, about, ocsr, privacy, terms) and on mobile (<768px).
 * Tool name cross-fades at ~150ms when switching tabs.
 */
import React from "react";
import { useLocation, useParams, Link } from "react-router-dom";
import { AnimatePresence, motion } from "motion/react";
import { ChevronRight } from "lucide-react";
import { tools, getSectionByPath } from "@/data/navigationRegistry";

// Route param names vary by page
interface RouteParams {
  toolId?: string;
  convertId?: string;
  depictId?: string;
}

// Tool pages where breadcrumbs should appear
const TOOL_ROUTES = ["/chem", "/convert", "/depict", "/tools"];

export function Breadcrumbs() {
  const { pathname } = useLocation();
  const params = useParams<Record<string, string>>();

  // Only show on tool page routes
  const basePath = "/" + pathname.split("/").filter(Boolean)[0];
  if (!TOOL_ROUTES.includes(basePath)) {
    return null;
  }

  const section = getSectionByPath(pathname);
  if (!section) return null;

  // Resolve the tool ID from whichever route param is active
  const toolId = params.toolId || params.convertId || params.depictId;

  // Look up the tool name from the navigation registry
  const tool = toolId ? tools.find((t) => t.id === toolId) : undefined;
  const toolName = tool?.name;

  return (
    <nav
      className="hidden md:block max-w-7xl mx-auto px-4 sm:px-6 lg:px-8 pt-4 pb-0"
      aria-label="Breadcrumb"
      data-testid="breadcrumbs"
    >
      <ol className="flex items-center gap-1.5 text-sm text-muted-foreground">
        <li>
          <Link to={section.path} className="hover:text-foreground transition-colors font-medium">
            {section.name}
          </Link>
        </li>
        {toolName && (
          <>
            <li>
              <ChevronRight className="h-3.5 w-3.5" />
            </li>
            <li>
              <AnimatePresence mode="wait">
                <motion.span
                  key={toolName}
                  initial={{ opacity: 0 }}
                  animate={{ opacity: 1 }}
                  exit={{ opacity: 0 }}
                  transition={{ duration: 0.15 }}
                  className="text-foreground font-medium"
                >
                  {toolName}
                </motion.span>
              </AnimatePresence>
            </li>
          </>
        )}
      </ol>
    </nav>
  );
}
