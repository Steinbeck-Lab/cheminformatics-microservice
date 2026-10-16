import { useEffect, useRef, useState, type ComponentType, type KeyboardEvent } from "react";
import { NavLink, useNavigate, useParams } from "react-router-dom";
import { ListTree } from "lucide-react";
import { ToolPageHead } from "./Bento";

export interface SectionTool {
  id: string;
  name: string;
  description: string;
  component: ComponentType;
  /** Sidebar layout only: group heading and icon. */
  category?: string;
  icon?: ComponentType<{ className?: string }>;
}

type ViewMode = "categories" | "all";

/** Glass sidebar: the tools grouped by category or as one list, with a show/hide button on phones. */
function Sidebar({
  label,
  basePath,
  tools,
}: {
  label: string;
  basePath: string;
  tools: SectionTool[];
}) {
  const [mode, setMode] = useState<ViewMode>("categories");
  const [open, setOpen] = useState(false);

  const groups =
    mode === "all"
      ? [{ name: null as string | null, items: tools }]
      : [...new Set(tools.map((t) => t.category ?? ""))].map((name) => ({
          name,
          items: tools.filter((t) => (t.category ?? "") === name),
        }));

  const onSegKey = (e: KeyboardEvent) => {
    if (!["ArrowLeft", "ArrowRight", "ArrowUp", "ArrowDown"].includes(e.key)) return;
    e.preventDefault();
    const next = mode === "categories" ? "all" : "categories";
    setMode(next);
    e.currentTarget.querySelector<HTMLElement>(`[data-v="${next}"]`)?.focus();
  };

  return (
    <aside className="pc-side-wrap">
      <div className="pc-side glass">
        <button
          type="button"
          className="clay pc-side-toggle"
          aria-expanded={open}
          onClick={() => setOpen((o) => !o)}
        >
          <ListTree className="pc-i" aria-hidden="true" />
          {open ? "Hide tools" : "Show tools"}
        </button>
        <div className={open ? "pc-side-body" : "pc-side-body pc-side-closed"}>
          <div
            className="pc-seg full"
            role="radiogroup"
            aria-label="Tool list view"
            onKeyDown={onSegKey}
          >
            {(["categories", "all"] as const).map((v) => (
              <button
                key={v}
                type="button"
                role="radio"
                data-v={v}
                aria-checked={mode === v}
                tabIndex={mode === v ? 0 : -1}
                onClick={() => setMode(v)}
              >
                {v === "categories" ? "Categories" : "All tools"}
              </button>
            ))}
          </div>
          <nav aria-label={label}>
            {groups.map((g) => (
              <div key={g.name ?? "all"} className="pc-side-group">
                {g.name && <h2 className="pc-side-h">{g.name}</h2>}
                {g.items.map(({ id, name, icon: Icon }) => (
                  <NavLink key={id} to={`${basePath}/${id}`} onClick={() => setOpen(false)}>
                    {Icon && <Icon className="pc-i" aria-hidden="true" />}
                    {name}
                  </NavLink>
                ))}
              </div>
            ))}
          </nav>
        </div>
      </div>
    </aside>
  );
}

/**
 * Frame for a section with several tools: a glass tab bar (saffron underline on the open tab)
 * or, with layout="sidebar", a glass sidebar; then a glass title cell and the open tool.
 * The open tool comes from the `:toolId` route param; a missing or unknown id goes to the first tool.
 */
export function SectionPage({
  label,
  basePath,
  tools,
  layout = "tabs",
}: {
  /** Accessible name of the tool navigation, e.g. "Depict tools". */
  label: string;
  basePath: string;
  tools: SectionTool[];
  layout?: "tabs" | "sidebar";
}) {
  const { toolId } = useParams<{ toolId?: string }>();
  const navigate = useNavigate();
  const tabsRef = useRef<HTMLDivElement>(null);
  const active = tools.find((t) => t.id === toolId);

  useEffect(() => {
    if (!active && tools[0]) navigate(`${basePath}/${tools[0].id}`, { replace: true });
  }, [active, basePath, navigate, tools]);

  // Keep the open tab visible when the bar scrolls sideways (many tools, narrow screens).
  useEffect(() => {
    tabsRef.current
      ?.querySelector('[aria-current="page"]')
      ?.scrollIntoView?.({ block: "nearest", inline: "nearest" });
  }, [toolId]);

  if (!active) return null;
  const Tool = active.component;
  const body = (
    <>
      <ToolPageHead title={active.name} description={active.description} />
      <Tool />
    </>
  );

  if (layout === "sidebar") {
    return (
      <div className="pc-wrap pc-side-layout">
        <Sidebar label={label} basePath={basePath} tools={tools} />
        <div className="min-w-0">{body}</div>
      </div>
    );
  }

  return (
    <div className="pc-wrap">
      <nav className="pc-subnav" aria-label={label}>
        <div ref={tabsRef} className="pc-tabs glass">
          {tools.map((t) => (
            <NavLink key={t.id} to={`${basePath}/${t.id}`}>
              {t.name}
            </NavLink>
          ))}
        </div>
      </nav>
      {body}
    </div>
  );
}
