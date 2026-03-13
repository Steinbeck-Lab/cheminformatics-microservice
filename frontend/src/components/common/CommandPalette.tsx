/**
 * CommandPalette -- Cmd+K / Ctrl+K global command palette.
 *
 * Searches across pages, tools, example molecules, and recent molecules.
 * Uses cmdk (via shadcn Command) for fuzzy search and keyboard navigation.
 */
import React, { useEffect, useState } from "react";
import { useNavigate } from "react-router-dom";
import {
  CommandDialog,
  CommandInput,
  CommandList,
  CommandEmpty,
  CommandGroup,
  CommandItem,
  CommandSeparator,
} from "@/components/ui/command";
import { pages, tools, exampleMolecules } from "@/data/navigationRegistry";
import { useAppContext } from "../../context/AppContext";
import { Clock, FlaskConical } from "lucide-react";

export function CommandPalette() {
  const [open, setOpen] = useState(false);
  const navigate = useNavigate();
  const { recentMolecules } = useAppContext();

  // Global keyboard shortcut: Cmd+K (Mac) or Ctrl+K (Windows)
  useEffect(() => {
    function handleKeyDown(e: KeyboardEvent) {
      if (e.key === "k" && (e.metaKey || e.ctrlKey)) {
        e.preventDefault();
        setOpen((prev) => !prev);
      }
    }
    document.addEventListener("keydown", handleKeyDown);
    return () => document.removeEventListener("keydown", handleKeyDown);
  }, []);

  function handleSelect(path: string) {
    navigate(path);
    setOpen(false);
  }

  function handleMoleculeSelect(smiles: string) {
    navigate("/chem?smiles=" + encodeURIComponent(smiles));
    setOpen(false);
  }

  return (
    <CommandDialog
      open={open}
      onOpenChange={setOpen}
      title="Command Palette"
      description="Search pages, tools, and molecules"
    >
      <CommandInput placeholder="Search pages, tools, molecules..." />
      <CommandList>
        <CommandEmpty>No results found.</CommandEmpty>

        <CommandGroup heading="Pages">
          {pages.map((page) => (
            <CommandItem key={page.id} value={page.name} onSelect={() => handleSelect(page.path)}>
              {page.icon && <page.icon className="mr-2 h-4 w-4" />}
              {page.name}
            </CommandItem>
          ))}
        </CommandGroup>

        <CommandSeparator />

        <CommandGroup heading="Tools">
          {tools.map((tool) => (
            <CommandItem
              key={tool.id}
              value={`${tool.name} ${tool.keywords?.join(" ") || ""}`}
              onSelect={() => handleSelect(tool.path)}
            >
              {tool.icon && <tool.icon className="mr-2 h-4 w-4" />}
              <span>{tool.name}</span>
              <span className="ml-auto text-xs text-muted-foreground">{tool.section}</span>
            </CommandItem>
          ))}
        </CommandGroup>

        <CommandSeparator />

        <CommandGroup heading="Example Molecules">
          {exampleMolecules.map((mol) => (
            <CommandItem
              key={mol.name}
              value={`${mol.name} ${mol.smiles}`}
              onSelect={() => handleMoleculeSelect(mol.smiles)}
            >
              <FlaskConical className="mr-2 h-4 w-4" />
              <span>{mol.name}</span>
              <span className="ml-auto text-xs text-muted-foreground font-mono truncate max-w-[200px]">
                {mol.smiles}
              </span>
            </CommandItem>
          ))}
        </CommandGroup>

        {recentMolecules.length > 0 && (
          <>
            <CommandSeparator />
            <CommandGroup heading="Recent">
              {recentMolecules.map((mol, i) => (
                <CommandItem
                  key={`recent-${i}`}
                  value={`${mol.name || "Recent"} ${mol.smiles}`}
                  onSelect={() => handleMoleculeSelect(mol.smiles)}
                >
                  <Clock className="mr-2 h-4 w-4" />
                  <span>{mol.name || "Recent Molecule"}</span>
                  <span className="ml-auto text-xs text-muted-foreground font-mono truncate max-w-[200px]">
                    {mol.smiles}
                  </span>
                </CommandItem>
              ))}
            </CommandGroup>
          </>
        )}
      </CommandList>
    </CommandDialog>
  );
}
