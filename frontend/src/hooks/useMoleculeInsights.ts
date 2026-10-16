import { useEffect, useState } from "react";
import { calculateDescriptors } from "@/services/chemService";
import { generateInChI, generateInChIKey } from "@/services/convertService";

export interface Insights {
  smiles: string;
  /** RDKit descriptors, including molecular_formula and nplikeness. */
  descriptors?: Record<string, string | number | boolean>;
  inchi?: string;
  inchikey?: string;
  failed: { descriptors: boolean; ids: boolean };
}

/** The convert endpoints answer with a JSON string, quotes included. */
const unquote = (v: unknown) => String(v).replace(/^"([\s\S]*)"$/, "$1");

/**
 * Identifiers, NP-likeness and descriptors for one structure (RDKit), fetched in parallel.
 * Returns null without a SMILES; data for a previous SMILES is never returned for a new one.
 */
export function useMoleculeInsights(smiles: string | undefined): Insights | null {
  const [insights, setInsights] = useState<Insights | null>(null);

  useEffect(() => {
    if (!smiles) return setInsights(null);
    let alive = true;
    setInsights({ smiles, failed: { descriptors: false, ids: false } });
    Promise.allSettled([
      calculateDescriptors(smiles, "rdkit"),
      generateInChI(smiles, "rdkit"),
      generateInChIKey(smiles, "rdkit"),
    ]).then(([d, i, k]) => {
      if (!alive) return;
      setInsights({
        smiles,
        descriptors: d.status === "fulfilled" ? d.value : undefined,
        inchi: i.status === "fulfilled" ? unquote(i.value) : undefined,
        inchikey: k.status === "fulfilled" ? unquote(k.value) : undefined,
        failed: {
          descriptors: d.status === "rejected",
          ids: i.status === "rejected" && k.status === "rejected",
        },
      });
    });
    return () => {
      alive = false;
    };
  }, [smiles]);

  return insights?.smiles === smiles ? insights : null;
}
