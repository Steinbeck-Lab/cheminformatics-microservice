/**
 * Molecule data type definitions.
 */

export interface RecentMolecule {
  smiles: string;
  name?: string;
  timestamp?: string;
  description?: string;
  _savedAt: number;
}

export interface MoleculeData {
  smiles: string;
  descriptors: Record<string, string | number | boolean>;
  npScore: number;
  timestamp: string;
}

export interface ApiConfig {
  baseUrl: string;
  timeout: number;
}

export interface AppContextValue {
  isDarkMode: boolean;
  toggleDarkMode: () => void;
  recentMolecules: RecentMolecule[];
  addRecentMolecule: (molecule: Omit<RecentMolecule, "_savedAt">) => void;
  clearRecentMolecules: () => void;
  apiConfig: ApiConfig;
  updateApiConfig: (newConfig: Partial<ApiConfig>) => void;
  isLoading: boolean;
  setIsLoading: React.Dispatch<React.SetStateAction<boolean>>;
  globalError: string | null;
  setGlobalError: React.Dispatch<React.SetStateAction<string | null>>;
}
