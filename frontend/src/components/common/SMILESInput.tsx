// Description: A React component for inputting SMILES notation with example and recent molecule dropdowns.
import { useState, useRef } from "react";
// Import necessary icons
import { useAppContext } from "../../context/AppContext"; // Assuming context provides recentMolecules
import { Clipboard, Clock, X } from "lucide-react";
import { Button } from "@/components/ui/button";
import { Textarea } from "@/components/ui/textarea";
import { Dialog, DialogContent, DialogTitle } from "@/components/ui/dialog";
import type { RecentMolecule } from "@/types/molecule";

// Example molecules data (static, moved outside component)
const EXAMPLE_MOLECULES = [
  {
    name: "Caffeine",
    smiles: "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
    description: "Stimulant found in coffee",
  },
  {
    name: "Aspirin",
    smiles: "CC(=O)OC1=CC=CC=C1C(=O)O",
    description: "Analgesic & anti-inflammatory",
  },
  {
    name: "Sucrose",
    smiles: "C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O[C@]2([C@H]([C@@H]([C@H](O2)CO)O)O)CO)O)O)O)O",
    description: "Common sugar",
  },
  {
    name: "Cholesterol",
    smiles: "C[C@H](CCCC(C)C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC=C4[C@@]3(CC[C@@H](C4)O)C)C",
    description: "Steroid in cell membranes",
  },
  {
    name: "Paracetamol",
    smiles: "CC(=O)NC1=CC=C(C=C1)O",
    description: "Pain reliever & fever reducer",
  },
  {
    name: "Ibuprofen",
    smiles: "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O",
    description: "Anti-inflammatory drug",
  },
];

// Dropdown List Item Component (for cleaner code)

interface PasteModalProps {
  isOpen: boolean;
  onClose: () => void;
  onPaste: (text: string) => void;
}

// Paste Modal Component
const PasteModal = ({ isOpen, onClose, onPaste }: PasteModalProps) => {
  const [text, setText] = useState("");
  const textareaRef = useRef<HTMLTextAreaElement>(null);

  return (
    <Dialog open={isOpen} onOpenChange={(open) => !open && onClose()}>
      <DialogContent
        showCloseButton={false}
        className="block gap-0 border-0 bg-card p-6 rounded-xl shadow-2xl max-w-lg w-full"
        onOpenAutoFocus={(e) => {
          e.preventDefault();
          textareaRef.current?.focus();
        }}
      >
        <div className="flex justify-between items-center mb-4">
          <DialogTitle className="text-lg font-bold text-foreground flex items-center">
            <Clipboard className="h-5 w-5 mr-2 text-primary " />
            Paste SMILES
          </DialogTitle>
          <Button
            variant="ghost"
            onClick={onClose}

            aria-label="Close"
          >
            <X className="w-5 h-5" />
          </Button>
        </div>
        <p className="mb-4 text-foreground ">Please paste your SMILES string in the area below:</p>
        <div className="mb-4">
          <Textarea
            ref={textareaRef}
            value={text}
            onChange={(e) => setText(e.target.value)}
            className="w-full font-mono text-sm h-24"
            placeholder="Paste SMILES here..."
            aria-label="SMILES string"
          />
        </div>
        <div className="flex justify-end gap-3">
          <Button variant="secondary" onClick={onClose}>
            Cancel
          </Button>
          <Button
            onClick={() => {
              onPaste(text);
              onClose();
            }}

            disabled={!text.trim()}
          >
            Use SMILES
          </Button>
        </div>
      </DialogContent>
    </Dialog>
  );
};

interface SMILESInputProps {
  value: string;
  onChange: (smiles: string) => void;
  placeholder?: string;
  showExamples?: boolean;
  showRecent?: boolean;
  label?: string;
  required?: boolean;
  id?: string;
}

const SMILESInput = ({
  value,
  onChange,
  placeholder = "Enter a SMILES string",
  showExamples = true,
  showRecent = true,
  label = "SMILES",
  required = false,
  id: inputId = "smiles-input",
}: SMILESInputProps) => {
  const [showPasteModal, setShowPasteModal] = useState(false);
  const { recentMolecules = [] } = useAppContext() || {};

  // Paste: the clipboard API first, then the legacy command, then a dialog to paste into by hand.
  const handlePaste = async () => {
    if (navigator.clipboard?.readText) {
      try {
        const text = await navigator.clipboard.readText();
        if (text) return onChange(text.trim());
      } catch {
        // permission refused: try the next way
      }
    }
    try {
      const area = Object.assign(document.createElement("textarea"), {
        style: "position:fixed;top:-999px;left:-999px",
      });
      document.body.appendChild(area);
      area.focus();
      document.execCommand("paste");
      const pasted = area.value;
      area.remove();
      if (pasted) return onChange(pasted.trim());
    } catch {
      // fall through to the dialog
    }
    setShowPasteModal(true);
  };

  // A molecule chip: an example or a recent molecule.
  const chip = (
    m: Pick<RecentMolecule, "smiles" | "name" | "description">,
    key: string | number
  ) => (
    <button
      key={key}
      type="button"
      className="clay clay-outline pc-chip"
      aria-pressed={value === m.smiles}
      title={m.description || m.smiles}
      onClick={() => onChange(m.smiles)}
    >
      {m.name || (m.smiles.length > 14 ? `${m.smiles.slice(0, 14)}…` : m.smiles)}
    </button>
  );

  return (
    <div className="pc-smiles">
      {label && (
        <label htmlFor={inputId} className="pc-lbl">
          {label} {required && <span className="text-madder-ink">*</span>}
        </label>
      )}
      <div className="pc-smiles-field">
        <input
          id={inputId}
          type="text"
          className="pc-well"
          value={value}
          onChange={(e) => onChange(e.target.value)}
          placeholder={placeholder}
          spellCheck={false}
          autoComplete="off"
          autoCapitalize="off"
          required={required}
          aria-required={required}
        />
        <Button variant="outline" type="button" onClick={handlePaste}>
          <Clipboard />
          Paste
        </Button>
        {value && (
          <Button
            variant="destructive"
            size="icon"
            type="button"
            onClick={() => onChange("")}
            aria-label="Clear input"
          >
            <X />
          </Button>
        )}
      </div>
      {showExamples && (
        <div className="pc-smiles-chips">
          <span className="pc-lbl2" id={`${inputId}-ex`}>
            Examples
          </span>
          <div className="pc-row" role="group" aria-labelledby={`${inputId}-ex`}>
            {EXAMPLE_MOLECULES.map((m) => chip(m, m.name))}
          </div>
        </div>
      )}
      {showRecent && recentMolecules.length > 0 && (
        <div className="pc-smiles-chips">
          <span className="pc-lbl2" id={`${inputId}-recent`}>
            <Clock className="pc-i" aria-hidden="true" /> Recent
          </span>
          <div className="pc-row" role="group" aria-labelledby={`${inputId}-recent`}>
            {recentMolecules.slice(0, 4).map((m, i) => chip(m, m.timestamp || i))}
          </div>
        </div>
      )}
      <p className="pc-note">
        Need a SMILES? Find one in the{" "}
        <a href="/depict/structureexplorer" className="text-link">
          Structure Explorer
        </a>
        .
      </p>
      <PasteModal
        isOpen={showPasteModal}
        onClose={() => setShowPasteModal(false)}
        onPaste={(text) => text?.trim() && onChange(text.trim())}
      />
    </div>
  );
};
export default SMILESInput;
