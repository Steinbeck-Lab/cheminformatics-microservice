// Optical Chemical Structure Recognition (OCSR) with DECIMER.
import { ExternalLink } from "lucide-react";
import OCRView from "../components/ocsr/OCRView";
import { ToolPageHead } from "@/components/common/Bento";
import { Cell } from "@/components/common/porcelain-ui";
import { Button } from "@/components/ui/button";

const REFERENCES = [
  {
    text: "Rajan, K., Brinkhaus, H.O., Agea, I.A., Zielesny, A., Steinbeck, C. (2023). DECIMER.ai: an open platform for automated optical chemical structure identification, segmentation and recognition in scientific publications. Nat Commun, 14, 5045.",
    doi: "10.1038/s41467-023-40782-0",
  },
  {
    text: "Rajan K, Weissenborn VK, Lederer L, Zielesny A, Steinbeck C (2025). MARCUS: Molecular annotation and recognition for curating unravelled structures. Digit Discovery.",
    doi: "10.1039/d5dd00313j",
  },
];

const OCSRPage = () => (
  <div className="pc-wrap">
    <ToolPageHead
      title="Optical Chemical Structure Recognition"
      description="Extract chemical structures from images using deep learning (DECIMER & MARCUS)."
    />
    <OCRView />
    <div className="pc-bento">
      <Cell id="ocsr-decimer" title="DECIMER" className="pc-t-input pc-cell-end">
        <p className="pc-tool-note">
          DECIMER (Deep lEarning for Chemical IMagE Recognition) is a deep learning model from the
          Steinbeck Lab that turns chemical structures in images into SMILES. It reads publications,
          patents and other documents.
        </p>
        <div className="pc-row">
          <Button asChild>
            <a href="https://decimer.ai" target="_blank" rel="noopener noreferrer">
              Visit DECIMER
              <ExternalLink />
            </a>
          </Button>
        </div>
      </Cell>
      <Cell id="ocsr-marcus" title="MARCUS" className="pc-t-side pc-cell-end">
        <p className="pc-tool-note">
          MARCUS (Molecular Annotation and Recognition for Curating Unravelled Structures) is a web
          platform for natural product literature curation. It combines text annotation, OCSR with
          DECIMER, MolNexTR and MolScribe, and database submission.
        </p>
        <div className="pc-row">
          <Button asChild variant="outline">
            <a href="https://marcus.decimer.ai" target="_blank" rel="noopener noreferrer">
              Visit MARCUS
              <ExternalLink />
            </a>
          </Button>
        </div>
      </Cell>
      <Cell id="ocsr-refs" title="References" className="pc-t-wide">
        <ol className="pc-ref-list">
          {REFERENCES.map((r) => (
            <li key={r.doi}>
              <p className="pc-ref-text">{r.text}</p>
              <p className="pc-ref-links">
                <a href={`https://doi.org/${r.doi}`} target="_blank" rel="noopener noreferrer">
                  doi.org/{r.doi}
                </a>
              </p>
            </li>
          ))}
        </ol>
      </Cell>
    </div>
  </div>
);

export default OCSRPage;
