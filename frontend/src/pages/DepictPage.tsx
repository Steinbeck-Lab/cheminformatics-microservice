import StructureVisualizerView from "../components/depict/StructureVisualizerView";
import Depict2DView from "../components/depict/Depict2DView";
import Depict3DView from "../components/depict/Depict3DView";
import StructureDrawView from "../components/depict/StructureDrawView";
import { SectionPage, type SectionTool } from "../components/common/SectionPage";

const tools: SectionTool[] = [
  {
    id: "structureexplorer",
    name: "Structure Explorer",
    description: "Find structures by name or identifier and visualize them in 2D and 3D.",
    component: StructureVisualizerView,
  },
  {
    id: "2ddepiction",
    name: "2D Depiction",
    description:
      "Render a structure from SMILES with CDK or RDKit, compare both toolkits, or draw many at once.",
    component: Depict2DView,
  },
  {
    id: "3ddepiction",
    name: "3D Depiction",
    description: "Create interactive 3D visualizations.",
    component: Depict3DView,
  },
  {
    id: "structuredraw",
    name: "Draw a Structure",
    description: "Draw and edit chemical structures in the Ketcher editor.",
    component: StructureDrawView,
  },
];

const DepictPage = () => <SectionPage label="Depict tools" basePath="/depict" tools={tools} />;

export default DepictPage;
