import SugarRemovalView from "../components/tools/SugarRemovalView";
import StructureGenView from "../components/tools/StructureGenView";
import InChIView from "../components/tools/InChIView";
import RInChIView from "../components/tools/RInChIView";
import { SectionPage, type SectionTool } from "../components/common/SectionPage";

const tools: SectionTool[] = [
  {
    id: "sugardetection",
    name: "Sugar Detection",
    description: "Detect and remove sugar moieties from complex molecular structures.",
    component: SugarRemovalView,
  },
  {
    id: "structuregeneration",
    name: "Structure Generation",
    description: "Generate every structure for a molecular formula with Surge.",
    component: StructureGenView,
  },
  {
    id: "inchiconverter",
    name: "InChI Converter",
    description: "Draw, edit and convert structures to InChI, with InChI version options.",
    component: InChIView,
  },
  {
    id: "rinchiconverter",
    name: "RInChI Converter",
    description: "Convert chemical reactions to RInChI notation.",
    component: RInChIView,
  },
];

const ToolsPage = () => <SectionPage label="Tools" basePath="/tools" tools={tools} />;

export default ToolsPage;
