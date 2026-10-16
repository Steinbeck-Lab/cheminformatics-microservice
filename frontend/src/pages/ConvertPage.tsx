import FormatConversionView from "../components/convert/FormatConversionView";
import Mol2DView from "../components/convert/Mol2DView";
import Mol3DView from "../components/convert/Mol3DView";
import { SectionPage, type SectionTool } from "../components/common/SectionPage";

const tools: SectionTool[] = [
  {
    id: "format-conversion",
    name: "Format Conversion",
    description: "Convert between chemical identifiers and file formats.",
    component: FormatConversionView,
  },
  {
    id: "2d-coordinates",
    name: "2D Coordinates",
    description: "Generate 2D coordinates for molecular structures.",
    component: Mol2DView,
  },
  {
    id: "3d-coordinates",
    name: "3D Coordinates",
    description: "Generate 3D coordinates for molecular structures.",
    component: Mol3DView,
  },
];

const ConvertPage = () => <SectionPage label="Convert tools" basePath="/convert" tools={tools} />;

export default ConvertPage;
