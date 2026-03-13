import { PAGE_GRADIENTS } from "../../config/gradients";
import type { GradientPageKey } from "../../config/gradients";
import { useAppContext } from "../../context/AppContext";
import { useMediaQuery } from "../../hooks/useMediaQuery";

interface GradientMeshProps {
  /** The page key determining which gradient palette to use. */
  page: GradientPageKey;
}

/**
 * Reusable gradient mesh background component.
 *
 * Renders an absolutely-positioned div with the per-page gradient mesh
 * background. Uses the app's dark mode context and a media query to
 * select the appropriate gradient variant (light/dark, desktop/mobile).
 *
 * Place this as a child of a relatively-positioned container.
 */
export function GradientMesh({ page }: GradientMeshProps) {
  const { isDarkMode } = useAppContext();
  const isMobile = useMediaQuery("(max-width: 767px)");

  const gradients = PAGE_GRADIENTS[page];
  let background: string;

  if (isMobile) {
    background = isDarkMode ? gradients.mobileDark : gradients.mobileLight;
  } else {
    background = isDarkMode ? gradients.dark : gradients.light;
  }

  return (
    <div className="absolute inset-0 -z-20 overflow-hidden" aria-hidden="true">
      <div className="absolute inset-0" style={{ background }} />
    </div>
  );
}
