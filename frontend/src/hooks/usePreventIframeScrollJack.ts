import { useEffect, type RefObject } from "react";

/**
 * Prevents the browser from scrolling the parent document when an iframe's
 * content calls focus() during initialization.
 *
 * When a same-origin iframe's content (e.g., Ketcher editor) focuses an
 * element, the browser scrolls the parent document to bring the iframe
 * into view. On pages where the iframe is below the fold, this causes a
 * jarring scroll jump.
 *
 * This hook saves the scroll position when the iframe loads and restores
 * it if a non-user-initiated scroll occurs within the initialization window.
 * User-initiated scrolls (wheel, touch, keyboard) are never blocked.
 */
export function usePreventIframeScrollJack(iframeRef: RefObject<HTMLIFrameElement | null>) {
  useEffect(() => {
    const iframe = iframeRef.current;
    if (!iframe) return;

    let savedY = 0;
    let guarding = false;
    let userScrolled = false;
    let timer: ReturnType<typeof setTimeout>;

    const markUserScroll = () => {
      userScrolled = true;
    };

    const startGuard = () => {
      savedY = window.scrollY;
      guarding = true;
      userScrolled = false;

      window.addEventListener("wheel", markUserScroll, {
        passive: true,
        once: true,
      });
      window.addEventListener("touchmove", markUserScroll, {
        passive: true,
        once: true,
      });
      window.addEventListener("keydown", markUserScroll, { once: true });

      // Stop guarding after Ketcher's init window (~3s after iframe load)
      timer = setTimeout(stopGuard, 5000);
    };

    const stopGuard = () => {
      guarding = false;
      window.removeEventListener("wheel", markUserScroll);
      window.removeEventListener("touchmove", markUserScroll);
      window.removeEventListener("keydown", markUserScroll);
      clearTimeout(timer);
    };

    const onScroll = () => {
      if (guarding && !userScrolled) {
        window.scrollTo({ top: savedY, behavior: "instant" as ScrollBehavior });
        stopGuard();
      }
    };

    iframe.addEventListener("load", startGuard);
    window.addEventListener("scroll", onScroll);

    return () => {
      iframe.removeEventListener("load", startGuard);
      window.removeEventListener("scroll", onScroll);
      stopGuard();
    };
  }, [iframeRef]);
}
