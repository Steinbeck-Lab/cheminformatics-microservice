/**
 * Per-page gradient mesh definitions for background atmospheres.
 *
 * Each page has a unique color palette visible through glass panels.
 * Light mode uses soft pastels, dark mode uses deep jewel tones.
 * Mobile variants use simplified 2-color linear gradients for GPU performance.
 */
export const PAGE_GRADIENTS = {
  home: {
    light: [
      "radial-gradient(ellipse at 20% 50%, rgba(147,197,253,0.5) 0%, transparent 50%)",
      "radial-gradient(ellipse at 80% 20%, rgba(196,181,253,0.4) 0%, transparent 50%)",
      "radial-gradient(ellipse at 50% 80%, rgba(165,180,252,0.3) 0%, transparent 50%)",
    ].join(", "),
    dark: [
      "radial-gradient(ellipse at 20% 50%, rgba(30,64,175,0.3) 0%, transparent 50%)",
      "radial-gradient(ellipse at 80% 20%, rgba(88,28,135,0.25) 0%, transparent 50%)",
      "radial-gradient(ellipse at 50% 80%, rgba(49,46,129,0.2) 0%, transparent 50%)",
    ].join(", "),
    mobileLight: "linear-gradient(135deg, rgba(147,197,253,0.3) 0%, rgba(196,181,253,0.2) 100%)",
    mobileDark: "linear-gradient(135deg, rgba(30,64,175,0.2) 0%, rgba(88,28,135,0.15) 100%)",
  },
  chem: {
    light: [
      "radial-gradient(ellipse at 30% 40%, rgba(125,211,252,0.5) 0%, transparent 50%)",
      "radial-gradient(ellipse at 70% 70%, rgba(94,234,212,0.4) 0%, transparent 50%)",
      "radial-gradient(ellipse at 50% 10%, rgba(186,230,253,0.3) 0%, transparent 50%)",
    ].join(", "),
    dark: [
      "radial-gradient(ellipse at 30% 40%, rgba(14,116,144,0.3) 0%, transparent 50%)",
      "radial-gradient(ellipse at 70% 70%, rgba(13,148,136,0.25) 0%, transparent 50%)",
      "radial-gradient(ellipse at 50% 10%, rgba(8,145,178,0.2) 0%, transparent 50%)",
    ].join(", "),
    mobileLight: "linear-gradient(135deg, rgba(125,211,252,0.3) 0%, rgba(94,234,212,0.2) 100%)",
    mobileDark: "linear-gradient(135deg, rgba(14,116,144,0.2) 0%, rgba(13,148,136,0.15) 100%)",
  },
  convert: {
    light: [
      "radial-gradient(ellipse at 25% 60%, rgba(165,180,252,0.5) 0%, transparent 50%)",
      "radial-gradient(ellipse at 75% 30%, rgba(147,197,253,0.4) 0%, transparent 50%)",
      "radial-gradient(ellipse at 50% 90%, rgba(129,140,248,0.3) 0%, transparent 50%)",
    ].join(", "),
    dark: [
      "radial-gradient(ellipse at 25% 60%, rgba(49,46,129,0.3) 0%, transparent 50%)",
      "radial-gradient(ellipse at 75% 30%, rgba(30,64,175,0.25) 0%, transparent 50%)",
      "radial-gradient(ellipse at 50% 90%, rgba(67,56,202,0.2) 0%, transparent 50%)",
    ].join(", "),
    mobileLight: "linear-gradient(135deg, rgba(165,180,252,0.3) 0%, rgba(147,197,253,0.2) 100%)",
    mobileDark: "linear-gradient(135deg, rgba(49,46,129,0.2) 0%, rgba(30,64,175,0.15) 100%)",
  },
  depict: {
    light: [
      "radial-gradient(ellipse at 35% 45%, rgba(196,181,253,0.5) 0%, transparent 50%)",
      "radial-gradient(ellipse at 65% 75%, rgba(249,168,212,0.4) 0%, transparent 50%)",
      "radial-gradient(ellipse at 50% 15%, rgba(216,180,254,0.3) 0%, transparent 50%)",
    ].join(", "),
    dark: [
      "radial-gradient(ellipse at 35% 45%, rgba(88,28,135,0.3) 0%, transparent 50%)",
      "radial-gradient(ellipse at 65% 75%, rgba(157,23,77,0.25) 0%, transparent 50%)",
      "radial-gradient(ellipse at 50% 15%, rgba(109,40,217,0.2) 0%, transparent 50%)",
    ].join(", "),
    mobileLight: "linear-gradient(135deg, rgba(196,181,253,0.3) 0%, rgba(249,168,212,0.2) 100%)",
    mobileDark: "linear-gradient(135deg, rgba(88,28,135,0.2) 0%, rgba(157,23,77,0.15) 100%)",
  },
  ocsr: {
    light: [
      "radial-gradient(ellipse at 40% 35%, rgba(110,231,183,0.5) 0%, transparent 50%)",
      "radial-gradient(ellipse at 60% 65%, rgba(103,232,249,0.4) 0%, transparent 50%)",
      "radial-gradient(ellipse at 50% 85%, rgba(167,243,208,0.3) 0%, transparent 50%)",
    ].join(", "),
    dark: [
      "radial-gradient(ellipse at 40% 35%, rgba(4,120,87,0.3) 0%, transparent 50%)",
      "radial-gradient(ellipse at 60% 65%, rgba(14,116,144,0.25) 0%, transparent 50%)",
      "radial-gradient(ellipse at 50% 85%, rgba(6,95,70,0.2) 0%, transparent 50%)",
    ].join(", "),
    mobileLight: "linear-gradient(135deg, rgba(110,231,183,0.3) 0%, rgba(103,232,249,0.2) 100%)",
    mobileDark: "linear-gradient(135deg, rgba(4,120,87,0.2) 0%, rgba(14,116,144,0.15) 100%)",
  },
  tools: {
    light: [
      "radial-gradient(ellipse at 30% 55%, rgba(253,230,138,0.5) 0%, transparent 50%)",
      "radial-gradient(ellipse at 70% 25%, rgba(253,186,116,0.4) 0%, transparent 50%)",
      "radial-gradient(ellipse at 50% 80%, rgba(254,215,170,0.3) 0%, transparent 50%)",
    ].join(", "),
    dark: [
      "radial-gradient(ellipse at 30% 55%, rgba(146,64,14,0.3) 0%, transparent 50%)",
      "radial-gradient(ellipse at 70% 25%, rgba(154,52,18,0.25) 0%, transparent 50%)",
      "radial-gradient(ellipse at 50% 80%, rgba(120,53,15,0.2) 0%, transparent 50%)",
    ].join(", "),
    mobileLight: "linear-gradient(135deg, rgba(253,230,138,0.3) 0%, rgba(253,186,116,0.2) 100%)",
    mobileDark: "linear-gradient(135deg, rgba(146,64,14,0.2) 0%, rgba(154,52,18,0.15) 100%)",
  },
  about: {
    light: [
      "radial-gradient(ellipse at 20% 40%, rgba(148,163,184,0.5) 0%, transparent 50%)",
      "radial-gradient(ellipse at 80% 60%, rgba(147,197,253,0.4) 0%, transparent 50%)",
      "radial-gradient(ellipse at 50% 80%, rgba(165,180,252,0.3) 0%, transparent 50%)",
    ].join(", "),
    dark: [
      "radial-gradient(ellipse at 20% 40%, rgba(51,65,85,0.3) 0%, transparent 50%)",
      "radial-gradient(ellipse at 80% 60%, rgba(30,64,175,0.25) 0%, transparent 50%)",
      "radial-gradient(ellipse at 50% 80%, rgba(49,46,129,0.2) 0%, transparent 50%)",
    ].join(", "),
    mobileLight: "linear-gradient(135deg, rgba(148,163,184,0.3) 0%, rgba(147,197,253,0.2) 100%)",
    mobileDark: "linear-gradient(135deg, rgba(51,65,85,0.2) 0%, rgba(30,64,175,0.15) 100%)",
  },
  privacy: {
    light: [
      "radial-gradient(ellipse at 25% 45%, rgba(148,163,184,0.4) 0%, transparent 50%)",
      "radial-gradient(ellipse at 75% 55%, rgba(156,163,175,0.3) 0%, transparent 50%)",
      "radial-gradient(ellipse at 50% 80%, rgba(163,163,163,0.25) 0%, transparent 50%)",
    ].join(", "),
    dark: [
      "radial-gradient(ellipse at 25% 45%, rgba(51,65,85,0.3) 0%, transparent 50%)",
      "radial-gradient(ellipse at 75% 55%, rgba(55,65,81,0.25) 0%, transparent 50%)",
      "radial-gradient(ellipse at 50% 80%, rgba(63,63,70,0.2) 0%, transparent 50%)",
    ].join(", "),
    mobileLight: "linear-gradient(135deg, rgba(148,163,184,0.3) 0%, rgba(156,163,175,0.2) 100%)",
    mobileDark: "linear-gradient(135deg, rgba(51,65,85,0.2) 0%, rgba(55,65,81,0.15) 100%)",
  },
  terms: {
    light: [
      "radial-gradient(ellipse at 25% 45%, rgba(148,163,184,0.4) 0%, transparent 50%)",
      "radial-gradient(ellipse at 75% 55%, rgba(156,163,175,0.3) 0%, transparent 50%)",
      "radial-gradient(ellipse at 50% 80%, rgba(163,163,163,0.25) 0%, transparent 50%)",
    ].join(", "),
    dark: [
      "radial-gradient(ellipse at 25% 45%, rgba(51,65,85,0.3) 0%, transparent 50%)",
      "radial-gradient(ellipse at 75% 55%, rgba(55,65,81,0.25) 0%, transparent 50%)",
      "radial-gradient(ellipse at 50% 80%, rgba(63,63,70,0.2) 0%, transparent 50%)",
    ].join(", "),
    mobileLight: "linear-gradient(135deg, rgba(148,163,184,0.3) 0%, rgba(156,163,175,0.2) 100%)",
    mobileDark: "linear-gradient(135deg, rgba(51,65,85,0.2) 0%, rgba(55,65,81,0.15) 100%)",
  },
} as const;

/** Type for page keys in the gradient config. */
export type GradientPageKey = keyof typeof PAGE_GRADIENTS;
