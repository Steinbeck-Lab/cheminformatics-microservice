/**
 * CaffeineMolecule3D — Immersive wireframe caffeine molecule background.
 *
 * Renders a glowing wireframe caffeine (C8H10N4O2) molecule using
 * @react-three/fiber. The molecule slowly auto-rotates and follows
 * the mouse cursor with a smooth parallax tilt effect.
 *
 * Element colors: C=cyan, N=indigo, O=rose — matching the site palette.
 * Lazy-loaded; mobile devices get a CSS-only fallback gradient.
 */
import React, { useRef, useMemo, useEffect, useState } from "react";
import { Canvas, useFrame, useThree } from "@react-three/fiber";
import { Float, Line } from "@react-three/drei";
import * as THREE from "three";

// ---------------------------------------------------------------------------
// Caffeine molecule data — heavy atoms only for clean wireframe aesthetic
// Coordinates approximate the real purine ring geometry with slight Z offsets
// ---------------------------------------------------------------------------
type AtomDef = { element: "C" | "N" | "O"; pos: [number, number, number] };

const SCALE = 1.6;

const ATOMS: AtomDef[] = [
  // Pyrimidine ring (6-membered)
  { element: "N", pos: [0, 1.4, 0.05] }, // N1
  { element: "C", pos: [1.21, 0.7, -0.03] }, // C2
  { element: "N", pos: [1.21, -0.7, 0.02] }, // N3
  { element: "C", pos: [0, -1.4, -0.05] }, // C4
  { element: "C", pos: [-1.21, -0.7, 0.03] }, // C5
  { element: "C", pos: [-1.21, 0.7, -0.02] }, // C6
  // Imidazole ring (5-membered, fused at C4–C5)
  { element: "N", pos: [-2.15, -1.55, 0.08] }, // N7
  { element: "C", pos: [-1.5, -2.75, -0.04] }, // C8
  { element: "N", pos: [-0.15, -2.65, 0.06] }, // N9
  // Carbonyl oxygens
  { element: "O", pos: [2.3, 1.25, 0.12] }, // O10 (=O on C2)
  { element: "O", pos: [-2.3, 1.25, -0.1] }, // O11 (=O on C6)
  // Methyl carbons
  { element: "C", pos: [0.45, 2.85, -0.08] }, // C12 (N1-CH3)
  { element: "C", pos: [2.4, -1.35, 0.07] }, // C13 (N3-CH3)
  { element: "C", pos: [-3.4, -1.05, -0.06] }, // C14 (N7-CH3)
].map((a) => ({
  ...a,
  pos: [a.pos[0] * SCALE, a.pos[1] * SCALE, a.pos[2] * SCALE * 3] as [number, number, number],
}));

// Bond pairs (indices into ATOMS)
const BONDS: [number, number][] = [
  // Pyrimidine ring
  [0, 1],
  [1, 2],
  [2, 3],
  [3, 4],
  [4, 5],
  [5, 0],
  // Imidazole ring
  [3, 8],
  [8, 7],
  [7, 6],
  [6, 4],
  // Carbonyls
  [1, 9],
  [5, 10],
  // Methyl
  [0, 11],
  [2, 12],
  [6, 13],
];

// Element → color mapping (site palette)
const ELEMENT_COLORS: Record<string, string> = {
  C: "#06b6d4", // cyan-500
  N: "#818cf8", // indigo-400
  O: "#f472b6", // pink-400
};

const ELEMENT_EMISSIVE: Record<string, string> = {
  C: "#0891b2", // cyan-600
  N: "#6366f1", // indigo-500
  O: "#ec4899", // pink-500
};

// ---------------------------------------------------------------------------
// Atom sphere component
// ---------------------------------------------------------------------------
const Atom = ({ position, element }: { position: [number, number, number]; element: string }) => {
  const color = ELEMENT_COLORS[element] || "#06b6d4";
  const emissive = ELEMENT_EMISSIVE[element] || "#0891b2";
  const radius = element === "O" ? 0.18 : element === "N" ? 0.16 : 0.14;

  return (
    <mesh position={position}>
      <sphereGeometry args={[radius, 16, 16]} />
      <meshStandardMaterial
        color={color}
        emissive={emissive}
        emissiveIntensity={0.8}
        transparent
        opacity={0.9}
        roughness={0.3}
        metalness={0.1}
      />
    </mesh>
  );
};

// ---------------------------------------------------------------------------
// Bond line component
// ---------------------------------------------------------------------------
const Bond = ({
  start,
  end,
}: {
  start: [number, number, number];
  end: [number, number, number];
}) => {
  return <Line points={[start, end]} color="#38bdf8" lineWidth={1.5} transparent opacity={0.5} />;
};

// ---------------------------------------------------------------------------
// Glow ring — decorative orbit ring around the molecule
// ---------------------------------------------------------------------------
const GlowRing = ({
  radius,
  color,
  rotationSpeed,
}: {
  radius: number;
  color: string;
  rotationSpeed: number;
}) => {
  const ref = useRef<THREE.Mesh>(null!);
  useFrame((_, delta) => {
    ref.current.rotation.z += delta * rotationSpeed;
  });
  return (
    <mesh ref={ref} rotation={[Math.PI / 2, 0, 0]}>
      <torusGeometry args={[radius, 0.015, 8, 100]} />
      <meshStandardMaterial
        color={color}
        emissive={color}
        emissiveIntensity={0.6}
        transparent
        opacity={0.2}
      />
    </mesh>
  );
};

// ---------------------------------------------------------------------------
// Molecule group — handles rotation + mouse parallax
// ---------------------------------------------------------------------------
const MoleculeScene = () => {
  const groupRef = useRef<THREE.Group>(null!);
  const mouseTarget = useRef({ x: 0, y: 0 });
  const { size } = useThree();

  useEffect(() => {
    const handleMouseMove = (e: MouseEvent) => {
      mouseTarget.current.x = (e.clientX / size.width - 0.5) * 2;
      mouseTarget.current.y = -(e.clientY / size.height - 0.5) * 2;
    };
    window.addEventListener("mousemove", handleMouseMove);
    return () => window.removeEventListener("mousemove", handleMouseMove);
  }, [size]);

  useFrame((state, delta) => {
    if (!groupRef.current) return;

    // Slow auto-rotation
    groupRef.current.rotation.y += delta * 0.08;
    groupRef.current.rotation.x += delta * 0.03;

    // Smooth parallax toward mouse position
    const targetRotX = mouseTarget.current.y * 0.3;
    const targetRotY = mouseTarget.current.x * 0.4;

    groupRef.current.rotation.x += (targetRotX - groupRef.current.rotation.x) * 0.02;
    groupRef.current.rotation.y += (targetRotY - groupRef.current.rotation.y) * 0.02;

    // Subtle position shift for depth
    groupRef.current.position.x +=
      (mouseTarget.current.x * 0.5 - groupRef.current.position.x) * 0.02;
    groupRef.current.position.y +=
      (mouseTarget.current.y * 0.3 - groupRef.current.position.y) * 0.02;
  });

  // Particle field — small floating dots for atmosphere
  const particles = useMemo(() => {
    const points: [number, number, number][] = [];
    for (let i = 0; i < 60; i++) {
      points.push([
        (Math.random() - 0.5) * 14,
        (Math.random() - 0.5) * 14,
        (Math.random() - 0.5) * 8,
      ]);
    }
    return points;
  }, []);

  return (
    <>
      {/* Ambient + directional lighting */}
      <ambientLight intensity={0.15} />
      <directionalLight position={[5, 5, 5]} intensity={0.4} color="#e0f2fe" />
      <pointLight position={[-3, 2, 4]} intensity={0.6} color="#06b6d4" distance={20} />
      <pointLight position={[3, -2, -4]} intensity={0.3} color="#818cf8" distance={15} />

      <Float speed={0.8} rotationIntensity={0.1} floatIntensity={0.3}>
        <group ref={groupRef}>
          {/* Molecule atoms */}
          {ATOMS.map((atom, i) => (
            <Atom key={`atom-${i}`} position={atom.pos} element={atom.element} />
          ))}

          {/* Molecule bonds */}
          {BONDS.map(([a, b], i) => (
            <Bond key={`bond-${i}`} start={ATOMS[a].pos} end={ATOMS[b].pos} />
          ))}

          {/* Decorative glow rings */}
          <GlowRing radius={5.5} color="#06b6d4" rotationSpeed={0.15} />
          <GlowRing radius={4.2} color="#818cf8" rotationSpeed={-0.1} />
        </group>
      </Float>

      {/* Floating particles */}
      {particles.map((pos, i) => (
        <mesh key={`particle-${i}`} position={pos}>
          <sphereGeometry args={[0.03 + Math.random() * 0.03, 6, 6]} />
          <meshStandardMaterial
            color="#38bdf8"
            emissive="#0ea5e9"
            emissiveIntensity={0.5}
            transparent
            opacity={0.3 + Math.random() * 0.3}
          />
        </mesh>
      ))}
    </>
  );
};

// ---------------------------------------------------------------------------
// Reduced motion check
// ---------------------------------------------------------------------------
const useReducedMotion = () => {
  const [reduced, setReduced] = useState(false);
  useEffect(() => {
    const mq = window.matchMedia("(prefers-reduced-motion: reduce)");
    setReduced(mq.matches);
    const handler = (e: MediaQueryListEvent) => setReduced(e.matches);
    mq.addEventListener("change", handler);
    return () => mq.removeEventListener("change", handler);
  }, []);
  return reduced;
};

// ---------------------------------------------------------------------------
// Mobile detection (skip WebGL on small screens)
// ---------------------------------------------------------------------------
const useIsMobile = () => {
  const [mobile, setMobile] = useState(false);
  useEffect(() => {
    const check = () => setMobile(window.innerWidth < 768);
    check();
    window.addEventListener("resize", check);
    return () => window.removeEventListener("resize", check);
  }, []);
  return mobile;
};

// ---------------------------------------------------------------------------
// Exported component — full-viewport 3D canvas background
// ---------------------------------------------------------------------------
export const CaffeineMolecule3D: React.FC<{ className?: string }> = ({ className = "" }) => {
  const isMobile = useIsMobile();
  const reducedMotion = useReducedMotion();

  // Mobile or reduced motion: skip WebGL, show nothing (GradientMesh handles bg)
  if (isMobile || reducedMotion) return null;

  return (
    <div
      className={`absolute inset-0 z-0 pointer-events-none ${className}`}
      style={{ opacity: 0.6 }}
    >
      <Canvas
        dpr={[1, 1.5]}
        camera={{ position: [0, 0, 10], fov: 50 }}
        style={{ background: "transparent" }}
        gl={{ alpha: true, antialias: true, powerPreference: "high-performance" }}
      >
        <MoleculeScene />
      </Canvas>
    </div>
  );
};

export default CaffeineMolecule3D;
