/**
 * CaffeineMolecule3D — Minimal wireframe caffeine molecule background.
 *
 * Subtle, elegant molecule that adapts to light/dark mode.
 * Dark mode: soft cyan/indigo glow with gentle bloom.
 * Light mode: muted slate/blue wireframe, very understated.
 * No orbit rings, minimal particles — molecule is the focus.
 */
import React, { useRef, useMemo, useEffect, useState, useCallback } from "react";
import { Canvas, useFrame, useThree } from "@react-three/fiber";
import { Float } from "@react-three/drei";
import { EffectComposer, Bloom } from "@react-three/postprocessing";
import * as THREE from "three";

// ---------------------------------------------------------------------------
// Caffeine molecule data (C8H10N4O2)
// ---------------------------------------------------------------------------
type AtomDef = { element: "C" | "N" | "O"; pos: [number, number, number] };
const SCALE = 1.6;

const RAW_ATOMS: AtomDef[] = [
  { element: "N", pos: [0, 1.4, 0.05] },
  { element: "C", pos: [1.21, 0.7, -0.03] },
  { element: "N", pos: [1.21, -0.7, 0.02] },
  { element: "C", pos: [0, -1.4, -0.05] },
  { element: "C", pos: [-1.21, -0.7, 0.03] },
  { element: "C", pos: [-1.21, 0.7, -0.02] },
  { element: "N", pos: [-2.15, -1.55, 0.08] },
  { element: "C", pos: [-1.5, -2.75, -0.04] },
  { element: "N", pos: [-0.15, -2.65, 0.06] },
  { element: "O", pos: [2.3, 1.25, 0.12] },
  { element: "O", pos: [-2.3, 1.25, -0.1] },
  { element: "C", pos: [0.45, 2.85, -0.08] },
  { element: "C", pos: [2.4, -1.35, 0.07] },
  { element: "C", pos: [-3.4, -1.05, -0.06] },
];

const ATOMS = RAW_ATOMS.map((a) => ({
  ...a,
  pos: [a.pos[0] * SCALE, a.pos[1] * SCALE, a.pos[2] * SCALE * 3] as [number, number, number],
}));

const BONDS: [number, number][] = [
  [0, 1],
  [1, 2],
  [2, 3],
  [3, 4],
  [4, 5],
  [5, 0],
  [3, 8],
  [8, 7],
  [7, 6],
  [6, 4],
  [1, 9],
  [5, 10],
  [0, 11],
  [2, 12],
  [6, 13],
];

// Theme-aware colors
const DARK_COLORS: Record<string, [number, number, number]> = {
  C: [0.05, 0.55, 0.65],
  N: [0.3, 0.28, 0.7],
  O: [0.65, 0.2, 0.42],
};
const LIGHT_COLORS: Record<string, [number, number, number]> = {
  C: [0.25, 0.45, 0.55],
  N: [0.3, 0.3, 0.55],
  O: [0.55, 0.25, 0.4],
};

// ---------------------------------------------------------------------------
// Atom
// ---------------------------------------------------------------------------
const Atom = ({
  position,
  element,
  dark,
}: {
  position: [number, number, number];
  element: string;
  dark: boolean;
}) => {
  const colors = dark ? DARK_COLORS : LIGHT_COLORS;
  const [r, g, b] = colors[element] || colors.C;
  const radius = element === "O" ? 0.16 : element === "N" ? 0.14 : 0.12;
  const mult = dark ? 2.0 : 0.9;

  return (
    <mesh position={position}>
      <sphereGeometry args={[radius, 20, 20]} />
      <meshBasicMaterial
        color={[r * mult, g * mult, b * mult]}
        toneMapped={false}
        transparent
        opacity={dark ? 0.95 : 0.6}
      />
    </mesh>
  );
};

// ---------------------------------------------------------------------------
// Bond — thin cylinder
// ---------------------------------------------------------------------------
const Bond = ({
  start,
  end,
  dark,
}: {
  start: [number, number, number];
  end: [number, number, number];
  dark: boolean;
}) => {
  const s = new THREE.Vector3(...start);
  const e = new THREE.Vector3(...end);
  const mid = s.clone().add(e).multiplyScalar(0.5);
  const dir = e.clone().sub(s);
  const len = dir.length();
  const quat = new THREE.Quaternion().setFromUnitVectors(
    new THREE.Vector3(0, 1, 0),
    dir.normalize()
  );

  const c = dark ? [0.12, 0.45, 0.6] : [0.35, 0.45, 0.55];
  const mult = dark ? 1.3 : 0.7;

  return (
    <mesh position={[mid.x, mid.y, mid.z]} quaternion={quat}>
      <cylinderGeometry args={[0.025, 0.025, len, 6, 1]} />
      <meshBasicMaterial
        color={[c[0] * mult, c[1] * mult, c[2] * mult]}
        toneMapped={false}
        transparent
        opacity={dark ? 0.5 : 0.3}
      />
    </mesh>
  );
};

// ---------------------------------------------------------------------------
// Sparse particle dust
// ---------------------------------------------------------------------------
const ParticleDust = ({ dark }: { dark: boolean }) => {
  const ref = useRef<THREE.Points>(null!);
  const count = 80;

  const positions = useMemo(() => {
    const pos = new Float32Array(count * 3);
    for (let i = 0; i < count; i++) {
      pos[i * 3] = (Math.random() - 0.5) * 18;
      pos[i * 3 + 1] = (Math.random() - 0.5) * 18;
      pos[i * 3 + 2] = (Math.random() - 0.5) * 10;
    }
    return pos;
  }, []);

  useFrame((state) => {
    if (!ref.current) return;
    const posAttr = ref.current.geometry.attributes.position as THREE.BufferAttribute;
    const t = state.clock.elapsedTime * 0.1;
    for (let i = 0; i < count; i++) {
      posAttr.array[i * 3 + 1] += Math.sin(t + i * 0.3) * 0.001;
    }
    posAttr.needsUpdate = true;
  });

  const col: [number, number, number] = dark ? [0.6, 1.2, 1.8] : [0.3, 0.4, 0.55];

  return (
    <points ref={ref}>
      <bufferGeometry>
        <bufferAttribute attach="attributes-position" args={[positions, 3]} />
      </bufferGeometry>
      <pointsMaterial
        size={0.04}
        color={col}
        toneMapped={false}
        transparent
        opacity={dark ? 0.4 : 0.2}
        sizeAttenuation
        depthWrite={false}
        blending={dark ? THREE.AdditiveBlending : THREE.NormalBlending}
      />
    </points>
  );
};

// ---------------------------------------------------------------------------
// Scene
// ---------------------------------------------------------------------------
const MoleculeScene = ({ dark }: { dark: boolean }) => {
  const groupRef = useRef<THREE.Group>(null!);
  const mouse = useRef({ x: 0, y: 0 });

  useEffect(() => {
    const onMove = (e: MouseEvent) => {
      mouse.current.x = (e.clientX / window.innerWidth - 0.5) * 2;
      mouse.current.y = -(e.clientY / window.innerHeight - 0.5) * 2;
    };
    window.addEventListener("mousemove", onMove);
    return () => window.removeEventListener("mousemove", onMove);
  }, []);

  useFrame((_, delta) => {
    if (!groupRef.current) return;
    groupRef.current.rotation.y += delta * 0.05;
    groupRef.current.rotation.x += delta * 0.015;
    groupRef.current.rotation.x += (mouse.current.y * 0.2 - groupRef.current.rotation.x) * 0.01;
    groupRef.current.rotation.y += (mouse.current.x * 0.25 - groupRef.current.rotation.y) * 0.01;
    groupRef.current.position.x += (mouse.current.x * 0.4 - groupRef.current.position.x) * 0.01;
    groupRef.current.position.y += (mouse.current.y * 0.25 - groupRef.current.position.y) * 0.01;
  });

  return (
    <>
      <Float speed={0.5} rotationIntensity={0.05} floatIntensity={0.15}>
        <group ref={groupRef}>
          {ATOMS.map((atom, i) => (
            <Atom key={`a${i}`} position={atom.pos} element={atom.element} dark={dark} />
          ))}
          {BONDS.map(([a, b], i) => (
            <Bond key={`b${i}`} start={ATOMS[a].pos} end={ATOMS[b].pos} dark={dark} />
          ))}
        </group>
      </Float>

      <ParticleDust dark={dark} />

      {dark && (
        <EffectComposer>
          <Bloom luminanceThreshold={0.6} luminanceSmoothing={0.4} intensity={0.8} mipmapBlur />
        </EffectComposer>
      )}
    </>
  );
};

// ---------------------------------------------------------------------------
// Hooks
// ---------------------------------------------------------------------------
const useReducedMotion = () => {
  const [r, setR] = useState(false);
  useEffect(() => {
    const mq = window.matchMedia("(prefers-reduced-motion: reduce)");
    setR(mq.matches);
    const h = (e: MediaQueryListEvent) => setR(e.matches);
    mq.addEventListener("change", h);
    return () => mq.removeEventListener("change", h);
  }, []);
  return r;
};

const useIsMobile = () => {
  const [m, setM] = useState(false);
  useEffect(() => {
    const c = () => setM(window.innerWidth < 768);
    c();
    window.addEventListener("resize", c);
    return () => window.removeEventListener("resize", c);
  }, []);
  return m;
};

const useDarkMode = () => {
  const [dark, setDark] = useState(false);
  const check = useCallback(() => {
    setDark(document.documentElement.classList.contains("dark"));
  }, []);
  useEffect(() => {
    check();
    const obs = new MutationObserver(check);
    obs.observe(document.documentElement, { attributes: true, attributeFilter: ["class"] });
    return () => obs.disconnect();
  }, [check]);
  return dark;
};

// ---------------------------------------------------------------------------
// Export
// ---------------------------------------------------------------------------
export const CaffeineMolecule3D: React.FC<{ className?: string }> = ({ className = "" }) => {
  const isMobile = useIsMobile();
  const reduced = useReducedMotion();
  const dark = useDarkMode();

  if (isMobile || reduced) return null;

  return (
    <div
      className={`absolute inset-0 z-0 pointer-events-none ${className}`}
      style={{ opacity: dark ? 0.7 : 0.4 }}
    >
      <Canvas
        dpr={[1, 1.5]}
        camera={{ position: [0, 0, 12], fov: 45 }}
        gl={{ alpha: true, antialias: true, powerPreference: "high-performance" }}
        style={{ background: "transparent" }}
      >
        <MoleculeScene dark={dark} />
      </Canvas>
    </div>
  );
};

export default CaffeineMolecule3D;
