/**
 * CaffeineMolecule3D — Neon-glow wireframe caffeine molecule background.
 *
 * Uses @react-three/postprocessing Bloom to create a real neon glow effect.
 * Atoms are emissive spheres (intensity > 1 triggers bloom). Bonds are
 * glowing tubes. A BufferGeometry point cloud adds atmospheric depth.
 *
 * The molecule auto-rotates and tracks the mouse with smooth parallax.
 * Lazy-loaded; mobile/reduced-motion gets CSS fallback.
 */
import React, { useRef, useMemo, useEffect, useState } from "react";
import { Canvas, useFrame, useThree } from "@react-three/fiber";
import { Float } from "@react-three/drei";
import { EffectComposer, Bloom } from "@react-three/postprocessing";
import * as THREE from "three";

// ---------------------------------------------------------------------------
// Caffeine molecule data (C8H10N4O2) — heavy atoms, approximate 3D coords
// ---------------------------------------------------------------------------
type AtomDef = { element: "C" | "N" | "O"; pos: [number, number, number] };

const SCALE = 1.8;

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
  pos: [a.pos[0] * SCALE, a.pos[1] * SCALE, a.pos[2] * SCALE * 4] as [number, number, number],
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

// Colors — high-intensity for bloom trigger
const ELEMENT_COLOR: Record<string, [number, number, number]> = {
  C: [0.1, 0.85, 0.95], // cyan
  N: [0.45, 0.4, 1.0], // indigo
  O: [1.0, 0.35, 0.65], // rose
};

// ---------------------------------------------------------------------------
// Glowing atom — emissiveIntensity > 1 + toneMapped=false triggers bloom
// ---------------------------------------------------------------------------
const Atom = ({
  position,
  element,
  pulse,
}: {
  position: [number, number, number];
  element: string;
  pulse: number;
}) => {
  const [r, g, b] = ELEMENT_COLOR[element] || ELEMENT_COLOR.C;
  const radius = element === "O" ? 0.22 : element === "N" ? 0.19 : 0.16;
  const intensity = 2.5 + Math.sin(pulse) * 0.8;

  return (
    <mesh position={position}>
      <sphereGeometry args={[radius, 24, 24]} />
      <meshBasicMaterial color={[r * intensity, g * intensity, b * intensity]} toneMapped={false} />
    </mesh>
  );
};

// ---------------------------------------------------------------------------
// Glowing bond — cylinder tube between two atoms
// ---------------------------------------------------------------------------
const Bond = ({
  start,
  end,
  pulse,
}: {
  start: [number, number, number];
  end: [number, number, number];
  pulse: number;
}) => {
  const ref = useRef<THREE.Mesh>(null!);

  useMemo(() => {
    if (!ref.current) return;
    const s = new THREE.Vector3(...start);
    const e = new THREE.Vector3(...end);
    const mid = s.clone().add(e).multiplyScalar(0.5);
    const dir = e.clone().sub(s);
    const len = dir.length();
    ref.current.position.copy(mid);
    ref.current.scale.set(1, len, 1);
    ref.current.quaternion.setFromUnitVectors(new THREE.Vector3(0, 1, 0), dir.normalize());
  }, [start, end]);

  const s = new THREE.Vector3(...start);
  const e = new THREE.Vector3(...end);
  const mid = s.clone().add(e).multiplyScalar(0.5);
  const dir = e.clone().sub(s);
  const len = dir.length();
  const quat = new THREE.Quaternion().setFromUnitVectors(
    new THREE.Vector3(0, 1, 0),
    dir.normalize()
  );

  const intensity = 1.5 + Math.sin(pulse * 0.7) * 0.3;

  return (
    <mesh ref={ref} position={[mid.x, mid.y, mid.z]} quaternion={quat}>
      <cylinderGeometry args={[0.04, 0.04, len, 8, 1]} />
      <meshBasicMaterial
        color={[0.2 * intensity, 0.75 * intensity, 1.0 * intensity]}
        toneMapped={false}
        transparent
        opacity={0.7}
      />
    </mesh>
  );
};

// ---------------------------------------------------------------------------
// Particle field — BufferGeometry Points for atmosphere
// ---------------------------------------------------------------------------
const ParticleField = () => {
  const ref = useRef<THREE.Points>(null!);
  const count = 200;

  const [positions, sizes] = useMemo(() => {
    const pos = new Float32Array(count * 3);
    const sz = new Float32Array(count);
    for (let i = 0; i < count; i++) {
      pos[i * 3] = (Math.random() - 0.5) * 20;
      pos[i * 3 + 1] = (Math.random() - 0.5) * 20;
      pos[i * 3 + 2] = (Math.random() - 0.5) * 12;
      sz[i] = Math.random() * 3 + 1;
    }
    return [pos, sz];
  }, []);

  useFrame((state) => {
    if (!ref.current) return;
    const geo = ref.current.geometry;
    const posAttr = geo.attributes.position as THREE.BufferAttribute;
    const t = state.clock.elapsedTime * 0.15;
    for (let i = 0; i < count; i++) {
      posAttr.array[i * 3 + 1] += Math.sin(t + i * 0.5) * 0.002;
    }
    posAttr.needsUpdate = true;
  });

  return (
    <points ref={ref}>
      <bufferGeometry>
        <bufferAttribute attach="attributes-position" args={[positions, 3]} />
        <bufferAttribute attach="attributes-size" args={[sizes, 1]} />
      </bufferGeometry>
      <pointsMaterial
        size={0.06}
        color={[1.5, 2.5, 4.0]}
        toneMapped={false}
        transparent
        opacity={0.6}
        sizeAttenuation
        depthWrite={false}
        blending={THREE.AdditiveBlending}
      />
    </points>
  );
};

// ---------------------------------------------------------------------------
// Orbit ring — glowing torus that slowly spins
// ---------------------------------------------------------------------------
const OrbitRing = ({
  radius,
  color,
  speed,
  tilt,
}: {
  radius: number;
  color: [number, number, number];
  speed: number;
  tilt: [number, number, number];
}) => {
  const ref = useRef<THREE.Mesh>(null!);
  useFrame((_, delta) => {
    ref.current.rotation.z += delta * speed;
  });
  return (
    <mesh ref={ref} rotation={tilt}>
      <torusGeometry args={[radius, 0.02, 16, 120]} />
      <meshBasicMaterial color={color} toneMapped={false} transparent opacity={0.35} />
    </mesh>
  );
};

// ---------------------------------------------------------------------------
// Main scene — molecule + effects + mouse tracking
// ---------------------------------------------------------------------------
const MoleculeScene = () => {
  const groupRef = useRef<THREE.Group>(null!);
  const mouse = useRef({ x: 0, y: 0 });
  const { size } = useThree();
  const pulseRef = useRef(0);

  useEffect(() => {
    const onMove = (e: MouseEvent) => {
      mouse.current.x = (e.clientX / window.innerWidth - 0.5) * 2;
      mouse.current.y = -(e.clientY / window.innerHeight - 0.5) * 2;
    };
    window.addEventListener("mousemove", onMove);
    return () => window.removeEventListener("mousemove", onMove);
  }, [size]);

  useFrame((state, delta) => {
    if (!groupRef.current) return;
    pulseRef.current = state.clock.elapsedTime;

    // Slow auto-rotation
    groupRef.current.rotation.y += delta * 0.06;
    groupRef.current.rotation.x += delta * 0.02;

    // Mouse parallax (smooth lerp)
    const tx = mouse.current.x * 0.35;
    const ty = mouse.current.y * 0.25;
    groupRef.current.rotation.x += (ty - groupRef.current.rotation.x) * 0.015;
    groupRef.current.rotation.y += (tx - groupRef.current.rotation.y) * 0.015;
    groupRef.current.position.x += (mouse.current.x * 0.6 - groupRef.current.position.x) * 0.015;
    groupRef.current.position.y += (mouse.current.y * 0.4 - groupRef.current.position.y) * 0.015;
  });

  return (
    <>
      <color attach="background" args={["#000000"]} />
      <fog attach="fog" args={["#000308", 8, 25]} />

      <Float speed={0.6} rotationIntensity={0.08} floatIntensity={0.2}>
        <group ref={groupRef}>
          {ATOMS.map((atom, i) => (
            <Atom
              key={`a${i}`}
              position={atom.pos}
              element={atom.element}
              pulse={pulseRef.current + i * 0.4}
            />
          ))}
          {BONDS.map(([a, b], i) => (
            <Bond key={`b${i}`} start={ATOMS[a].pos} end={ATOMS[b].pos} pulse={pulseRef.current} />
          ))}
          <OrbitRing radius={6} color={[0.3, 0.8, 1.2]} speed={0.12} tilt={[Math.PI / 2.5, 0, 0]} />
          <OrbitRing
            radius={5}
            color={[0.5, 0.4, 1.4]}
            speed={-0.08}
            tilt={[Math.PI / 3, 0.3, 0]}
          />
          <OrbitRing
            radius={7.5}
            color={[0.2, 0.6, 0.9]}
            speed={0.05}
            tilt={[Math.PI / 2, -0.2, 0.4]}
          />
        </group>
      </Float>

      <ParticleField />

      {/* Bloom postprocessing — the key to neon glow */}
      <EffectComposer>
        <Bloom luminanceThreshold={0.8} luminanceSmoothing={0.3} intensity={1.8} mipmapBlur />
      </EffectComposer>
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

// ---------------------------------------------------------------------------
// Export — full-viewport 3D canvas
// ---------------------------------------------------------------------------
export const CaffeineMolecule3D: React.FC<{ className?: string }> = ({ className = "" }) => {
  const isMobile = useIsMobile();
  const reduced = useReducedMotion();

  if (isMobile || reduced) return null;

  return (
    <div className={`absolute inset-0 z-0 pointer-events-none ${className}`}>
      <Canvas
        dpr={[1, 1.5]}
        camera={{ position: [0, 0, 12], fov: 45 }}
        gl={{
          alpha: false,
          antialias: true,
          powerPreference: "high-performance",
          toneMapping: THREE.ACESFilmicToneMapping,
          toneMappingExposure: 1.2,
        }}
      >
        <MoleculeScene />
      </Canvas>
    </div>
  );
};

export default CaffeineMolecule3D;
