import { useEffect, useRef } from "react";
import * as THREE from "three";
import { RoomEnvironment } from "three/examples/jsm/environments/RoomEnvironment.js";
import { OrbitControls } from "three/examples/jsm/controls/OrbitControls.js";
import { generate3DCoordinates } from "@/services/convertService";
import { assembly, parseMolblock, type Mol3 } from "@/lib/molblock";

// Ball-and-stick radii in ångström-ish scene units.
const ATOM_R: Record<string, number> = { H: 0.2, C: 0.34, N: 0.33, O: 0.33, S: 0.4, P: 0.4 };
const BOND_R = 0.085;

/** Element colours from the site tokens, so light and dark mode both match the page. */
function elementColour(el: string): string {
  const css = getComputedStyle(document.documentElement);
  const tok = (name: string, fallback: string) => css.getPropertyValue(name).trim() || fallback;
  switch (el) {
    case "C":
      return tok("--head", "#4a3d5e");
    case "N":
      return tok("--primary", "#2b2f7f");
    case "O":
      return tok("--madder", "#b23a65");
    case "S":
    case "P":
      return tok("--saffron", "#f2b31a");
    case "H":
      return tok("--el-h", "#e8e4ea");
    default:
      return tok("--mauve", "#5b4a72");
  }
}
const bondColour = () =>
  getComputedStyle(document.documentElement).getPropertyValue("--bar4").trim() || "#a39c95";

/** Seeded scatter: where each atom starts before it flies in. */
function scatter(i: number, r: number): THREE.Vector3 {
  const a = Math.sin(i * 12.9898) * 43758.5453;
  const u = a - Math.floor(a);
  const b = Math.sin(i * 78.233) * 12345.6789;
  const v = b - Math.floor(b);
  const theta = 2 * Math.PI * u;
  const phi = Math.acos(2 * v - 1);
  return new THREE.Vector3(
    Math.sin(phi) * Math.cos(theta),
    Math.sin(phi) * Math.sin(theta),
    Math.cos(phi)
  ).multiplyScalar(r);
}

/**
 * A molecule in 3D with three.js, in the site palette: glossy clear-coated atoms, silver bonds.
 * It assembles when shown (atoms fly in, then bonds grow), then turns slowly; drag to turn.
 * Under reduced motion it appears assembled and still. Without WebGL it renders nothing visual.
 */
export default function Molecule3D({
  smiles,
  label,
  className,
}: {
  smiles: string;
  label?: string;
  className?: string;
}) {
  const host = useRef<HTMLDivElement>(null);
  const api = useRef<{ show(m: Mol3): void } | null>(null);

  // the scene: built once
  useEffect(() => {
    const el = host.current;
    if (!el) return;
    let gl: THREE.WebGLRenderer;
    try {
      gl = new THREE.WebGLRenderer({ alpha: true, antialias: true });
    } catch {
      return; // no WebGL: the labelled container stays, empty
    }
    // read every frame (and the loop woken on change), so turning reduced motion on or off
    // applies at once
    const motionPref = window.matchMedia?.("(prefers-reduced-motion: reduce)");
    const reduce = () => motionPref?.matches ?? false;
    gl.setPixelRatio(Math.min(window.devicePixelRatio || 1, 2));
    gl.toneMapping = THREE.ACESFilmicToneMapping;
    gl.outputColorSpace = THREE.SRGBColorSpace;
    el.appendChild(gl.domElement);

    const scene = new THREE.Scene();
    const pmrem = new THREE.PMREMGenerator(gl);
    const env = pmrem.fromScene(new RoomEnvironment(), 0.04).texture;
    scene.environment = env;
    const key = new THREE.DirectionalLight(0xffffff, 1.4);
    key.position.set(3, 4, 5);
    const rim = new THREE.DirectionalLight(0xf2b31a, 0.9); // saffron rim from behind
    rim.position.set(-4, 2, -5);
    scene.add(key, rim, new THREE.AmbientLight(0xffffff, 0.35));

    const camera = new THREE.PerspectiveCamera(35, 1, 0.1, 200);
    const controls = new OrbitControls(camera, gl.domElement);
    controls.enableZoom = false;
    controls.enablePan = false;
    controls.enableDamping = true;
    controls.autoRotateSpeed = 1.1;

    const group = new THREE.Group();
    scene.add(group);
    const sphere = new THREE.SphereGeometry(1, 40, 28);
    const cylinder = new THREE.CylinderGeometry(1, 1, 1, 20, 1);
    const materials = new Map<string, THREE.MeshPhysicalMaterial>();
    const material = (colourOf: () => string, keyName: string) => {
      let m = materials.get(keyName);
      if (!m) {
        m = new THREE.MeshPhysicalMaterial({
          color: colourOf(),
          roughness: 0.3,
          metalness: 0.05,
          clearcoat: 1,
          clearcoatRoughness: 0.12,
        });
        materials.set(keyName, m);
      }
      return m;
    };
    const recolour = () => {
      for (const [k, m] of materials) m.color.set(k === "bond" ? bondColour() : elementColour(k));
      wake();
    };
    const themeWatch = new MutationObserver(recolour);
    themeWatch.observe(document.documentElement, { attributes: true, attributeFilter: ["class"] });

    // per-molecule state
    let atoms: { mesh: THREE.Mesh; from: THREE.Vector3; to: THREE.Vector3; r: number }[] = [];
    let bonds: {
      meshes: THREE.Mesh[];
      a: THREE.Vector3;
      dir: THREE.Vector3;
      len: number;
      offsets: THREE.Vector3[];
    }[] = [];
    let started = 0;
    // Render on demand: the loop runs while the molecule assembles, turns on its own or eases
    // out of a drag, and stops when idle, off-screen or in a hidden tab.
    let raf = 0;
    let onScreen = true;
    const up = new THREE.Vector3(0, 1, 0);

    const show = (mol: Mol3) => {
      group.clear();
      atoms = [];
      bonds = [];
      const radius = Math.max(1.5, ...mol.atoms.map((a) => Math.hypot(...a.pos))) + 0.6;
      // 1.35: the view is a circle, so the molecule must clear the round edge, not only the square
      camera.position.set(0, 0, (1.35 * radius) / Math.sin(THREE.MathUtils.degToRad(17.5)));
      camera.lookAt(0, 0, 0);
      controls.update();
      mol.atoms.forEach((a, i) => {
        const r = ATOM_R[a.el] ?? 0.38;
        const mesh = new THREE.Mesh(
          sphere,
          material(() => elementColour(a.el), a.el)
        );
        group.add(mesh);
        atoms.push({
          mesh,
          from: scatter(i + 1, radius * 2.6),
          to: new THREE.Vector3(...a.pos),
          r,
        });
      });
      for (const bd of mol.bonds) {
        const a = atoms[bd.a]!.to;
        const b = atoms[bd.b]!.to;
        const dir = new THREE.Vector3().subVectors(b, a).normalize();
        let side = new THREE.Vector3().crossVectors(dir, new THREE.Vector3(0, 0, 1));
        if (side.lengthSq() < 1e-6)
          side = new THREE.Vector3().crossVectors(dir, new THREE.Vector3(1, 0, 0));
        side.normalize().multiplyScalar(0.11);
        const offsets = bd.order === 2 ? [side, side.clone().negate()] : [new THREE.Vector3()];
        const meshes = offsets.map(() => {
          const m = new THREE.Mesh(cylinder, material(bondColour, "bond"));
          m.quaternion.setFromUnitVectors(up, dir);
          group.add(m);
          return m;
        });
        bonds.push({ meshes, a, dir, len: a.distanceTo(b), offsets });
      }
      started = performance.now();
      wake();
    };
    api.current = { show };

    const fit = () => {
      const w = el.clientWidth || 1;
      const h = el.clientHeight || 1;
      gl.setSize(w, h, false);
      gl.domElement.style.width = "100%";
      gl.domElement.style.height = "100%";
      camera.aspect = w / h;
      camera.updateProjectionMatrix();
      wake(); // resizing clears the canvas
    };
    fit();
    const resize = new ResizeObserver(fit);
    resize.observe(el);

    function wake() {
      if (!raf && onScreen && !document.hidden) raf = requestAnimationFrame(frame);
    }
    const sleep = () => {
      cancelAnimationFrame(raf);
      raf = 0;
    };
    function frame() {
      raf = 0;
      if (document.hidden || !onScreen) return;
      const still = reduce();
      const t = still ? Infinity : (performance.now() - started) / 1000;
      const s = assembly(t, { atoms: atoms.length, bonds: bonds.length });
      atoms.forEach((x, i) => {
        const p = s.atom(i);
        x.mesh.position.lerpVectors(x.from, x.to, p);
        x.mesh.scale.setScalar(Math.max(0.0001, x.r * p));
      });
      bonds.forEach((bd, j) => {
        const p = s.bond(j);
        bd.meshes.forEach((m, k) => {
          // the bond grows out of its first atom
          m.position
            .copy(bd.a)
            .addScaledVector(bd.dir, (bd.len * p) / 2)
            .add(bd.offsets[k]!);
          m.scale.set(BOND_R, Math.max(0.0001, bd.len * p), BOND_R);
        });
      });
      controls.autoRotate = s.finished && !still;
      // a moving camera fires "change", which wakes the next frame (damping settles that way)
      controls.update();
      gl.render(scene, camera);
      if (!s.finished || controls.autoRotate) wake();
    }
    // a drag, its damped ease-out, a tab coming back or the motion setting flipping: draw again
    controls.addEventListener("start", wake);
    controls.addEventListener("change", wake);
    const onVisibility = () => (document.hidden ? sleep() : wake());
    document.addEventListener("visibilitychange", onVisibility);
    motionPref?.addEventListener?.("change", wake);
    const seen =
      typeof IntersectionObserver === "undefined"
        ? null
        : new IntersectionObserver(([entry]) => {
            onScreen = entry?.isIntersecting ?? true;
            if (onScreen) wake();
            else sleep();
          });
    seen?.observe(el);
    wake();

    return () => {
      sleep();
      seen?.disconnect();
      controls.removeEventListener("start", wake);
      controls.removeEventListener("change", wake);
      document.removeEventListener("visibilitychange", onVisibility);
      motionPref?.removeEventListener?.("change", wake);
      resize.disconnect();
      themeWatch.disconnect();
      controls.dispose();
      sphere.dispose();
      cylinder.dispose();
      materials.forEach((m) => m.dispose());
      env.dispose();
      pmrem.dispose();
      gl.dispose();
      gl.domElement.remove();
      api.current = null;
    };
  }, []);

  // the molecule: fetched and assembled whenever the SMILES changes
  useEffect(() => {
    let alive = true;
    generate3DCoordinates(smiles, "rdkit")
      .then((block) => {
        const mol = parseMolblock(String(block));
        if (alive && mol.atoms.length) api.current?.show(mol);
      })
      .catch(() => {});
    return () => {
      alive = false;
    };
  }, [smiles]);

  return (
    <div
      ref={host}
      className={className}
      role={label ? "img" : undefined}
      aria-label={label}
      aria-hidden={label ? undefined : true}
    />
  );
}
