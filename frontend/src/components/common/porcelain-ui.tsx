import { useState, type KeyboardEvent, type ReactNode } from "react";
import { Check, Copy, RotateCw } from "lucide-react";
import { Button } from "@/components/ui/button";

// Porcelain building blocks shared by the tool pages: a bento cell, a segmented radio group,
// a switch, a labelled select, a rotation slider and a copy button. Styles live in styles/porcelain.css.

/** Clay segmented radio group with arrow-key movement. */
export function Seg<T extends string>({
  label,
  labelId,
  value,
  options,
  onChange,
  disabled,
  full,
}: {
  label?: string;
  labelId?: string;
  value: T;
  options: { value: T; label: string }[];
  onChange: (v: T) => void;
  disabled?: boolean;
  full?: boolean;
}) {
  const onKey = (e: KeyboardEvent<HTMLDivElement>) => {
    const step = { ArrowRight: 1, ArrowDown: 1, ArrowLeft: -1, ArrowUp: -1 }[e.key];
    if (!step) return;
    e.preventDefault();
    const i = options.findIndex((o) => o.value === value);
    const next = options[(i + step + options.length) % options.length];
    if (!next) return;
    onChange(next.value);
    e.currentTarget.querySelector<HTMLElement>(`[data-v="${next.value}"]`)?.focus();
  };
  return (
    <div
      className={full ? "pc-seg full" : "pc-seg"}
      role="radiogroup"
      aria-label={label}
      aria-labelledby={labelId}
      onKeyDown={onKey}
    >
      {options.map((o) => (
        <button
          key={o.value}
          type="button"
          role="radio"
          data-v={o.value}
          aria-checked={o.value === value}
          tabIndex={o.value === value ? 0 : -1}
          disabled={disabled}
          onClick={() => onChange(o.value)}
        >
          {o.label}
        </button>
      ))}
    </div>
  );
}

export function Switch({
  checked,
  onChange,
  label,
  hint,
  disabled,
}: {
  checked: boolean;
  onChange: (v: boolean) => void;
  label: string;
  hint?: string;
  disabled?: boolean;
}) {
  return (
    <button
      type="button"
      role="switch"
      aria-checked={checked}
      className="pc-switch"
      disabled={disabled}
      onClick={() => onChange(!checked)}
    >
      <span className="pc-txt">
        {label}
        {hint && <small>{hint}</small>}
      </span>
      <span className="pc-track" aria-hidden="true" />
    </button>
  );
}

export function Select({
  label,
  value,
  options,
  onChange,
  disabled,
}: {
  label: string;
  value: string;
  options: [string, string][];
  onChange: (v: string) => void;
  disabled?: boolean;
}) {
  return (
    <label className="pc-field">
      <span>{label}</span>
      <span className="pc-select">
        <select value={value} disabled={disabled} onChange={(e) => onChange(e.target.value)}>
          {options.map(([v, l]) => (
            <option key={v} value={v}>
              {l}
            </option>
          ))}
        </select>
      </span>
    </label>
  );
}

/**
 * Rotation slider, 0–359°. With an `id` the label text is visible (`<label htmlFor>`); without
 * one, `label` becomes the slider's aria-label next to the icon.
 */
export function RotationSlider({
  id,
  value,
  onChange,
  label,
  disabled,
  className,
  onCommit,
}: {
  id?: string;
  value: number;
  onChange: (deg: number) => void;
  /** Fires on pointer release, so a debounced render can settle at once. */
  onCommit?: () => void;
  label: string;
  disabled?: boolean;
  className?: string;
}) {
  const icon = <RotateCw className="pc-i" aria-hidden="true" />;
  return (
    <div className={className ? `pc-rot ${className}` : "pc-rot"}>
      {id ? (
        <label htmlFor={id}>
          {icon}
          {label}
        </label>
      ) : (
        icon
      )}
      <input
        id={id}
        type="range"
        min={0}
        max={359}
        step={1}
        value={value}
        disabled={disabled}
        aria-label={id ? undefined : label}
        onChange={(e) => onChange(Number(e.target.value))}
        onPointerUp={onCommit}
      />
      <output htmlFor={id}>{value}°</output>
    </div>
  );
}

export function Cell({
  id,
  title,
  meta,
  className,
  children,
}: {
  id: string;
  title: string;
  meta?: ReactNode;
  className: string;
  children: ReactNode;
}) {
  return (
    <section className={`pc-cell glass ${className}`} aria-labelledby={id}>
      <div className="pc-cell-h">
        <h2 id={id}>{title}</h2>
        {meta && <div className="pc-meta">{meta}</div>}
      </div>
      {children}
    </section>
  );
}

export function CopyButton({
  text,
  label,
  icon = false,
}: {
  text: string;
  label: string;
  icon?: boolean;
}) {
  const [done, setDone] = useState(false);
  const copy = async () => {
    await navigator.clipboard?.writeText(text);
    setDone(true);
    setTimeout(() => setDone(false), 1500);
  };
  const Icon = done ? Check : Copy;
  return icon ? (
    <Button variant="outline" size="icon-sm" aria-label={label} onClick={copy} disabled={!text}>
      <Icon />
    </Button>
  ) : (
    <Button variant="outline" onClick={copy} disabled={!text}>
      <Icon />
      {done ? "Copied" : label}
    </Button>
  );
}
