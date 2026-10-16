/** A molecular formula with element colours and real subscripts (C8H10N4O2 → C₈H₁₀N₄O₂). */
export function Formula({ value, as: Tag = "dd" }: { value: string; as?: "dd" | "span" }) {
  return (
    <Tag className="pc-formula pc-val" data-testid="formula">
      {[...value.matchAll(/([A-Z][a-z]?)(\d*)|([^A-Za-z\d]+)/g)].map((m, i) =>
        m[1] ? (
          <span key={i} className={`pc-el-${m[1]}`}>
            {m[1]}
            {m[2] && <sub>{m[2]}</sub>}
          </span>
        ) : (
          <span key={i}>{m[3]}</span>
        )
      )}
    </Tag>
  );
}
