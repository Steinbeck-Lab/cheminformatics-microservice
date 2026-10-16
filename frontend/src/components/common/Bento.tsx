import { useId, type ComponentProps, type ReactNode } from "react";
import { cn } from "@/lib/utils";

/** Glass title cell at the top of a tool page. */
export function ToolPageHead({
  title,
  description,
  status,
}: {
  title: string;
  description?: ReactNode;
  status?: ReactNode;
}) {
  return (
    <header className="glass my-3 flex flex-wrap items-end justify-between gap-4 rounded-cell px-4 py-[18px] sm:px-6">
      <div className="min-w-0">
        <h1 className="pc-display m-0 text-[30px] leading-[1.1] text-head sm:text-[36px]">
          {title}
        </h1>
        {description && (
          <p className="mt-1.5 mb-0 max-w-[96ch] text-[15px] text-muted-foreground">
            {description}
          </p>
        )}
      </div>
      {status && <div className="text-[13px] text-muted-foreground">{status}</div>}
    </header>
  );
}

/** 12-column bento grid. Cells are full width until a page places them with `lg:` classes. */
export function Bento({ className, ...props }: ComponentProps<"div">) {
  return <div className={cn("grid grid-cols-12 gap-4 pb-4", className)} {...props} />;
}

interface BentoCellProps extends Omit<ComponentProps<"section">, "title"> {
  /** Uppercase cell label; also the region's accessible name. */
  title: string;
  /** Small text or controls at the right end of the label row. */
  meta?: ReactNode;
  /** Madder dot at the end of the label rule (input and NP-likeness cells). */
  accent?: "madder";
}

export function BentoCell({ title, meta, accent, className, children, ...props }: BentoCellProps) {
  const id = useId();
  return (
    <section
      aria-labelledby={id}
      className={cn("glass col-span-full min-w-0 rounded-cell p-4 sm:p-5", className)}
      {...props}
    >
      <div className="mb-4 flex min-h-8 items-center gap-3.5">
        <h2 id={id} className="bento-label" data-accent={accent}>
          {title}
        </h2>
        {meta && (
          <div className="ml-auto flex items-center gap-2 text-[12.5px] text-muted-foreground">
            {meta}
          </div>
        )}
      </div>
      {children}
    </section>
  );
}
