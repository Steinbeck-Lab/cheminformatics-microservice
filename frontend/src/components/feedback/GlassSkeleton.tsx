import { cn } from "@/lib/utils";

const variantClasses: Record<string, string> = {
  text: "h-4 w-full",
  image: "h-48 w-full",
  row: "h-6 w-3/4",
};

interface GlassSkeletonProps {
  className?: string;
  variant?: "text" | "image" | "row";
}

/** A calm placeholder block; a silk sheen glides across it only when the viewer allows motion. */
export function GlassSkeleton({ className, variant = "text" }: GlassSkeletonProps) {
  return (
    <div
      data-testid="glass-skeleton"
      aria-hidden="true"
      className={cn(
        "pc-sheen rounded-lg border border-[var(--line)] bg-[var(--inset)]",
        variantClasses[variant],
        className
      )}
    />
  );
}
