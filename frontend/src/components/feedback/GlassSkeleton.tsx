import { cn } from "@/lib/utils";

const variantClasses: Record<string, string> = {
  text: "h-4 w-full",
  image: "h-48 w-full",
  card: "h-32 w-full",
  row: "h-6 w-3/4",
  circle: "h-12 w-12 rounded-full",
};

interface GlassSkeletonProps {
  className?: string;
  variant?: "text" | "image" | "card" | "row" | "circle";
}

export function GlassSkeleton({ className, variant = "text" }: GlassSkeletonProps) {
  return (
    <div
      data-testid="glass-skeleton"
      className={cn(
        "relative overflow-hidden rounded-lg bg-white/10 dark:bg-slate-700/20 backdrop-blur-sm",
        variantClasses[variant],
        className
      )}
    >
      <div className="absolute inset-0 animate-shimmer bg-gradient-to-r from-transparent via-white/20 dark:via-white/5 to-transparent" />
    </div>
  );
}
