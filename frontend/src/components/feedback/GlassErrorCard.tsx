import { AlertTriangle, RotateCcw } from "lucide-react";
import { Button } from "@/components/ui/button";
import { cn } from "@/lib/utils";

interface GlassErrorCardProps {
  message: string;
  onRetry?: () => void;
  className?: string;
}

/** Madder-tinted error panel: icon, message and an optional outline "Try again". */
export function GlassErrorCard({ message, onRetry, className }: GlassErrorCardProps) {
  return (
    <div
      role="alert"
      className={cn(
        "flex items-start gap-3 rounded-2xl border border-[color-mix(in_srgb,var(--madder)_35%,transparent)] bg-[color-mix(in_srgb,var(--madder)_8%,transparent)] p-4",
        className
      )}
    >
      <AlertTriangle className="mt-0.5 size-5 shrink-0 text-[var(--madder)]" aria-hidden="true" />
      <div className="grid flex-1 justify-items-start gap-3">
        <p className="m-0 text-sm text-foreground">{message}</p>
        {onRetry && (
          <Button variant="outline" size="sm" onClick={onRetry}>
            <RotateCcw aria-hidden="true" />
            Try again
          </Button>
        )}
      </div>
    </div>
  );
}
