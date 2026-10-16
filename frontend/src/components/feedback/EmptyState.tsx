import { Beaker } from "lucide-react";
import { cn } from "@/lib/utils";

interface EmptyStateProps {
  message?: string;
  icon?: React.ReactNode;
  className?: string;
}

/** Dashed empty panel, the component form of .pc-empty. */
export function EmptyState({
  message = "Enter input above and submit to see results",
  icon,
  className,
}: EmptyStateProps) {
  return (
    <div
      data-testid="empty-state"
      className={cn("pc-empty grid justify-items-center gap-2", className)}
    >
      {icon ?? <Beaker className="size-6 text-muted-foreground" aria-hidden="true" />}
      <p className="m-0">{message}</p>
    </div>
  );
}
