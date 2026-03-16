import { Beaker } from "lucide-react";
import { cn } from "@/lib/utils";

interface EmptyStateProps {
  message?: string;
  icon?: React.ReactNode;
  className?: string;
}

export function EmptyState({
  message = "Enter input above and submit to see results",
  icon,
  className,
}: EmptyStateProps) {
  return (
    <div
      data-testid="empty-state"
      className={cn("glass-bold rounded-xl p-8 text-center", className)}
    >
      <div className="flex flex-col items-center gap-3">
        {icon ?? <Beaker className="h-12 w-12 text-muted-foreground" />}
        <p className="text-sm text-muted-foreground">{message}</p>
      </div>
    </div>
  );
}
