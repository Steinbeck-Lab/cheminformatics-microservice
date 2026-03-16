import { AlertTriangle, RotateCcw } from "lucide-react";
import { motion } from "motion/react";
import { Button } from "@/components/ui/button";
import { cn } from "@/lib/utils";

interface GlassErrorCardProps {
  message: string;
  onRetry?: () => void;
  className?: string;
}

export function GlassErrorCard({ message, onRetry, className }: GlassErrorCardProps) {
  return (
    <motion.div
      initial={{ opacity: 0, y: 8 }}
      animate={{ opacity: 1, y: 0 }}
      transition={{ duration: 0.3 }}
      role="alert"
      className={cn(
        "glass-bold rounded-xl p-6 border-l-4 border-l-red-500 dark:border-l-red-400",
        className
      )}
    >
      <div className="flex items-start gap-3">
        <AlertTriangle className="h-5 w-5 shrink-0 text-red-500 dark:text-red-400 mt-0.5" />
        <div className="flex-1 space-y-3">
          <p className="text-sm text-foreground">{message}</p>
          {onRetry && (
            <Button variant="ghost" size="sm" onClick={onRetry}>
              <RotateCcw className="h-3.5 w-3.5" />
              Try again
            </Button>
          )}
        </div>
      </div>
    </motion.div>
  );
}
