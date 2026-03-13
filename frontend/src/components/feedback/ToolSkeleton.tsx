import { GlassSkeleton } from "./GlassSkeleton";
import { cn } from "@/lib/utils";

interface ToolSkeletonProps {
  variant?: "descriptors" | "molecule" | "conversion" | "general";
  className?: string;
}

function DescriptorsSkeleton() {
  return (
    <div className="space-y-4">
      <GlassSkeleton variant="image" />
      <div className="space-y-2">
        <GlassSkeleton variant="row" />
        <GlassSkeleton variant="text" />
        <GlassSkeleton variant="row" />
        <GlassSkeleton variant="text" />
        <GlassSkeleton variant="row" />
        <GlassSkeleton variant="text" />
      </div>
    </div>
  );
}

function MoleculeSkeleton() {
  return (
    <div className="space-y-4">
      <GlassSkeleton variant="image" className="h-64" />
      <GlassSkeleton variant="text" />
      <GlassSkeleton variant="row" />
      <div className="space-y-2">
        <GlassSkeleton variant="text" />
        <GlassSkeleton variant="row" />
        <GlassSkeleton variant="text" />
        <GlassSkeleton variant="row" />
      </div>
    </div>
  );
}

function ConversionSkeleton() {
  return (
    <div className="space-y-3">
      <GlassSkeleton variant="text" />
      <GlassSkeleton variant="text" />
      <GlassSkeleton variant="text" />
      <GlassSkeleton variant="text" className="w-1/2" />
    </div>
  );
}

function GeneralSkeleton() {
  return (
    <div className="space-y-4">
      <GlassSkeleton variant="text" />
      <GlassSkeleton variant="row" />
      <GlassSkeleton variant="image" />
      <div className="space-y-2">
        <GlassSkeleton variant="text" />
        <GlassSkeleton variant="row" />
        <GlassSkeleton variant="text" />
      </div>
    </div>
  );
}

const skeletonVariants: Record<string, () => React.JSX.Element> = {
  descriptors: DescriptorsSkeleton,
  molecule: MoleculeSkeleton,
  conversion: ConversionSkeleton,
  general: GeneralSkeleton,
};

export function ToolSkeleton({ variant = "general", className }: ToolSkeletonProps) {
  const Skeleton = skeletonVariants[variant] ?? GeneralSkeleton;
  return (
    <div data-testid="tool-skeleton" className={cn("glass-bold rounded-xl p-6", className)}>
      <Skeleton />
    </div>
  );
}
