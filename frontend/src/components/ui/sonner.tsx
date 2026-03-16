import { Toaster as SonnerToaster } from "sonner";
import { useAppContext } from "@/context/AppContext";

export function Toaster() {
  const { isDarkMode } = useAppContext();

  return (
    <SonnerToaster
      position="bottom-right"
      theme={isDarkMode ? "dark" : "light"}
      toastOptions={{
        duration: 3000,
        classNames: {
          toast:
            "!backdrop-blur-xl !bg-white/70 dark:!bg-slate-800/70 !border !border-white/20 dark:!border-slate-600/30 !shadow-lg",
          error: "!border-l-4 !border-l-red-500",
          success: "!border-l-4 !border-l-green-500",
          description: "!text-muted-foreground",
        },
      }}
    />
  );
}
