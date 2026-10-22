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
          toast: "glass !rounded-2xl !text-foreground",
          error: "!border-destructive/50 [&_[data-icon]]:!text-madder-ink",
          success: "[&_[data-icon]]:!text-primary",
          description: "!text-muted-foreground",
        },
      }}
    />
  );
}
