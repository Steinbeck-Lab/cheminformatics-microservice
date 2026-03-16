import { useEffect, useRef, useState, useCallback, type RefObject } from "react";
import { usePreventIframeScrollJack } from "./usePreventIframeScrollJack";

/**
 * Shared hook for Ketcher editor iframe communication.
 *
 * Handles initialization, postMessage communication, command execution,
 * and cleanup for all components that embed Ketcher (StructureDrawView,
 * InChIView, RInChIView).
 */

const ALLOWED_COMMANDS = new Set([
  "getSmiles",
  "getMolfile",
  "getRxn",
  "setMolecule",
  "containsReaction",
]);

interface MessageHandler {
  resolve: (data: unknown) => void;
  reject: (err: Error) => void;
  timeoutId: ReturnType<typeof setTimeout>;
}

interface UseKetcherEditorReturn {
  ketcherFrame: RefObject<HTMLIFrameElement | null>;
  isEditorReady: boolean;
  executeCommand: (command: string, args?: unknown[]) => Promise<unknown>;
  resetEditor: () => void;
}

export function useKetcherEditor(): UseKetcherEditorReturn {
  const ketcherFrame = useRef<HTMLIFrameElement | null>(null);
  const [isEditorReady, setIsEditorReady] = useState(false);
  const messageHandlers = useRef<Record<number, MessageHandler>>({});
  const messageId = useRef(0);
  const [retryCount, setRetryCount] = useState(0);

  usePreventIframeScrollJack(ketcherFrame);

  // Listen for messages from the Ketcher iframe
  useEffect(() => {
    const handleMessage = (event: MessageEvent) => {
      if (event.origin !== window.location.origin) return;
      if (!ketcherFrame.current || event.source !== ketcherFrame.current.contentWindow) return;
      if (typeof event.data !== "object" || event.data === null) return;

      const { id, type, status, data, error } = event.data;

      if (id !== undefined && messageHandlers.current[id]) {
        const handler = messageHandlers.current[id];
        clearTimeout(handler.timeoutId);
        delete messageHandlers.current[id];

        if (status === "success") {
          handler.resolve(data);
        } else if (status === "error") {
          handler.reject(new Error(error || "Unknown Ketcher error"));
        }
      }

      if (type === "ketcher-ready") {
        setIsEditorReady(true);
      }
    };

    window.addEventListener("message", handleMessage);
    return () => {
      window.removeEventListener("message", handleMessage);
      // Clean up pending handlers
      for (const id of Object.keys(messageHandlers.current)) {
        clearTimeout(messageHandlers.current[Number(id)].timeoutId);
      }
      messageHandlers.current = {};
    };
  }, []);

  // Send a message to the Ketcher iframe and await response
  const sendMessage = useCallback((type: string, payload: Record<string, unknown> = {}) => {
    return new Promise<unknown>((resolve, reject) => {
      if (!ketcherFrame.current || !ketcherFrame.current.contentWindow) {
        reject(new Error("Ketcher frame not available"));
        return;
      }

      const id = messageId.current++;

      const timeoutId = setTimeout(() => {
        if (messageHandlers.current[id]) {
          delete messageHandlers.current[id];
          reject(new Error("Timeout waiting for Ketcher response"));
        }
      }, 10000);

      messageHandlers.current[id] = {
        resolve: (data: unknown) => {
          clearTimeout(timeoutId);
          resolve(data);
        },
        reject: (err: Error) => {
          clearTimeout(timeoutId);
          reject(err);
        },
        timeoutId,
      };

      const message = { id, type, payload };
      try {
        ketcherFrame.current.contentWindow.postMessage(message, window.location.origin);
      } catch (err) {
        clearTimeout(timeoutId);
        delete messageHandlers.current[id];
        reject(
          new Error(`Failed to send message: ${err instanceof Error ? err.message : String(err)}`)
        );
      }
    });
  }, []);

  // Execute a Ketcher command with allowlist enforcement
  const executeCommand = useCallback(
    async (command: string, args: unknown[] = []): Promise<unknown> => {
      if (!ALLOWED_COMMANDS.has(command)) {
        throw new Error(`Ketcher command not allowed: ${command}`);
      }

      // Try direct access first (faster, no postMessage overhead)
      try {
        if (ketcherFrame.current && ketcherFrame.current.contentWindow) {
          const ketcher = (
            ketcherFrame.current.contentWindow as Window & { ketcher?: Record<string, unknown> }
          ).ketcher;
          if (ketcher && typeof ketcher[command] === "function") {
            return await (ketcher[command] as (...a: unknown[]) => Promise<unknown>)(...args);
          }
        }
      } catch {
        // Direct access failed (cross-origin or not ready), fall through to postMessage
      }

      // Fall back to message passing
      return sendMessage("ketcher-command", { command, args });
    },
    [sendMessage]
  );

  // Inject the communication bridge script into the Ketcher iframe
  const initializeKetcher = useCallback(() => {
    const injectScript = () => {
      try {
        if (!ketcherFrame.current || !ketcherFrame.current.contentWindow) return;

        const iframeWindow = ketcherFrame.current.contentWindow;
        const iframeDocument = iframeWindow.document;

        const script = iframeDocument.createElement("script");
        script.textContent = `
          window.addEventListener('message', async function(event) {
            if (event.source !== window.parent) return;
            const { id, type, payload } = event.data;
            if (type === 'ketcher-command') {
              try {
                const { command, args } = payload;
                if (!window.ketcher) throw new Error('Ketcher not initialized');
                if (typeof window.ketcher[command] !== 'function') {
                  throw new Error('Command ' + command + ' not available');
                }
                const result = await window.ketcher[command](...(args || []));
                event.source.postMessage({ id, status: 'success', data: result }, event.origin);
              } catch (error) {
                event.source.postMessage({ id, status: 'error', error: error.message }, event.origin);
              }
            }
          });
          var checkKetcher = function() {
            if (window.ketcher) {
              window.parent.postMessage({ type: 'ketcher-ready' }, window.location.origin);
            } else {
              setTimeout(checkKetcher, 100);
            }
          };
          checkKetcher();
        `;
        iframeDocument.head.appendChild(script);
      } catch (err) {
        if (err instanceof DOMException && err.name === "SecurityError") {
          console.debug("Ketcher iframe cross-origin, using postMessage mode");
        } else {
          console.error("Failed to inject Ketcher communication script:", err);
        }
      }
    };

    if (ketcherFrame.current) {
      ketcherFrame.current.onload = () => {
        setTimeout(injectScript, 500);
      };
    }
  }, []);

  // Initialize on mount and on retry
  useEffect(() => {
    initializeKetcher();
  }, [retryCount, initializeKetcher]);

  // Reset editor state and re-initialize
  const resetEditor = useCallback(() => {
    setIsEditorReady(false);
    setRetryCount((prev) => prev + 1);
  }, []);

  return {
    ketcherFrame,
    isEditorReady,
    executeCommand,
    resetEditor,
  };
}
