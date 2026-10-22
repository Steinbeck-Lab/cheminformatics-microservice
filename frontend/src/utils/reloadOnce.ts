const FLAG = "chunk-reload";

/**
 * A page's code failed to load, usually because a new deploy replaced the hashed file an open tab
 * still asks for. Reload once to pick up the new build; if it fails again, let the error show.
 */
export function reloadOnce(
  err: unknown,
  store: Pick<Storage, "getItem" | "setItem"> = sessionStorage,
  reload: () => void = () => window.location.reload()
): Promise<never> {
  try {
    if (store.getItem(FLAG)) return Promise.reject(err);
    store.setItem(FLAG, "1");
  } catch {
    return Promise.reject(err); // storage blocked: never risk a reload loop
  }
  reload();
  return new Promise<never>(() => {}); // the reload replaces the page
}

/** A page loaded fine: allow one reload again for the next failure. */
export function clearReloadFlag(store: Pick<Storage, "removeItem"> = sessionStorage) {
  try {
    store.removeItem(FLAG);
  } catch {
    // storage blocked: nothing to clear
  }
}
