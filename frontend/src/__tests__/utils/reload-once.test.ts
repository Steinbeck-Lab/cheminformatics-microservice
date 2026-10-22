import { describe, it, expect, vi } from "vitest";
import { reloadOnce } from "@/utils/reloadOnce";

const memory = () => {
  const m = new Map<string, string>();
  return {
    getItem: (k: string) => m.get(k) ?? null,
    setItem: (k: string, v: string) => void m.set(k, v),
  };
};

describe("reloadOnce", () => {
  it("reloads on the first failed chunk, then lets the next failure through", async () => {
    const store = memory();
    const reload = vi.fn();
    void reloadOnce(new Error("chunk 404"), store, reload);
    expect(reload).toHaveBeenCalledTimes(1);
    await expect(reloadOnce(new Error("chunk 404"), store, reload)).rejects.toThrow("chunk 404");
    expect(reload).toHaveBeenCalledTimes(1);
  });
});
