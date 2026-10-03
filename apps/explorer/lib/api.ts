const API_BASE = process.env.NEXT_PUBLIC_MCL_API_URL || "http://127.0.0.1:8000";

export async function apiGet<T>(path: string): Promise<T> {
  const response = await fetch(`${API_BASE}${path}`, { cache: "no-store" });
  if (!response.ok) {
    const text = await response.text();
    throw new Error(`${response.status}: ${text}`);
  }
  return response.json() as Promise<T>;
}

export function formatNumber(value: unknown, digits = 3): string {
  if (value === null || value === undefined || value === "") return "—";
  const n = Number(value);
  if (!Number.isFinite(n)) return String(value);
  if (Math.abs(n) >= 1000) return new Intl.NumberFormat("ru-RU").format(n);
  return n.toFixed(digits).replace(/\.0+$/, "");
}
