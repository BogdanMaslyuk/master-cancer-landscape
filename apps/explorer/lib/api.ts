const API_BASE = process.env.NEXT_PUBLIC_MCL_API_URL || "http://127.0.0.1:8000";

// MCL Explorer reads immutable/slow-changing processed research outputs.
// Cache server-side API reads between route transitions; after rebuilding MCL data,
// restart the frontend (or wait for revalidation) to refresh the view.
const API_REVALIDATE_SECONDS = Number(process.env.MCL_API_REVALIDATE_SECONDS || 300);

export async function apiGet<T>(path: string): Promise<T> {
  const response = await fetch(`${API_BASE}${path}`, {
    cache: "force-cache",
    next: { revalidate: API_REVALIDATE_SECONDS },
  });
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
