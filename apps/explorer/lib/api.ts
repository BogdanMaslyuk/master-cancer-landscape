const API_BASE = process.env.NEXT_PUBLIC_MCL_API_URL || "http://127.0.0.1:8000";

// MCL Explorer reads immutable/slow-changing processed research outputs.
// Cache ordinary server-side API reads between route transitions. A small number
// of research payloads are larger than Next.js' 2 MB Data Cache item limit; those
// are fetched without the Next data cache while FastAPI still serves them from its
// in-process cache.
const API_REVALIDATE_SECONDS = Number(process.env.MCL_API_REVALIDATE_SECONDS || 300);

function isOversizedPayload(path: string): boolean {
  return path.startsWith("/api/pathways/stability");
}

export async function apiGet<T>(path: string): Promise<T> {
  const response = await fetch(
    `${API_BASE}${path}`,
    isOversizedPayload(path)
      ? { cache: "no-store" }
      : {
          cache: "force-cache",
          next: { revalidate: API_REVALIDATE_SECONDS },
        },
  );
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
