type IconName = "lung" | "pancreas" | "brain" | "blood" | "tissue" | string;

export default function AtlasIcon({ name }: { name: IconName }) {
  if (name === "lung") return <svg viewBox="0 0 48 48" aria-hidden="true"><path d="M22 8v14c-3-7-6-11-9-12-3 5-5 12-5 19 0 7 4 11 9 11 4 0 6-3 6-8V8h-1Z"/><path d="M26 8v14c3-7 6-11 9-12 3 5 5 12 5 19 0 7-4 11-9 11-4 0-6-3-6-8V8h1Z"/></svg>;
  if (name === "pancreas") return <svg viewBox="0 0 48 48" aria-hidden="true"><path d="M9 26c5-8 12-11 20-9 4 1 7 0 10-2 2-1 4 1 3 3-3 8-10 13-19 13-6 0-10 3-14 7-2 2-5 0-4-3 1-4 2-6 4-9Z"/><circle cx="31" cy="22" r="2.5"/></svg>;
  if (name === "brain") return <svg viewBox="0 0 48 48" aria-hidden="true"><path d="M24 8c-5-5-12-1-11 5-5 0-7 7-3 10-5 4-1 12 5 11 1 6 8 8 11 3 3 5 10 3 10-3 6 1 10-7 5-11 4-3 2-10-3-10 1-6-6-10-11-5Z"/><path d="M24 10v27M16 18c4 0 6 2 8 5M33 16c-4 0-7 3-9 6M15 29c4-1 7 1 9 4M34 29c-4-1-7 1-10 4"/></svg>;
  if (name === "blood") return <svg viewBox="0 0 48 48" aria-hidden="true"><path d="M24 7c6 10 13 17 13 25a13 13 0 1 1-26 0c0-8 7-15 13-25Z"/><circle cx="20" cy="30" r="3"/><circle cx="29" cy="26" r="2.5"/><circle cx="28" cy="35" r="2"/></svg>;
  return <svg viewBox="0 0 48 48" aria-hidden="true"><circle cx="24" cy="24" r="15"/><circle cx="24" cy="24" r="5"/><path d="M24 9v10M24 29v10M9 24h10M29 24h10"/></svg>;
}
