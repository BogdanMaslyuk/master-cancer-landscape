import "./globals.css";
import "./research.css";
import "./atlas.css";
import "./atlas-polish.css";
import AppShell from "../components/AppShell";

// MCL Explorer is a runtime-backed research UI over the local FastAPI service.
// Do not prerender API-dependent pages during `next build`; the scientific data
// are materialized separately and served by the backend at request time.
export const dynamic = "force-dynamic";

export const metadata = {
  title: "MCL Explorer",
  description: "Master Cancer Landscape — исследовательский интерфейс",
};

export default function RootLayout({ children }: { children: React.ReactNode }) {
  return (
    <html lang="ru">
      <body>
        <AppShell>{children}</AppShell>
      </body>
    </html>
  );
}
