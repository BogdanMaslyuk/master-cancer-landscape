import "./globals.css";
import "./research.css";
import AppShell from "../components/AppShell";

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
