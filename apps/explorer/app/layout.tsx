import "./globals.css";
import Link from "next/link";

export const metadata = {
  title: "MCL Explorer",
  description: "Master Cancer Landscape Explorer",
};

const links = [
  ["/", "Обзор"],
  ["/comparisons", "Опухолевые контексты"],
  ["/genes", "Гены-кандидаты"],
  ["/pathways", "Функциональные модули"],
  ["/network", "Карта связей"],
  ["/qc", "Контроль качества"],
  ["/methodology", "Методика и термины"],
];

export default function RootLayout({ children }: { children: React.ReactNode }) {
  return (
    <html lang="ru">
      <body>
        <div className="shell">
          <header className="topbar">
            <div className="topbar-inner">
              <Link href="/" className="brand">MASTER CANCER LANDSCAPE</Link>
              <nav className="nav">
                {links.map(([href, label]) => <Link key={href} href={href}>{label}</Link>)}
              </nav>
            </div>
          </header>
          <main>{children}</main>
        </div>
      </body>
    </html>
  );
}
