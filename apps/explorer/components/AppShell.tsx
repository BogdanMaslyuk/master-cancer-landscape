"use client";

import Link from "next/link";
import { usePathname } from "next/navigation";

const researchLinks = [
  { href: "/", label: "Обзор", index: "01", note: "Сводка проекта" },
  { href: "/atlas", label: "Атлас опухолей", index: "02", note: "Орган → опухоль → контекст" },
  { href: "/models", label: "Клеточные линии", index: "03", note: "Модели DepMap" },
  { href: "/genes", label: "Гены-кандидаты", index: "04", note: "Устойчивые зависимости" },
  { href: "/pathways", label: "Функциональные модули", index: "05", note: "Пути и комплексы" },
  { href: "/network", label: "Карта связей", index: "06", note: "Гены ↔ модули" },
];

const systemLinks = [
  { href: "/qc", label: "Контроль качества", index: "07", note: "Ошибки и ограничения" },
  { href: "/methodology", label: "Методика и термины", index: "08", note: "Как читать MCL" },
];

function isActive(pathname: string, href: string) {
  if (href === "/") return pathname === "/";
  return pathname === href || pathname.startsWith(`${href}/`);
}

function NavGroup({ title, items, pathname }: { title: string; items: typeof researchLinks; pathname: string }) {
  return (
    <div className="side-group">
      <div className="side-group-title">{title}</div>
      <nav className="side-nav">
        {items.map((item) => {
          const active = isActive(pathname, item.href);
          return (
            <Link key={item.href} href={item.href} className={`side-link ${active ? "active" : ""}`}>
              <span className="side-link-index">{item.index}</span>
              <span className="side-link-copy">
                <strong>{item.label}</strong>
                <small>{item.note}</small>
              </span>
            </Link>
          );
        })}
      </nav>
    </div>
  );
}

export default function AppShell({ children }: { children: React.ReactNode }) {
  const pathname = usePathname();
  const inAtlas = pathname.startsWith("/atlas") || pathname.startsWith("/models") || pathname.startsWith("/comparisons");
  return (
    <div className="app-shell">
      <aside className="sidebar">
        <Link href="/" className="brand-block">
          <span className="brand-mark">M</span>
          <span>
            <strong>Master Cancer Landscape</strong>
            <small>Исследовательская система</small>
          </span>
        </Link>

        <div className="sidebar-scroll">
          <NavGroup title="ИССЛЕДОВАНИЕ" items={researchLinks} pathname={pathname} />
          <NavGroup title="СИСТЕМА" items={systemLinks} pathname={pathname} />
        </div>

        <div className="sidebar-footer">
          <div className="status-dot" />
          <div>
            <strong>MCL Explorer v0.2</strong>
            <span>Локальная исследовательская среда</span>
          </div>
        </div>
      </aside>

      <div className="workspace">
        <header className="workspace-bar">
          <div className="workspace-context">
            <span className="workspace-kicker">MASTER CANCER LANDSCAPE</span>
            <span className="workspace-divider" />
            <span>{inAtlas ? "От заболевания к клеточной модели" : "Поиск противоопухолевых зависимостей"}</span>
          </div>
          <div className="workspace-badge">M3.3.1 · исследовательский режим</div>
        </header>
        <main className="workspace-main">{children}</main>
      </div>
    </div>
  );
}
