"use client";

import Link from "next/link";
import { usePathname } from "next/navigation";
import GeneSearch from "./GeneSearch";

const researchLinks = [
  { href: "/", label: "Обзор", index: "01", note: "Сводка проекта" },
  { href: "/atlas", label: "Атлас опухолей", index: "02", note: "Орган → опухоль → модели" },
  { href: "/models", label: "Исследователь моделей", index: "03", note: "Поиск · фильтры · карточки" },
  { href: "/genes", label: "Гены и мишени", index: "04", note: "Любой ген → функции → контексты" },
  { href: "/genes/matrix", label: "Ген × опухоль", index: "05", note: "Сравнение зависимостей" },
  { href: "/pathways", label: "Функциональные модули", index: "06", note: "Пути и комплексы" },
  { href: "/network", label: "Карта связей", index: "07", note: "Гены ↔ модули" },
];

const systemLinks = [
  { href: "/qc", label: "Контроль качества", index: "08", note: "Ошибки и ограничения" },
  { href: "/methodology", label: "Методика и термины", index: "09", note: "Как читать MCL" },
];

function isActive(pathname: string, href: string) {
  if (href === "/") return pathname === "/";
  if (href === "/atlas" && pathname.startsWith("/comparisons")) return true;
  if (href === "/genes" && pathname.startsWith("/genes/matrix")) return false;
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
  const inGenes = pathname.startsWith("/genes");
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
            <strong>MCL Explorer v0.6</strong>
            <span>Локальная исследовательская среда</span>
          </div>
        </div>
      </aside>

      <div className="workspace">
        <header className="workspace-bar">
          <div className="workspace-context">
            <span className="workspace-kicker">MASTER CANCER LANDSCAPE</span>
            <span className="workspace-divider" />
            <span>{inGenes ? "От функции и гена к опухолевому контексту" : inAtlas ? "От опухоли к модели и функциональным зависимостям" : "Поиск противоопухолевых зависимостей"}</span>
          </div>
          <GeneSearch compact />
          <div className="workspace-badge">M3.3.1 · исследовательский режим</div>
        </header>
        <main className="workspace-main">{children}</main>
      </div>
    </div>
  );
}
