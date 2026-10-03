import Link from "next/link";

const steps = [
  { n: 1, href: "/atlas", title: "Выбрать орган, опухоль и молекулярный контекст", short: "Заболевание" },
  { n: 2, href: "/atlas", title: "Проверить клеточные модели и выбрать сравнение", short: "Модели" },
  { n: 3, href: "/genes", title: "Оценить устойчивые гены", short: "Гены" },
  { n: 4, href: "/pathways", title: "Понять функциональные модули", short: "Модули" },
  { n: 5, href: "/methodology", title: "Фармакологическая оценка — следующий этап M3.4", short: "M3.4" },
];

export default function ResearchTrail({ current }: { current: number }) {
  return (
    <nav className="research-trail" aria-label="Путь исследователя">
      <div className="trail-title">Путь исследователя</div>
      <div className="trail-steps">
        {steps.map((step) => {
          const state = step.n < current ? "done" : step.n === current ? "current" : "future";
          return (
            <Link key={step.n} href={step.href} className={`trail-step ${state}`} title={step.title}>
              <span className="trail-number">{step.n}</span>
              <span><b>{step.short}</b><small>{step.title}</small></span>
            </Link>
          );
        })}
      </div>
    </nav>
  );
}
