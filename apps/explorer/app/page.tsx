import Link from "next/link";
import { apiGet, formatNumber } from "../lib/api";
import ResearchTrail from "../components/ResearchTrail";

type Summary = {
  genes_analyzed_n: number;
  comparisons_n: number;
  stable_recurrent_genes_n: number;
  stable_pathways_n: number;
  thresholds: number[];
  per_threshold: Record<string, { query_genes_n?: number; significant_terms_n?: number }>;
  data_release?: string;
};

type Term = Record<string, any>;

const russianModule: Record<string, { title: string; meaning: string }> = {
  "CORUM:4200": {
    title: "Митохондриальная динамика",
    meaning: "Устойчивый сигнал комплекса DNM1L–FIS1, связанного с делением и перестройкой митохондрий.",
  },
  "CORUM:1181": {
    title: "Сплайсинг и процессинг РНК",
    meaning: "Устойчивый сигнал C-комплекса сплайсосомы; в пересечении повторяются HNRNPA1, SNRPB2 и SYNCRIP.",
  },
};

export default async function OverviewPage() {
  const [s, stability] = await Promise.all([
    apiGet<Summary>("/api/summary"),
    apiGet<Term[]>("/api/pathways/stability"),
  ]);
  const stableTerms = stability.filter((x:any) => String(x.significant_all_thresholds).toLowerCase() === "true");

  return <>
    <section className="hero">
      <div className="eyebrow hero-eyebrow">MCL EXPLORER v0.1</div>
      <h1>От опухолевого контекста к лекарственной гипотезе</h1>
      <p>Интерфейс помогает пройти от полногеномных CRISPR-данных DepMap к устойчивым генетическим зависимостям и функциональным модулям. Он не выбирает «лучший препарат», а показывает, какие гипотезы уже поддерживаются текущими данными и что ещё нужно проверить.</p>
    </section>

    <ResearchTrail current={1} />

    <section className="grid cards">
      <div className="card"><div className="label">Генов в одном полногеномном сравнении</div><div className="value">{formatNumber(s.genes_analyzed_n,0)}</div><div className="section-copy">Почти весь доступный CRISPR-набор DepMap</div></div>
      <div className="card"><div className="label">Основных сравнений</div><div className="value">{s.comparisons_n}</div><div className="section-copy">Независимые молекулярные контексты для поиска повторяющихся зависимостей</div></div>
      <div className="card"><div className="label">Устойчивых повторяющихся генов</div><div className="value">{s.stable_recurrent_genes_n}</div><div className="section-copy">Сохраняются при Top-50, Top-100 и Top-200</div></div>
      <div className="card"><div className="label">Устойчивых функциональных сигналов</div><div className="value">{s.stable_pathways_n}</div><div className="section-copy">Значимы при всех трёх порогах отбора</div></div>
    </section>

    <section className="section">
      <div className="section-header"><div><div className="eyebrow">ТЕКУЩИЙ БИОЛОГИЧЕСКИЙ РЕЗУЛЬТАТ</div><h2>Что уже выделяется из данных</h2><div className="section-copy">Сначала — биологический смысл. Точные p/q-значения, идентификаторы и технические детали доступны глубже по страницам.</div></div></div>
      <div className="pathway-grid">
        {stableTerms.map((term:any) => {
          const ru = russianModule[String(term.term_id)] || { title: term.term_name, meaning: "Функциональный сигнал устойчив к изменению порога отбора кандидатов." };
          return <article className="pathway-card" key={`${term.source}-${term.term_id}`}>
            <div className="eyebrow">УСТОЙЧИВЫЙ СИГНАЛ · 3/3 ПОРОГА</div>
            <h3>{ru.title}</h3>
            <p className="section-copy">{ru.meaning}</p>
            <div className="technical-note">Исходный термин: <b>{term.term_name}</b> · {term.source} · {term.term_id}</div>
          </article>
        })}
      </div>
    </section>

    <section className="section">
      <div className="section-header"><div><div className="eyebrow">ПРОВЕРКА УСТОЙЧИВОСТИ</div><h2>Зависит ли вывод от размера списка кандидатов?</h2><div className="section-copy">Мы повторили анализ при трёх порогах отбора. Чем чаще ген или модуль сохраняется, тем меньше результат зависит от произвольного выбора Top-N.</div></div></div>
      <div className="grid cards">
        {s.thresholds.map(t => {
          const d = s.per_threshold?.[String(t)] || {};
          return <div className="card" key={t}><div className="label">Порог Top-{t}</div><div className="value">{d.query_genes_n ?? "—"}</div><div className="section-copy">повторяющихся генов · {d.significant_terms_n ?? "—"} значимых функциональных терминов</div></div>
        })}
      </div>
    </section>

    <section className="section split">
      <div className="callout">
        <div className="eyebrow">ЧТО ЭТО ЗНАЧИТ</div>
        <h3>Мы нашли не готовые лекарственные мишени, а устойчивые биологические зависимости</h3>
        <p className="section-copy">CRISPR-нокаут показывает, насколько клетка зависит от гена. Следующий этап должен проверить лекарственную достижимость, селективность, известные лиганды, структуры и риски безопасности.</p>
      </div>
      <div className="card next-action">
        <div className="eyebrow">С ЧЕГО НАЧАТЬ</div>
        <h3>Выберите опухолевый контекст</h3>
        <p className="section-copy">Посмотрите, какие зависимости отличают целевую группу клеточных моделей от группы сравнения, а затем переходите к устойчивым генам и модулям.</p>
        <Link href="/comparisons" className="primary-link">Перейти к опухолевым контекстам →</Link>
      </div>
    </section>

    <section className="section technical-note"><b>Релиз данных DepMap:</b> {s.data_release || "—"}. Английские названия генов, баз данных и стандартных метрик сохраняются там, где это необходимо для сопоставления с исходными данными и литературой.</section>
  </>;
}
