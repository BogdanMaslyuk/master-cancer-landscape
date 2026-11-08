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
  const [s, stableTerms, pyz] = await Promise.all([
    apiGet<Summary>("/api/summary"),
    apiGet<Term[]>("/api/pathways/stability?stable_only=true&limit=50"),
    apiGet<Record<string,any>>("/api/pyz/summary"),
  ]);

  return <>
    <section className="hero">
      <div className="eyebrow hero-eyebrow">MCL EXPLORER v1.0</div>
      <h1>От опухоли и собственной молекулы к проверяемому эксперименту</h1>
      <p>MCL связывает орган, тип опухоли, клеточную модель и CRISPR-зависимость с фармакологией и собственной серией PYZ. Теперь исследование можно начинать как от болезни, так и от конкретной молекулы: PYZ → предполагаемая мишень → опухолевый контекст → клеточные культуры → лабораторная проверка.</p>
    </section>

    <ResearchTrail current={1} />

    <section className="grid cards">
      <div className="card"><div className="label">Собственных молекул PYZ</div><div className="value">{formatNumber(pyz.compounds_n,0)}</div><div className="section-copy">Единый реестр структур, ADMET-снимка и target-гипотез</div></div>
      <div className="card"><div className="label">Core-мишеней PYZ-контура</div><div className="value">{formatNumber(pyz.targets_n,0)}</div><div className="section-copy">Текущие мишени Candidate v2 со статусом priority_for_in_vitro</div></div>
      <div className="card"><div className="label">Пар PYZ × мишень ≥0,35</div><div className="value">{formatNumber(pyz.shortlist_rows_n,0)}</div><div className="section-copy">Исследовательский shortlist, а не подтверждённые механизмы</div></div>
      <div className="card"><div className="label">Устойчивых функциональных сигналов</div><div className="value">{s.stable_pathways_n}</div><div className="section-copy">Значимы при всех трёх порогах отбора</div></div>
    </section>

    <section className="section split">
      <div className="callout">
        <div className="eyebrow">НОВЫЙ МАРШРУТ · ОТ НАШЕЙ МОЛЕКУЛЫ</div>
        <h3>PYZ → мишень → орган → опухоль → клеточная культура</h3>
        <p className="section-copy">Для каждой PYZ теперь можно увидеть химическое сходство с экспериментальными лигандами 13 мишеней, строгие онкологические эталоны, PRISM/CRISPR-контекст известных лигандов, релевантные опухоли и конкретные линии для проверки. Отдельно отмечается, какие культуры физически доступны в нашей лаборатории.</p>
        <Link href="/pyz" className="primary-link">Открыть молекулы PYZ →</Link>
      </div>
      <div className="card next-action">
        <div className="eyebrow">МАТРИЦА ДОКАЗАТЕЛЬСТВ</div>
        <h3>65 молекул × 13 мишеней</h3>
        <p className="section-copy">Матрица показывает лучший Tanimoto до экспериментального лиганда каждой мишени и быстро отделяет 10 текущих исследовательских пар от 835 слабых по 2D-сходству комбинаций.</p>
        <Link href="/pyz/matrix" className="primary-link">Открыть матрицу PYZ × мишень →</Link>
      </div>
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
        <p className="section-copy">CRISPR-нокаут показывает, насколько клетка зависит от гена. Химическое сходство, PRISM, ADMET и докинг добавляются как отдельные оси и не превращаются в непрозрачный итоговый балл.</p>
      </div>
      <div className="card next-action">
        <div className="eyebrow">ВТОРОЙ МАРШРУТ · ОТ БОЛЕЗНИ</div>
        <h3>Откройте атлас опухолей</h3>
        <p className="section-copy">Выберите орган, опухоль и молекулярный контекст, проверьте клеточные линии и зависимости, а затем переходите к веществам и собственным PYZ.</p>
        <Link href="/atlas" className="primary-link">Открыть атлас опухолей →</Link>
      </div>
    </section>

    <section className="section technical-note"><b>Релиз данных DepMap:</b> {s.data_release || "—"}. Английские названия генов, баз данных и стандартных метрик сохраняются там, где это необходимо для сопоставления с исходными данными и литературой.</section>
  </>;
}
