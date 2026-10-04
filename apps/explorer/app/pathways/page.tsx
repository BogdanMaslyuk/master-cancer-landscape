import Link from "next/link";
import { apiGet, formatNumber } from "../../lib/api";
import ResearchTrail from "../../components/ResearchTrail";

type Term=Record<string,any>;

function genes(value: unknown): string[] {
  if (!value) return [];
  try { return JSON.parse(String(value)); } catch { return []; }
}

const russianMeaning: Record<string, {title:string; text:string}> = {
  "CORUM:4200": { title:"Митохондриальная динамика", text:"Повторяющийся сигнал комплекса DNM1L–FIS1 указывает на устойчивую зависимость от механизмов деления и перестройки митохондрий." },
  "CORUM:1181": { title:"Сплайсинг и процессинг РНК", text:"Повторяющийся сигнал C-комплекса сплайсосомы указывает на устойчивую зависимость от обработки пре-мРНК." },
};

export default async function PathwaysPage(){
  const stability=await apiGet<Term[]>("/api/pathways/stability?min_significant_thresholds=2&limit=500");
  const all=await apiGet<Term[]>("/api/pathways?significant_only=true&limit=500");
  const stable=stability.filter((x:any)=>String(x.significant_all_thresholds).toLowerCase()==="true");
  const recurrent=stability.filter((x:any)=>Number(x.significant_thresholds_n)>=2 && String(x.significant_all_thresholds).toLowerCase()!=="true");

  return <>
    <div className="section-header" style={{marginBottom:18}}>
      <div>
        <div className="eyebrow">ШАГ 4 · ФУНКЦИОНАЛЬНАЯ ИНТЕРПРЕТАЦИЯ</div>
        <h1 style={{margin:"4px 0 8px"}}>Функциональные модули</h1>
        <div className="section-copy">Здесь мы переходим от отдельных генов к биологическим процессам и белковым комплексам. Источники GO, Reactome, KEGG и CORUM сохраняются как техническое доказательство, но основной акцент сделан на понятном биологическом смысле.</div>
      </div>
    </div>

    <ResearchTrail current={4} />

    <section className="kpi-strip" style={{marginTop:18}}>
      <div className="kpi"><strong>{stable.length}</strong><span>устойчивых сигналов при 50/100/200</span></div>
      <div className="kpi"><strong>{recurrent.length}</strong><span>сигналов, значимых минимум при двух порогах</span></div>
      <div className="kpi"><strong>{all.filter((x:any)=>Number(x.top_n)===50).length}</strong><span>значимых терминов при Top-50</span></div>
      <div className="kpi"><strong>{all.filter((x:any)=>Number(x.top_n)===200).length}</strong><span>значимых терминов при Top-200</span></div>
    </section>

    <section className="section">
      <div className="section-header"><div><div className="eyebrow">НАИБОЛЕЕ УСТОЙЧИВЫЕ СИГНАЛЫ</div><h2>Что сохраняется при всех трёх порогах</h2><div className="section-copy">Каждый функциональный сигнал показан один раз. Три ячейки ниже отражают его значимость при Top-50, Top-100 и Top-200.</div></div></div>
      <div className="pathway-grid">
        {stable.map((x:any)=>{
          const union=[...new Set([...genes(x.intersecting_gene_symbols_top50),...genes(x.intersecting_gene_symbols_top100),...genes(x.intersecting_gene_symbols_top200)])];
          const ru=russianMeaning[String(x.term_id)] || {title:x.term_name,text:"Функциональный сигнал устойчив к изменению порога отбора кандидатов."};
          return <article className="pathway-card" key={`${x.source}-${x.term_id}`}>
            <div className="eyebrow">{ru.title}</div>
            <h3>{x.term_name}</h3>
            <p className="section-copy">{ru.text}</p>
            <div className="gene-cloud">{union.map((g:string)=><Link className="gene-pill" href={`/genes/${g}`} key={g}>{g}</Link>)}</div>
            <div className="pathway-thresholds">
              {[50,100,200].map(t=>{
                const on=String(x[`significant_top${t}`]).toLowerCase()==="true";
                return <div key={t} className={`threshold-cell ${on?"on":"off"}`}><b>Top-{t}</b><br/>{on?`скорр. p = ${formatNumber(x[`p_value_adjusted_top${t}`],4)}`:"не значим"}</div>
              })}
            </div>
            <div className="technical-note">Источник: <b>{x.source}</b> · {x.term_id}</div>
          </article>
        })}
      </div>
    </section>

    <section className="section">
      <div className="section-header"><div><h2>Воспроизводимые, но менее устойчивые сигналы</h2><div className="section-copy">Эти процессы значимы минимум при двух порогах, но не проходят критерий 3/3. Их следует рассматривать как гипотезы второго уровня.</div></div></div>
      <div className="table-wrap"><table><thead><tr><th>Источник</th><th>Функциональный термин</th><th>При каких порогах значим</th><th>Минимальное скорректированное p</th></tr></thead><tbody>{recurrent.slice(0,40).map((x:any)=><tr key={`${x.source}-${x.term_id}`}><td>{x.source}</td><td>{x.term_name}</td><td>{x.thresholds_significant}</td><td>{formatNumber(x.min_p_value_adjusted,6)}</td></tr>)}</tbody></table></div>
    </section>

    <section className="section callout">
      <div className="eyebrow">ВАЖНО ДЛЯ РАЗРАБОТКИ ПРЕПАРАТА</div>
      <h3>Функциональный модуль — это направление для проверки, а не готовый механизм действия</h3>
      <p className="section-copy">Статистическое обогащение показывает, что определённые процессы встречаются среди кандидатов чаще ожидаемого. Дальше нужно проверять конкретные белки внутри модуля: лекарственную достижимость, селективность, известные лиганды, структуры и экспериментальную валидируемость.</p>
    </section>

    <section className="section next-step-banner">
      <div><div className="eyebrow">СЛЕДУЮЩИЙ ЭТАП M3.4</div><h3>От биологической зависимости к фармакологической мишени</h3><p className="section-copy">Следующий слой MCL должен добавить для каждого кандидата лекарственную достижимость, известные молекулы, структуры, данные Open Targets и ограничения безопасности.</p></div>
      <Link href="/methodology" className="primary-link">Посмотреть, что будет оцениваться дальше →</Link>
    </section>

    <section className="section">
      <details className="details-block"><summary>Показать все значимые результаты обогащения</summary><div className="table-wrap" style={{marginTop:12}}><table><thead><tr><th>Top-N</th><th>Источник</th><th>Термин</th><th>Скорректированное p</th><th>Гены пересечения</th></tr></thead><tbody>{all.map((x:any,i:number)=><tr key={`${x.term_id}-${x.top_n}-${i}`}><td>{x.top_n}</td><td>{x.source}</td><td>{x.term_name}</td><td>{formatNumber(x.p_value_adjusted,6)}</td><td>{x.intersecting_gene_symbols_json}</td></tr>)}</tbody></table></div></details>
    </section>
  </>;
}
