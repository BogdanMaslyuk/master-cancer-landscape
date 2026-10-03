import Link from "next/link";
import { Fragment } from "react";
import { apiGet } from "../../lib/api";
import ResearchTrail from "../../components/ResearchTrail";

type Gene = Record<string, any>;

export default async function GenesPage(){
  const items = await apiGet<Gene[]>("/api/genes/stable");
  return <>
    <div className="section-header" style={{marginBottom:18}}>
      <div>
        <div className="eyebrow">ШАГ 3 · ПРОВЕРКА УСТОЙЧИВОСТИ</div>
        <h1 style={{margin:"4px 0 8px"}}>Устойчивые гены-кандидаты</h1>
        <div className="section-copy">Здесь собраны гены, которые сохраняют статус повторяющейся зависимости при Top-50, Top-100 и Top-200. Это ещё не рейтинг лекарственных мишеней, а наиболее устойчивое ядро текущего CRISPR-анализа.</div>
      </div>
    </div>

    <ResearchTrail current={3} />

    <section className="split section">
      <div className="card">
        <div className="eyebrow">УСТОЙЧИВОЕ ЯДРО</div>
        <div className="value">{items.length}</div>
        <div className="section-copy">Генов, повторяющихся при всех трёх порогах отбора.</div>
      </div>
      <div className="callout">
        <div className="eyebrow">ПОЧЕМУ ЭТО ВАЖНО</div>
        <h3>Мы уменьшаем зависимость вывода от произвольного Top-N</h3>
        <p className="section-copy">Если кандидат остаётся повторяющимся при Top-50, Top-100 и Top-200, его присутствие меньше зависит от того, насколько широкий список мы решили анализировать.</p>
      </div>
    </section>

    <section className="section">
      <div className="section-header"><div><h2>Матрица устойчивости</h2><div className="section-copy">Точка означает, что ген сохраняет повторяющийся статус при соответствующем пороге отбора.</div></div></div>
      <div className="matrix" style={{marginTop:12}}>
        <div className="head">Ген</div><div className="head">Top-50</div><div className="head">Top-100</div><div className="head">Top-200</div>
        {items.map((g:any)=><Fragment key={g.gene_symbol}>
          <div className="name"><Link href={`/genes/${g.gene_symbol}`} style={{color:"var(--accent)"}}>{g.gene_symbol}</Link></div>
          <div className={String(g.recurrent_top50).toLowerCase()==="true"?"yes":"no"}>{String(g.recurrent_top50).toLowerCase()==="true"?"●":"—"}</div>
          <div className={String(g.recurrent_top100).toLowerCase()==="true"?"yes":"no"}>{String(g.recurrent_top100).toLowerCase()==="true"?"●":"—"}</div>
          <div className={String(g.recurrent_top200).toLowerCase()==="true"?"yes":"no"}>{String(g.recurrent_top200).toLowerCase()==="true"?"●":"—"}</div>
        </Fragment>)}
      </div>
    </section>

    <section className="section">
      <div className="section-header"><div><h2>Все устойчивые гены</h2><div className="section-copy">Клик по гену открывает его результаты по каждому сравнению, связанные функциональные сигналы и ограничения интерпретации.</div></div></div>
      <div className="gene-cloud" style={{marginTop:12}}>{items.map((g:any)=><Link className="gene-pill" key={g.gene_symbol} href={`/genes/${g.gene_symbol}`}>{g.gene_symbol}</Link>)}</div>
    </section>

    <section className="section next-step-banner">
      <div><div className="eyebrow">СЛЕДУЮЩИЙ ШАГ</div><h3>Понять, какие биологические системы объединяют эти гены</h3><p className="section-copy">Отдельный ген труднее интерпретировать, чем повторяющийся функциональный модуль. Поэтому дальше переходим к обогащению процессов и белковых комплексов.</p></div>
      <Link href="/pathways" className="primary-link">Перейти к функциональным модулям →</Link>
    </section>
  </>;
}
