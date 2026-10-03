import Link from "next/link";
import { apiGet } from "../../lib/api";

type Gene = Record<string, any>;

export default async function GenesPage(){
  const items = await apiGet<Gene[]>("/api/genes/stable");
  return <>
    <div className="section-header" style={{marginBottom:18}}>
      <div>
        <div className="eyebrow">ROBUSTNESS LAYER</div>
        <h1 style={{margin:"4px 0 8px"}}>Stable Gene Explorer</h1>
        <div className="section-copy">Гены, которые сохраняют recurrent-статус при Top-50, Top-100 и Top-200. Это не рейтинг лекарственных мишеней, а устойчивое dependency-ядро.</div>
      </div>
    </div>

    <section className="split">
      <div className="card">
        <div className="eyebrow">Stable recurrent set</div>
        <div className="value">{items.length}</div>
        <div className="section-copy">Кандидаты, пережившие изменение cut-off. Клик по гену открывает его comparison matrix и pathway associations.</div>
      </div>
      <div className="callout">
        <div className="eyebrow">Interpretation</div>
        <h3>Почему это важнее одного Top-100</h3>
        <p className="section-copy">Если ген остаётся recurrent при трёх разных cut-off, его присутствие меньше зависит от произвольного выбора размера списка кандидатов.</p>
      </div>
    </section>

    <section className="section">
      <h2>Threshold stability matrix</h2>
      <div className="matrix" style={{marginTop:12}}>
        <div className="head">Gene</div><div className="head">Top-50</div><div className="head">Top-100</div><div className="head">Top-200</div>
        {items.map((g:any)=><>
          <div className="name" key={`${g.gene_symbol}-name`}><Link href={`/genes/${g.gene_symbol}`} style={{color:"var(--accent)"}}>{g.gene_symbol}</Link></div>
          <div className={String(g.recurrent_top50).toLowerCase()==="true"?"yes":"no"}>{String(g.recurrent_top50).toLowerCase()==="true"?"●":"—"}</div>
          <div className={String(g.recurrent_top100).toLowerCase()==="true"?"yes":"no"}>{String(g.recurrent_top100).toLowerCase()==="true"?"●":"—"}</div>
          <div className={String(g.recurrent_top200).toLowerCase()==="true"?"yes":"no"}>{String(g.recurrent_top200).toLowerCase()==="true"?"●":"—"}</div>
        </>)}
      </div>
    </section>

    <section className="section">
      <h2>Stable gene set</h2>
      <div className="gene-cloud" style={{marginTop:12}}>{items.map((g:any)=><Link className="gene-pill" key={g.gene_symbol} href={`/genes/${g.gene_symbol}`}>{g.gene_symbol}</Link>)}</div>
    </section>
  </>;
}
