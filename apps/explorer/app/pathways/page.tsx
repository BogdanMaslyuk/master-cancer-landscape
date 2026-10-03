import { apiGet, formatNumber } from "../../lib/api";

type Term=Record<string,any>;

function genes(value: unknown): string[] {
  if (!value) return [];
  try { return JSON.parse(String(value)); } catch { return []; }
}

export default async function PathwaysPage(){
  const stability=await apiGet<Term[]>("/api/pathways/stability");
  const all=await apiGet<Term[]>("/api/pathways?significant_only=true&limit=500");
  const stable=stability.filter((x:any)=>String(x.significant_all_thresholds).toLowerCase()==="true");
  const recurrent=stability.filter((x:any)=>Number(x.significant_thresholds_n)>=2 && String(x.significant_all_thresholds).toLowerCase()!=="true");

  return <>
    <div className="section-header" style={{marginBottom:18}}>
      <div>
        <div className="eyebrow">FUNCTIONAL LANDSCAPE</div>
        <h1 style={{margin:"4px 0 8px"}}>Pathway Explorer</h1>
        <div className="section-copy">GO Biological Process, Reactome, KEGG и CORUM. Главная задача этой страницы — отделить устойчивые функциональные сигналы от cut-off-sensitive enrichment.</div>
      </div>
    </div>

    <section className="kpi-strip" style={{marginTop:0}}>
      <div className="kpi"><strong>{stable.length}</strong><span>stable across 50/100/200</span></div>
      <div className="kpi"><strong>{recurrent.length}</strong><span>significant in ≥2 thresholds</span></div>
      <div className="kpi"><strong>{all.filter((x:any)=>Number(x.top_n)===50).length}</strong><span>significant at Top-50</span></div>
      <div className="kpi"><strong>{all.filter((x:any)=>Number(x.top_n)===200).length}</strong><span>significant at Top-200</span></div>
    </section>

    <section className="section">
      <div className="section-header"><div><h2>Threshold-stable signals</h2><div className="section-copy">Каждый термин показан один раз. Три ячейки отражают его значимость при Top-50, Top-100 и Top-200.</div></div></div>
      <div className="pathway-grid">
        {stable.map((x:any)=>{
          const union=[...new Set([...
            genes(x.intersecting_gene_symbols_top50),
            ...genes(x.intersecting_gene_symbols_top100),
            ...genes(x.intersecting_gene_symbols_top200),
          ])];
          return <article className="pathway-card" key={`${x.source}-${x.term_id}`}>
            <div className="eyebrow">{x.source} · {x.term_id}</div>
            <h3>{x.term_name}</h3>
            <div className="gene-cloud">{union.map((g:string)=><span className="gene-pill" key={g}>{g}</span>)}</div>
            <div className="pathway-thresholds">
              {[50,100,200].map(t=>{
                const on=String(x[`significant_top${t}`]).toLowerCase()==="true";
                return <div key={t} className={`threshold-cell ${on?"on":"off"}`}><b>Top-{t}</b><br/>{on?`p=${formatNumber(x[`p_value_adjusted_top${t}`],4)}`:"not significant"}</div>
              })}
            </div>
          </article>
        })}
      </div>
    </section>

    <section className="section">
      <div className="section-header"><div><h2>Signals significant in two thresholds</h2><div className="section-copy">Более широкие модули, которые воспроизводятся, но не проходят самый строгий критерий стабильности 3/3.</div></div></div>
      <div className="table-wrap"><table><thead><tr><th>Source</th><th>Term</th><th>Thresholds</th><th>Minimum adjusted p</th></tr></thead><tbody>{recurrent.slice(0,40).map((x:any)=><tr key={`${x.source}-${x.term_id}`}><td>{x.source}</td><td>{x.term_name}</td><td>{x.thresholds_significant}</td><td>{formatNumber(x.min_p_value_adjusted,6)}</td></tr>)}</tbody></table></div>
    </section>

    <section className="section">
      <div className="section-header"><div><h2>All significant enrichment observations</h2><div className="section-copy">Подробный слой: одна строка соответствует одному термину при конкретном Top-N.</div></div></div>
      <div className="table-wrap"><table><thead><tr><th>Top-N</th><th>Source</th><th>Term</th><th>Adjusted p</th><th>Intersection</th></tr></thead><tbody>{all.map((x:any,i:number)=><tr key={`${x.term_id}-${x.top_n}-${i}`}><td>{x.top_n}</td><td>{x.source}</td><td>{x.term_name}</td><td>{formatNumber(x.p_value_adjusted,6)}</td><td>{x.intersecting_gene_symbols_json}</td></tr>)}</tbody></table></div>
    </section>
  </>;
}
