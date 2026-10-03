import { apiGet, formatNumber } from "../../lib/api";

type Term=Record<string,any>;
export default async function PathwaysPage(){
  const stable=await apiGet<Term[]>("/api/pathways?stable_only=true&significant_only=true&limit=500");
  const all=await apiGet<Term[]>("/api/pathways?significant_only=true&limit=500");
  return <><h1>Pathway Explorer</h1><p className="muted">GO Biological Process, Reactome, KEGG и CORUM. Устойчивость оценивается между Top-50/100/200.</p>
    <section className="section"><h2>Stable across all thresholds</h2><div className="grid cards">{stable.map((x:any)=><div className="card" key={`${x.source}-${x.term_id}-${x.top_n}`}><div className="label">{x.source} · Top-{x.top_n}</div><h3>{x.term_name}</h3><div>p<sub>adj</sub>: <b>{formatNumber(x.p_value_adjusted,6)}</b></div><div className="muted" style={{marginTop:8}}>{x.intersecting_gene_symbols_json}</div></div>)}</div></section>
    <section className="section"><h2>All significant terms</h2><div className="table-wrap"><table><thead><tr><th>Top-N</th><th>Source</th><th>Term</th><th>Adjusted p</th><th>Intersection</th></tr></thead><tbody>{all.map((x:any,i:number)=><tr key={`${x.term_id}-${x.top_n}-${i}`}><td>{x.top_n}</td><td>{x.source}</td><td>{x.term_name}</td><td>{formatNumber(x.p_value_adjusted,6)}</td><td>{x.intersecting_gene_symbols_json}</td></tr>)}</tbody></table></div></section>
  </>;
}
