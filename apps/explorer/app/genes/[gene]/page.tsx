import { apiGet, formatNumber } from "../../../lib/api";

type Payload = { identity:Record<string,any>; stability:Record<string,any>|null; comparisons:Record<string,any>[]; pathways:Record<string,any>[] };
export default async function GenePage({params}:{params:Promise<{gene:string}>}){
  const {gene}=await params;
  const data=await apiGet<Payload>(`/api/genes/${encodeURIComponent(gene)}`);
  return <>
    <section className="hero"><div className="muted" style={{color:"#bad5f1"}}>GENE</div><h1>{data.identity.gene_symbol || gene}</h1><p>Контекстные CRISPR-зависимости, устойчивость и связанные enrichment-сигналы.</p></section>
    <section className="grid cards">
      <div className="card"><div className="label">Stability</div><div className="value" style={{fontSize:22}}>{data.stability?.present_all_thresholds ? "Stable recurrent" : "Context-dependent"}</div></div>
      <div className="card"><div className="label">HGNC</div><div className="value" style={{fontSize:22}}>{data.identity.hgnc_id || data.identity.HGNC_ID || "—"}</div></div>
      <div className="card"><div className="label">Ensembl</div><div className="value" style={{fontSize:18}}>{data.identity.ensembl_gene_id || data.identity.ensembl_id || "—"}</div></div>
      <div className="card"><div className="label">Thresholds</div><div className="value" style={{fontSize:22}}>{data.stability?.thresholds_present || "—"}</div></div>
    </section>
    <section className="section"><h2>Comparison matrix</h2><div className="table-wrap"><table><thead><tr><th>Comparison</th><th>Δ Gene Effect</th><th>Context</th><th>Comparator</th><th>Cliff's delta</th><th>q-value</th><th>Broad</th><th>Low sample</th></tr></thead><tbody>{data.comparisons.map((r:any)=><tr key={r.comparison_id}><td>{r.comparison_label}</td><td>{formatNumber(r.delta_gene_effect)}</td><td>{formatNumber(r.context_median_gene_effect)}</td><td>{formatNumber(r.comparator_median_gene_effect)}</td><td>{formatNumber(r.cliffs_delta)}</td><td>{formatNumber(r.q_value,4)}</td><td>{String(r.broad_dependency_warning)}</td><td>{String(r.low_sample_size)}</td></tr>)}</tbody></table></div></section>
    <section className="section"><h2>Pathway associations</h2><div className="table-wrap"><table><thead><tr><th>Top-N</th><th>Source</th><th>Term</th><th>Adjusted p</th><th>Intersection</th></tr></thead><tbody>{data.pathways.filter((p:any)=>String(p.significant).toLowerCase()==="true").slice(0,80).map((p:any,i:number)=><tr key={`${p.term_id}-${p.top_n}-${i}`}><td>{p.top_n}</td><td>{p.source}</td><td>{p.term_name}</td><td>{formatNumber(p.p_value_adjusted,5)}</td><td>{p.intersecting_gene_symbols_json}</td></tr>)}</tbody></table></div></section>
  </>;
}
