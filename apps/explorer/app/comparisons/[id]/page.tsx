import Link from "next/link";
import { apiGet, formatNumber } from "../../../lib/api";

type Gene = Record<string, any>;
type GenePage = { comparison_id:string; total:number; items:Gene[] };
type Comparison = Record<string, any>;

export default async function ComparisonPage({ params }: { params: Promise<{id:string}> }) {
  const {id} = await params;
  const comparisonId = decodeURIComponent(id);
  const [meta, genes] = await Promise.all([
    apiGet<Comparison>(`/api/comparisons/${encodeURIComponent(comparisonId)}`),
    apiGet<GenePage>(`/api/comparisons/${encodeURIComponent(comparisonId)}/genes?page_size=100&exclude_broad=true&exclude_low_sample=true`)
  ]);
  return <>
    <div className="hero"><div className="muted" style={{color:"#bad5f1"}}>{meta.cancer_id} · {meta.depmap_release}</div><h1>{meta.label}</h1><p>{meta.context_definition || "Target context"} vs {meta.comparator_definition || "Comparator"}</p></div>
    <section className="grid cards">
      <div className="card"><div className="label">Context models</div><div className="value">{meta.context_models_n ?? "—"}</div></div>
      <div className="card"><div className="label">Comparator models</div><div className="value">{meta.comparator_models_n ?? "—"}</div></div>
      <div className="card"><div className="label">Genes analyzed</div><div className="value">{formatNumber(meta.genes_analyzed_n,0)}</div></div>
      <div className="card"><div className="label">QC</div><div style={{marginTop:12}}><span className={`badge ${(meta.qc_status || "PASS").toLowerCase()}`}>{meta.qc_status || "PASS"}</span></div></div>
    </section>
    <section className="section"><h2>Top genome-wide dependencies</h2><p className="muted">Показаны первые 100 после исключения broad dependency и low sample size.</p>
      <div className="table-wrap"><table><thead><tr><th>Gene</th><th>Δ Gene Effect</th><th>Context median</th><th>Comparator median</th><th>Cliff's delta</th><th>q-value</th><th>Stable</th></tr></thead><tbody>
      {genes.items.map((g:any) => <tr key={g.gene_symbol}><td><Link href={`/genes/${g.gene_symbol}`} style={{color:"var(--accent)",fontWeight:800}}>{g.gene_symbol}</Link></td><td>{formatNumber(g.delta_gene_effect)}</td><td>{formatNumber(g.context_median_gene_effect)}</td><td>{formatNumber(g.comparator_median_gene_effect)}</td><td>{formatNumber(g.cliffs_delta)}</td><td>{formatNumber(g.q_value,4)}</td><td>{String(g.present_all_thresholds).toLowerCase()==="true"?"● stable":"—"}</td></tr>)}
      </tbody></table></div>
    </section>
  </>;
}
