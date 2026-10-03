import Link from "next/link";
import DependencyScatter from "../../../components/DependencyScatter";
import { apiGet, formatNumber } from "../../../lib/api";

type Gene = Record<string, any>;
type GenePage = { comparison_id:string; total:number; items:Gene[] };
type Comparison = Record<string, any>;

export default async function ComparisonPage({ params }: { params: Promise<{id:string}> }) {
  const {id} = await params;
  const comparisonId = decodeURIComponent(id);
  const [meta, genes, scatter] = await Promise.all([
    apiGet<Comparison>(`/api/comparisons/${encodeURIComponent(comparisonId)}`),
    apiGet<GenePage>(`/api/comparisons/${encodeURIComponent(comparisonId)}/genes?page_size=100&exclude_broad=true&exclude_low_sample=true`),
    apiGet<GenePage>(`/api/comparisons/${encodeURIComponent(comparisonId)}/genes?page_size=2000&exclude_broad=true&exclude_low_sample=true`),
  ]);

  const stableN = scatter.items.filter((g:any)=>String(g.present_all_thresholds).toLowerCase()==="true").length;

  return <>
    <section className="hero">
      <div className="eyebrow" style={{color:"#bad5f1"}}>{meta.cancer_id} · DepMap {meta.depmap_release}</div>
      <h1>{meta.label}</h1>
      <p>{meta.context_definition || "Target context"} vs {meta.comparator_definition || "Comparator"}</p>
    </section>

    <section className="kpi-strip">
      <div className="kpi"><strong>{meta.context_models_n ?? "—"}</strong><span>target models</span></div>
      <div className="kpi"><strong>{meta.comparator_models_n ?? "—"}</strong><span>comparator models</span></div>
      <div className="kpi"><strong>{formatNumber(meta.genes_analyzed_n,0)}</strong><span>genes analyzed</span></div>
      <div className="kpi"><strong>{stableN}</strong><span>stable genes in displayed eligible set</span></div>
    </section>

    <section className="section split">
      <div>
        <div className="section-header">
          <div>
            <h2>Dependency landscape</h2>
            <div className="section-copy">Каждая точка — ген. По оси X медианный Gene Effect в контроле, по оси Y — в целевом контексте. Точки ниже диагонали сильнее необходимы целевой группе.</div>
          </div>
        </div>
        <DependencyScatter points={scatter.items} />
      </div>
      <aside className="callout">
        <div className="eyebrow">Как читать</div>
        <h3>Контекстная зависимость ≠ лекарственная мишень</h3>
        <p className="section-copy">На этом этапе график показывает генетическую CRISPR-зависимость. Broad dependency и low sample size исключены из визуального слоя, но остаются доступны в исходных данных и QC.</p>
        <div className="metric-row">
          <span className={`badge ${(meta.qc_status || "PASS").toLowerCase()}`}>QC {meta.qc_status || "PASS"}</span>
          <span className="metric-chip">displayed {formatNumber(scatter.items.length,0)}</span>
          <span className="metric-chip">eligible total {formatNumber(scatter.total,0)}</span>
        </div>
      </aside>
    </section>

    <section className="section">
      <div className="section-header">
        <div>
          <h2>Top genome-wide dependencies</h2>
          <div className="section-copy">Первые 100 генов после исключения broad dependency и low sample size. Сортировка — по наиболее отрицательной Δ Gene Effect.</div>
        </div>
      </div>
      <div className="table-wrap"><table><thead><tr><th>Gene</th><th>Δ Gene Effect</th><th>Target median</th><th>Comparator median</th><th>Cliff's δ</th><th>q-value</th><th>Robustness</th></tr></thead><tbody>
      {genes.items.map((g:any) => <tr key={g.gene_symbol}><td><Link href={`/genes/${g.gene_symbol}`} style={{color:"var(--accent)",fontWeight:900}}>{g.gene_symbol}</Link></td><td>{formatNumber(g.delta_gene_effect)}</td><td>{formatNumber(g.context_median_gene_effect)}</td><td>{formatNumber(g.comparator_median_gene_effect)}</td><td>{formatNumber(g.cliffs_delta)}</td><td>{formatNumber(g.q_value,4)}</td><td>{String(g.present_all_thresholds).toLowerCase()==="true"?<span className="badge stable">stable 50/100/200</span>:"—"}</td></tr>)}
      </tbody></table></div>
    </section>
  </>;
}
