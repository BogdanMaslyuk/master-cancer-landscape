import Link from "next/link";
import { apiGet, formatNumber } from "../lib/api";

type Summary = {
  genes_analyzed_n: number;
  comparisons_n: number;
  stable_recurrent_genes_n: number;
  stable_pathways_n: number;
  thresholds: number[];
  per_threshold: Record<string, { query_genes_n?: number; significant_terms_n?: number }>;
  data_release?: string;
};

export default async function OverviewPage() {
  const s = await apiGet<Summary>("/api/summary");
  return <>
    <section className="hero">
      <div className="muted" style={{color:"#bad5f1", fontWeight:700}}>MCL EXPLORER v0.1</div>
      <h1>Карта функциональных зависимостей опухоли</h1>
      <p>Интерактивный слой над Master Cancer Landscape: от genome-wide CRISPR-зависимостей DepMap до устойчивых recurrent-генов и функциональных модулей.</p>
    </section>
    <section className="grid cards">
      <div className="card"><div className="label">Генов на genome-wide сравнение</div><div className="value">{formatNumber(s.genes_analyzed_n,0)}</div></div>
      <div className="card"><div className="label">Основных сравнений</div><div className="value">{s.comparisons_n}</div></div>
      <div className="card"><div className="label">Stable recurrent genes</div><div className="value">{s.stable_recurrent_genes_n}</div></div>
      <div className="card"><div className="label">Stable pathway signals</div><div className="value">{s.stable_pathways_n}</div></div>
    </section>
    <section className="section">
      <h2>Воронка анализа</h2>
      <div className="funnel">
        <div className="step"><b>Genome-wide</b><br/>≈ {formatNumber(s.genes_analyzed_n,0)} генов</div>
        <div className="step"><b>Eligibility</b><br/>QC + broad dependency</div>
        <div className="step"><b>Top-N</b><br/>50 / 100 / 200</div>
        <div className="step"><b>Recurrent</b><br/>≥2 из 3 сравнений</div>
        <div className="step"><b>Stable</b><br/>{s.stable_recurrent_genes_n} генов</div>
        <div className="step"><b>Modules</b><br/>{s.stable_pathways_n} устойчивых сигнала</div>
      </div>
    </section>
    <section className="section">
      <h2>Sensitivity Top-N</h2>
      <div className="grid cards">
        {s.thresholds.map(t => {
          const d = s.per_threshold?.[String(t)] || {};
          return <div className="card" key={t}><div className="label">Top-{t}</div><div className="value">{d.query_genes_n ?? "—"}</div><div className="muted">recurrent genes · {d.significant_terms_n ?? "—"} significant terms</div></div>
        })}
      </div>
    </section>
    <section className="section card">
      <b>Текущий релиз данных:</b> {s.data_release || "—"}. <Link href="/comparisons" style={{color:"var(--accent)", fontWeight:700}}>Перейти к сравнению →</Link>
    </section>
  </>;
}
