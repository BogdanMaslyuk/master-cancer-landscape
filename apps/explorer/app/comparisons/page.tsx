import Link from "next/link";
import { apiGet } from "../../lib/api";

type Comparison = { id:string; label:string; cancer_id:string; comparison:string; context_models_n?:number; comparator_models_n?:number; genes_analyzed_n?:number; depmap_release?:string; qc_status:string; context_definition?:string; comparator_definition?:string };

export default async function ComparisonsPage(){
  const items = await apiGet<Comparison[]>("/api/comparisons");
  return <>
    <h1>Cancer Context Explorer</h1>
    <p className="muted">Выберите молекулярное сравнение и перейдите к genome-wide таблице зависимостей.</p>
    <div className="grid cards">
      {items.map(x => <Link key={x.id} href={`/comparisons/${encodeURIComponent(x.id)}`} className="card">
        <div className="label">{x.cancer_id} · {x.comparison}</div>
        <h3>{x.label}</h3>
        <p className="muted">{x.context_definition || "Target context"}<br/>vs<br/>{x.comparator_definition || "Comparator"}</p>
        <div>Models: <b>{x.context_models_n ?? "—"}</b> vs <b>{x.comparator_models_n ?? "—"}</b></div>
        <div style={{marginTop:8}}><span className={`badge ${x.qc_status.toLowerCase()}`}>{x.qc_status}</span></div>
      </Link>)}
    </div>
  </>;
}
