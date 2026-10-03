import Link from "next/link";
import { apiGet } from "../../lib/api";
import ResearchTrail from "../../components/ResearchTrail";

type Comparison = { id:string; label:string; cancer_id:string; comparison:string; context_models_n?:number; comparator_models_n?:number; genes_analyzed_n?:number; depmap_release?:string; qc_status:string; context_definition?:string; comparator_definition?:string };

function statusLabel(status:string){
  const s=status.toUpperCase();
  if(s==="ERROR") return "Критическая ошибка";
  if(s==="WARNING") return "Есть ограничения";
  return "Проверки пройдены";
}

export default async function ComparisonsPage(){
  const items = await apiGet<Comparison[]>("/api/comparisons");
  return <>
    <div className="section-header" style={{marginBottom:18}}>
      <div>
        <div className="eyebrow">ШАГ 1</div>
        <h1 style={{margin:"4px 0 8px"}}>Опухолевые контексты</h1>
        <div className="section-copy">Выберите научный вопрос: какую молекулярную подгруппу опухолевых клеток мы сравниваем и с какой контрольной группой.</div>
      </div>
    </div>
    <ResearchTrail current={1} />
    <section className="callout section">
      <div className="eyebrow">КАК ЧИТАТЬ ЭТУ СТРАНИЦУ</div>
      <h3>Каждая карточка — отдельный сравнительный эксперимент</h3>
      <p className="section-copy">Мы спрашиваем: какие гены сильнее нужны клеткам с конкретным молекулярным признаком, чем похожим клеткам из группы сравнения? Клик по карточке откроет полногеномную карту зависимостей.</p>
    </section>
    <div className="grid cards section">
      {items.map(x => <Link key={x.id} href={`/comparisons/${encodeURIComponent(x.id)}`} className="card comparison-card">
        <div className="label">{x.cancer_id}</div>
        <h3>{x.label}</h3>
        <div className="comparison-groups">
          <div><span>Целевая группа</span><b>{x.context_definition || "Молекулярный контекст"}</b></div>
          <div className="vs">против</div>
          <div><span>Группа сравнения</span><b>{x.comparator_definition || "Контрольная группа"}</b></div>
        </div>
        <div className="comparison-meta"><span>Клеточные модели: <b>{x.context_models_n ?? "—"} / {x.comparator_models_n ?? "—"}</b></span><span className={`badge ${x.qc_status.toLowerCase()}`}>{statusLabel(x.qc_status)}</span></div>
      </Link>)}
    </div>
  </>;
}
