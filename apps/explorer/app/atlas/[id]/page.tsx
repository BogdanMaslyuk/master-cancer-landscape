import Link from "next/link";
import CellModelTable from "../../../components/CellModelTable";
import { apiGet } from "../../../lib/api";

type Context = Record<string, any>;
type ModelsPayload = { total:number; items:Record<string,any>[] };

function statusCopy(ctx: Context){
  if(ctx.analysis_available) return {cls:"ready", title:"Полногеномный анализ доступен", text:"Можно переходить к уже рассчитанным сравнениям зависимостей."};
  return {cls:"pending", title:"Клеточные модели собраны", text:"Контекст уже описан в DepMap-аудите, но полногеномное сравнение для него пока не опубликовано в Explorer."};
}

export default async function AtlasContextPage({params}:{params:Promise<{id:string}>}){
  const {id}=await params;
  const cancerId=decodeURIComponent(id);
  const [ctx, models]=await Promise.all([
    apiGet<Context>(`/api/atlas/${encodeURIComponent(cancerId)}`),
    apiGet<ModelsPayload>(`/api/models?cancer_id=${encodeURIComponent(cancerId)}&limit=5000`),
  ]);
  const status=statusCopy(ctx);
  const sequencedPct=ctx.models_n ? Math.round((Number(ctx.sequenced_models_n||0)/Number(ctx.models_n))*100) : 0;

  return <>
    <div className="breadcrumbs">
      <Link href="/atlas">Атлас опухолей</Link><span>›</span><span>{ctx.organ_ru}</span><span>›</span><strong>{ctx.cancer_ru}</strong>
    </div>

    <section className="context-hero">
      <div className="context-hero-main">
        <div className="eyebrow">{ctx.organ_ru} · {ctx.id}</div>
        <h1>{ctx.cancer_ru}</h1>
        <div className="molecular-context-pill">Молекулярный контекст <b>{ctx.molecular_ru}</b></div>
        <p>{ctx.context_definition || ctx.name}</p>
      </div>
      <div className={`context-status-panel ${status.cls}`}>
        <span className="context-status-dot"/>
        <div><strong>{status.title}</strong><p>{status.text}</p></div>
      </div>
    </section>

    <section className="context-flow">
      <div className="context-flow-item done"><span>1</span><div><small>Орган</small><b>{ctx.organ_ru}</b></div></div>
      <i>→</i>
      <div className="context-flow-item done"><span>2</span><div><small>Опухоль</small><b>{ctx.cancer_ru}</b></div></div>
      <i>→</i>
      <div className="context-flow-item current"><span>3</span><div><small>Молекулярный контекст</small><b>{ctx.molecular_ru}</b></div></div>
      <i>→</i>
      <div className="context-flow-item"><span>4</span><div><small>Дальше</small><b>Клеточные линии и сравнение</b></div></div>
    </section>

    <section className="kpi-strip context-kpis">
      <div className="kpi"><strong>{ctx.models_n}</strong><span>клеточных моделей соответствуют опухолевому типу</span></div>
      <div className="kpi"><strong>{ctx.context_models_n}</strong><span>моделей в целевой молекулярной группе</span></div>
      <div className="kpi"><strong>{ctx.comparator_models_n}</strong><span>моделей назначены в базовую группу сравнения</span></div>
      <div className="kpi"><strong>{sequencedPct}%</strong><span>моделей имеют доступное мутационное профилирование</span></div>
    </section>

    <section className="section context-science-grid">
      <div className="card context-definition-card">
        <div className="eyebrow">ЦЕЛЕВАЯ ГРУППА</div>
        <h3>{ctx.molecular_ru}</h3>
        <p className="section-copy">{ctx.context_definition || "Определение целевой группы берётся из конфигурации MCL."}</p>
        <div className="context-count"><b>{ctx.context_models_n}</b><span>моделей</span></div>
      </div>
      <div className="versus-mark">VS</div>
      <div className="card context-definition-card comparator">
        <div className="eyebrow">БАЗОВАЯ ГРУППА СРАВНЕНИЯ</div>
        <h3>Контроль внутри того же заболевания</h3>
        <p className="section-copy">{ctx.comparator_definition || "Определение контрольной группы берётся из конфигурации MCL."}</p>
        <div className="context-count"><b>{ctx.comparator_models_n}</b><span>моделей</span></div>
      </div>
    </section>

    <section className="section">
      <div className="section-header">
        <div>
          <div className="eyebrow">КЛЕТОЧНЫЕ МОДЕЛИ</div>
          <h2>Какие линии формируют этот контекст?</h2>
          <div className="section-copy">Здесь виден состав выборки до перехода к результатам. Можно проверить название линии, DepMap ID, назначенную роль, доступность секвенирования и зарегистрированные варианты.</div>
        </div>
      </div>
      <CellModelTable models={models.items as any[]}/>
    </section>

    <section className="section">
      <div className="section-header">
        <div>
          <div className="eyebrow">НАУЧНЫЕ СРАВНЕНИЯ</div>
          <h2>Что уже можно исследовать внутри этого контекста?</h2>
          <div className="section-copy">Каждое сравнение — отдельный научный вопрос. Размеры групп для конкретного сравнения могут отличаться от базового аудита из-за дополнительных критериев отбора.</div>
        </div>
      </div>
      {ctx.comparisons?.length ? <div className="analysis-choice-grid">
        {ctx.comparisons.map((cmp:any)=><Link key={cmp.id} href={`/comparisons/${encodeURIComponent(cmp.id)}`} className="analysis-choice-card">
          <div className="analysis-choice-top"><span className="analysis-state ready">рассчитано</span><small>DepMap {cmp.depmap_release || "—"}</small></div>
          <h3>{cmp.label}</h3>
          <div className="analysis-question"><span>Целевая группа</span><b>{cmp.context_definition || ctx.context_definition}</b></div>
          <div className="analysis-vs">сравниваем с</div>
          <div className="analysis-question"><span>Группа сравнения</span><b>{cmp.comparator_definition || ctx.comparator_definition}</b></div>
          <div className="analysis-choice-footer"><span>{cmp.context_models_n ?? "—"} vs {cmp.comparator_models_n ?? "—"} моделей</span><b>Открыть карту зависимостей →</b></div>
        </Link>)}
      </div> : <div className="empty-state atlas-empty">
        <strong>Полногеномное сравнение пока не рассчитано</strong>
        <span>Клеточные модели и молекулярный контекст уже доступны для изучения. Этот контекст можно подключить к следующей волне анализа MCL.</span>
      </div>}
    </section>

    <section className="section technical-note context-provenance">
      <b>Происхождение данных:</b> иерархия заболевания — конфигурация MCL и OncoTree-аннотации DepMap; состав клеточных моделей и мутационный статус — локальный аудит DepMap. Внутренний идентификатор контекста: {ctx.id}.
    </section>
  </>;
}
