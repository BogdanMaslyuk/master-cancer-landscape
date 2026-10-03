import Link from "next/link";
import AtlasIcon from "../../components/AtlasIcon";
import { apiGet } from "../../lib/api";

type Context = Record<string, any>;
type Organ = { id:string; name_ru:string; name_en?:string; icon:string; contexts_n:number; models_n:number; analyses_n:number; contexts:Context[] };
type Atlas = { organs_n:number; contexts_n:number; models_n:number; analyses_n:number; organs:Organ[] };

export default async function AtlasPage(){
  const atlas = await apiGet<Atlas>("/api/atlas");
  return <>
    <section className="page-intro atlas-intro">
      <div>
        <div className="eyebrow">ГЛАВНАЯ ТОЧКА ВХОДА</div>
        <h1>Атлас опухолей</h1>
        <p>Начните с заболевания, а не с технического сравнения. Выберите орган, тип опухоли и молекулярный контекст — затем изучите доступные клеточные линии и только после этого переходите к генетическим зависимостям.</p>
      </div>
      <div className="atlas-route-card">
        <span>Путь исследования</span>
        <strong>Орган → опухоль → молекулярный контекст → клеточные линии → сравнение</strong>
      </div>
    </section>

    <section className="kpi-strip atlas-kpis">
      <div className="kpi"><strong>{atlas.organs_n}</strong><span>органа / системы представлены</span></div>
      <div className="kpi"><strong>{atlas.contexts_n}</strong><span>молекулярных контекста настроено</span></div>
      <div className="kpi"><strong>{atlas.models_n}</strong><span>уникальных клеточных моделей в аудите</span></div>
      <div className="kpi"><strong>{atlas.analyses_n}</strong><span>полногеномных сравнений доступны</span></div>
    </section>

    <section className="section">
      <div className="section-header">
        <div>
          <div className="eyebrow">ШАГ 1</div>
          <h2>Выберите орган или систему</h2>
          <div className="section-copy">Внутри каждого органа показаны только те опухолевые и молекулярные контексты, которые уже заведены в MCL. Контекст может иметь клеточные модели, даже если полногеномный анализ для него ещё не выполнен.</div>
        </div>
      </div>

      <div className="organ-grid">
        {atlas.organs.map((organ) => <article className="organ-card" key={organ.id}>
          <div className="organ-card-head">
            <div className={`organ-icon ${organ.icon}`}><AtlasIcon name={organ.icon}/></div>
            <div>
              <div className="eyebrow">ОРГАН / СИСТЕМА</div>
              <h3>{organ.name_ru}</h3>
              {organ.name_en && <div className="technical-note">OncoTree lineage: {organ.name_en}</div>}
            </div>
          </div>
          <div className="organ-metrics">
            <span><b>{organ.contexts_n}</b> контекст</span>
            <span><b>{organ.models_n}</b> моделей</span>
            <span><b>{organ.analyses_n}</b> анализов</span>
          </div>
          <div className="organ-contexts">
            {organ.contexts.map((ctx:any)=><Link href={`/atlas/${ctx.id}`} className="organ-context-link" key={ctx.id}>
              <div>
                <strong>{ctx.cancer_ru}</strong>
                <span>{ctx.molecular_ru}</span>
              </div>
              <div className="context-link-meta">
                <small>{ctx.models_n} клеточных моделей</small>
                <span className={`analysis-state ${ctx.analysis_available?"ready":"pending"}`}>{ctx.analysis_available?"анализ доступен":"модели собраны"}</span>
              </div>
              <b className="context-arrow">→</b>
            </Link>)}
          </div>
        </article>)}
      </div>
    </section>

    <section className="section callout atlas-callout">
      <div className="eyebrow">ПОЧЕМУ ТАК</div>
      <h3>Клеточная модель — не техническая деталь, а основа интерпретации</h3>
      <p className="section-copy">Одна и та же зависимость может быть специфична для ткани, мутации или конкретного набора моделей. Поэтому MCL теперь сначала показывает, какие именно линии формируют исследуемую группу и чем они отличаются от контроля.</p>
    </section>
  </>;
}
