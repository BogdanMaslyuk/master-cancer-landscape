import Link from "next/link";
import AtlasIcon from "../../components/AtlasIcon";
import { apiGet } from "../../lib/api";

export const dynamic = "force-dynamic";

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
        <p>Начните с заболевания. MCL последовательно разделяет данные опухолей у пациентов, молекулярный подтип и конкретные экспериментальные модели — и только затем ведёт к функциональным зависимостям.</p>
      </div>
      <div className="atlas-route-card">
        <span>Путь исследования</span>
        <strong>Опухоль → молекулярный подтип → клеточные модели → CRISPR-зависимости → кандидаты</strong>
      </div>
    </section>

    <section className="kpi-strip atlas-kpis">
      <div className="kpi"><strong>{atlas.organs_n}</strong><span>органа / системы представлены</span></div>
      <div className="kpi"><strong>{atlas.contexts_n}</strong><span>молекулярных контекста настроено</span></div>
      <div className="kpi"><strong>{atlas.models_n}</strong><span>уникальных клеточных моделей в аудите</span></div>
      <div className="kpi"><strong>{atlas.analyses_n}</strong><span>полногеномных сравнений доступны</span></div>
    </section>

    <section className="section callout atlas-callout">
      <div className="eyebrow">КАК УСТРОЕНЫ ДАННЫЕ</div>
      <h3>Пациентская опухоль и клеточная линия — не один уровень</h3>
      <p className="section-copy">Частота мутации в опухолях пациентов, критерий молекулярного подтипа и генетика конкретной клеточной линии отвечают на разные вопросы. MCL показывает их раздельно, чтобы не переносить свойства небольшого набора моделей на всех пациентов.</p>
    </section>

    <section className="section">
      <div className="section-header">
        <div>
          <div className="eyebrow">ШАГ 1</div>
          <h2>Выберите орган и заболевание</h2>
          <div className="section-copy">Внутри каждого раздела показаны молекулярные контексты, для которых уже описаны экспериментальные модели. Пациентская геномика будет подключаться отдельным доказательным слоем и не подменяется статистикой DepMap.</div>
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
                <span>Подтип: {ctx.molecular_ru}</span>
              </div>
              <div className="context-link-meta">
                <small>{ctx.models_n} клеточных моделей</small>
                <span className={`analysis-state ${ctx.analysis_available?"ready":"pending"}`}>{ctx.analysis_available?"функциональный анализ доступен":"модели собраны"}</span>
              </div>
              <b className="context-arrow">→</b>
            </Link>)}
          </div>
        </article>)}
      </div>
    </section>

    <section className="section callout atlas-callout">
      <div className="eyebrow">ЗАЧЕМ НУЖНЫ КЛЕТОЧНЫЕ МОДЕЛИ</div>
      <h3>Групповая зависимость должна быть объяснима через конкретные линии</h3>
      <p className="section-copy">Даже модели с одной определяющей мутацией отличаются сопутствующей генетикой. Поэтому из любой групповой зависимости MCL должен позволять спуститься к клеточным линиям и проверить индивидуальный молекулярный фон каждой из них.</p>
    </section>
  </>;
}
