import Link from "next/link";
import CellModelTable from "../../../components/CellModelTable";
import { apiGet } from "../../../lib/api";
import styles from "./page.module.css";

type Context = Record<string, any>;
type ModelsPayload = { total:number; items:Record<string,any>[] };

function statusCopy(ctx: Context){
  if(ctx.analysis_available) return {cls:"ready", title:"Функциональный анализ доступен", text:"Для этого контекста уже рассчитаны полногеномные сравнения CRISPR-зависимостей."};
  return {cls:"pending", title:"Экспериментальные модели описаны", text:"Модели уже собраны, но полногеномное сравнение зависимостей для этого контекста ещё не опубликовано."};
}

function pct(value:any){
  const n=Number(value);
  if(!Number.isFinite(n)) return 0;
  return Math.max(0,Math.min(100,Math.round(n*100)));
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
  const patient=ctx.patient_layer || {};
  const genetics=ctx.model_genetics || {};

  return <>
    <div className="breadcrumbs">
      <Link href="/atlas">Атлас опухолей</Link><span>›</span><span>{ctx.organ_ru}</span><span>›</span><strong>{ctx.cancer_ru}</strong>
    </div>

    <section className="context-hero">
      <div className="context-hero-main">
        <div className="eyebrow">{ctx.organ_ru} · {ctx.id}</div>
        <h1>{ctx.cancer_ru}</h1>
        <div className="molecular-context-pill">Исследуемый подтип <b>{ctx.molecular_ru}</b></div>
        <p>От заболевания у пациентов — через молекулярно определённую группу экспериментальных моделей — к функциональным зависимостям и кандидатам в мишени.</p>
      </div>
      <div className={`context-status-panel ${status.cls}`}>
        <span className="context-status-dot"/>
        <div><strong>{status.title}</strong><p>{status.text}</p></div>
      </div>
    </section>

    <section className={styles.researchRoute} aria-label="Маршрут исследования">
      <div className={styles.routeStep}><span>1</span><div><small>Заболевание</small><b>Опухоль у пациентов</b></div></div>
      <i>→</i>
      <div className={styles.routeStep}><span>2</span><div><small>Стратификация</small><b>{ctx.molecular_ru}</b></div></div>
      <i>→</i>
      <div className={`${styles.routeStep} ${styles.routeCurrent}`}><span>3</span><div><small>Эксперимент</small><b>{ctx.models_n} клеточных моделей</b></div></div>
      <i>→</i>
      <div className={styles.routeStep}><span>4</span><div><small>Функциональный слой</small><b>CRISPR-зависимости</b></div></div>
    </section>

    <section className={styles.layerGrid}>
      <article className={`${styles.layerCard} ${styles.patientLayer}`}>
        <div className={styles.layerTop}>
          <span className={styles.layerNumber}>01</span>
          <span className={styles.pendingBadge}>источник не подключён</span>
        </div>
        <div className="eyebrow">ОПУХОЛЬ У ПАЦИЕНТОВ</div>
        <h2>{ctx.cancer_ru}</h2>
        <p>{patient.description}</p>
        <div className={styles.layerFooter}>
          <b>Что должно появиться здесь</b>
          <span>Частоты генетических изменений, сочетания мутаций и пациентские подгруппы из реальных опухолевых когорт.</span>
          <small>Планируемые источники: {(patient.planned_sources || []).join(" · ") || "—"}</small>
        </div>
      </article>

      <article className={`${styles.layerCard} ${styles.molecularLayer}`}>
        <div className={styles.layerTop}>
          <span className={styles.layerNumber}>02</span>
          <span className={styles.readyBadge}>доступно</span>
        </div>
        <div className="eyebrow">МОЛЕКУЛЯРНЫЙ ПОДТИП</div>
        <h2>{ctx.molecular_ru}</h2>
        <p>{ctx.context_definition || ctx.name}</p>
        <div className={styles.layerFooter}>
          <b>Зачем нужен этот уровень</b>
          <span>Он определяет биологический признак, по которому MCL формирует целевую группу клеточных моделей.</span>
        </div>
      </article>

      <article className={`${styles.layerCard} ${styles.modelLayer}`}>
        <div className={styles.layerTop}>
          <span className={styles.layerNumber}>03</span>
          <span className={styles.readyBadge}>доступно</span>
        </div>
        <div className="eyebrow">ЭКСПЕРИМЕНТАЛЬНЫЕ МОДЕЛИ</div>
        <h2>{ctx.models_n} клеточных линий</h2>
        <p>Каждая модель имеет собственный молекулярный фон. Именно на этих моделях измеряются CRISPR-зависимости.</p>
        <div className={styles.miniMetrics}>
          <div><b>{ctx.context_models_n}</b><span>целевая группа</span></div>
          <div><b>{ctx.comparator_models_n}</b><span>группа сравнения</span></div>
          <div><b>{sequencedPct}%</b><span>с мутационными данными</span></div>
        </div>
      </article>
    </section>

    <section className={styles.guardrail}>
      <div className={styles.guardrailIcon}>!</div>
      <div><b>Не смешиваем уровни данных</b><span>Частота мутации среди клеточных линий ниже — это характеристика доступного набора моделей DepMap. Она не показывает распространённость мутации среди пациентов с {ctx.cancer_ru}.</span></div>
    </section>

    <section className="section">
      <div className="section-header">
        <div>
          <div className="eyebrow">ГЕНЕТИКА ЭКСПЕРИМЕНТАЛЬНЫХ МОДЕЛЕЙ</div>
          <h2>Какие изменения встречаются в доступных клеточных линиях?</h2>
          <div className="section-copy">Этот слой помогает увидеть молекулярный фон моделей до интерпретации функциональных зависимостей. Для каждой конкретной линии можно открыть индивидуальный мутационный профиль.</div>
        </div>
      </div>

      {genetics.available ? <div className={styles.geneticsPanel}>
        <div className={styles.geneticsSummary}>
          <div><strong>{genetics.profiled_models_n ?? 0}</strong><span>моделей имеют записи в мутационном индексе</span></div>
          <div><strong>{genetics.functional_variants_n ?? 0}</strong><span>функционально значимых записей в выбранном наборе</span></div>
          <div><strong>{genetics.top_genes?.length ?? 0}</strong><span>генов показано в сводке</span></div>
        </div>
        <div className={styles.geneBars}>
          {(genetics.top_genes || []).slice(0,15).map((gene:any)=>{
            const percent=pct(gene.model_fraction);
            return <div className={styles.geneBarRow} key={gene.gene}>
              <Link href={`/genes/${encodeURIComponent(gene.gene)}`} className={styles.geneName}>{gene.gene}</Link>
              <div className={styles.barTrack}><span className={styles.barFill} style={{width:`${percent}%`}}/></div>
              <div className={styles.geneCount}><b>{gene.models_n}</b><span>моделей · {percent}%</span></div>
            </div>;
          })}
        </div>
        <p className={styles.dataNote}>{genetics.note}</p>
      </div> : <div className={styles.indexPrompt}>
        <div>
          <span className={styles.indexIcon}>DNA</span>
          <div><b>Полный профиль клеточных линий ещё не проиндексирован</b><p>{genetics.note}</p></div>
        </div>
        <code>{genetics.build_command}</code>
        <small>Команда один раз создаст компактный локальный индекс только для моделей MCL из уже загруженного DepMap OmicsSomaticMutations.csv.</small>
      </div>}
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
          <h2>Из каких линий получаются наши выводы?</h2>
          <div className="section-copy">Откройте конкретную линию, чтобы увидеть её индивидуальную генетику и понять, чем она отличается от других моделей того же опухолевого контекста.</div>
        </div>
      </div>
      <CellModelTable models={models.items as any[]}/>
    </section>

    <section className="section">
      <div className="section-header">
        <div>
          <div className="eyebrow">ФУНКЦИОНАЛЬНЫЕ СРАВНЕНИЯ</div>
          <h2>Какие вопросы уже рассчитаны?</h2>
          <div className="section-copy">После проверки состава моделей можно переходить к CRISPR-зависимостям. Каждая карточка ниже — отдельный научный вопрос, а не просто технический файл.</div>
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
      <b>Происхождение данных:</b> заболевание и подтип — конфигурация MCL и OncoTree; клеточные модели — DepMap; индивидуальные мутации — DepMap OmicsSomaticMutations после локальной индексации. Пациентская геномика пока не подключена и не подменяется статистикой клеточных линий. Внутренний идентификатор: {ctx.id}.
    </section>
  </>;
}
