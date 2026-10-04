import Link from "next/link";
import { apiGet, formatNumber } from "../../../lib/api";
import GeneDependencyLandscape from "./GeneDependencyLandscape";
import GeneInsightCharts from "./GeneInsightCharts";
import styles from "./GeneWorkbench.module.css";

export const dynamic = "force-dynamic";

type DescriptivePayload={available:boolean;items:Record<string,any>[];guardrail?:string};
type CorrelationPayload={available:boolean;relationships:Record<string,any>[];points:Record<string,any>[];statistics_engine?:string};
type MutationPayload={available:boolean;deferred?:boolean;profiled_models_n?:number;items:Record<string,any>[];statistics_engine?:string;note?:string};
type BasePayload = {
  identity:Record<string,any>;
  summary?:Record<string,any>;
  stability?:Record<string,any>;
  comparisons?:Record<string,any>[];
  pathways?:Record<string,any>[];
  model_insights?:{descriptive_contexts?:DescriptivePayload;correlations?:CorrelationPayload;mutation_associations?:MutationPayload};
};
type ContextPayload = { gene_symbol:string; total:number; items:Record<string,any>[] };
type ModelsPayload = {gene_symbol:string;available:boolean;available_layers:Record<string,boolean>;total:number;page:number;page_size:number;pages?:number;items:Record<string,any>[];guardrails?:Record<string,string>;note?:string};
type AnnotationPayload = {
  gene_symbol:string;status:string;mcl_domains:string[];subdomains:string[];protein_classes:string[];compartments:string[];hallmarks:string[];
  mcl_domain_details?:Record<string,any>[];protein_class_details?:Record<string,any>[];compartment_details?:Record<string,any>[];hallmark_details?:Record<string,any>[];
  formal_annotations:Record<string,any>[];coverage?:{annotated_genes_n:number;gene_universe_n:number;coverage_fraction:number|null;status:string;reason:string};reference_available?:boolean;
};
type LandscapePayload={available:boolean;summary:Record<string,any>;cancers:Record<string,any>[];interpretation?:Record<string,string>};

function yes(value:unknown){ return String(value).toLowerCase()==="true"; }
function sourceDetails(items:Record<string,any>[]|undefined){ return (items||[]).slice(0,16); }
function pct(value:unknown,digits=0){const n=Number(value);return Number.isFinite(n)?`${(n*100).toFixed(digits)}%`:"—";}

export default async function GenePage({params}:{params:Promise<{gene:string}>}){
  const {gene}=await params;
  const symbol=decodeURIComponent(gene).trim().toUpperCase();
  const [base,landscape,contexts,models,annotations]=await Promise.all([
    apiGet<BasePayload>(`/api/genes/${encodeURIComponent(symbol)}`),
    apiGet<LandscapePayload>(`/api/genes/${encodeURIComponent(symbol)}/dependency-landscape`),
    apiGet<ContextPayload>(`/api/genes/${encodeURIComponent(symbol)}/contexts`),
    apiGet<ModelsPayload>(`/api/genes/${encodeURIComponent(symbol)}/models?page_size=100&sort_by=gene_effect&sort_order=asc`),
    apiGet<AnnotationPayload>(`/api/genes/${encodeURIComponent(symbol)}/annotations`),
  ]);

  const summary=landscape.summary||{};
  const stable=yes(base.stability?.present_all_thresholds);
  const formal=(annotations.formal_annotations||[]).slice(0,30);
  const domainDetails=sourceDetails(annotations.mcl_domain_details);
  const classDetails=sourceDetails(annotations.protein_class_details);
  const compartmentDetails=sourceDetails(annotations.compartment_details);
  const hallmarkDetails=sourceDetails(annotations.hallmark_details);
  const aliases=(base.identity?.aliases||[]) as string[];
  const correlations=base.model_insights?.correlations;
  const mutationAssociations=base.model_insights?.mutation_associations;
  const topModels=(models.items||[]).slice(0,30);

  return <>
    <div className="breadcrumbs"><Link href="/genes">Гены и мишени</Link><span>›</span><strong>{symbol}</strong></div>

    <section className={styles.hero}>
      <div>
        <div className="eyebrow" style={{color:"#b9d5ee"}}>ГЕН · ФУНКЦИЯ → ОПУХОЛЬ → МОДЕЛЬ</div>
        <h1>{symbol}</h1>
        {base.identity?.gene_name&&<div style={{fontSize:18,fontWeight:700,color:"#e5eff9",marginTop:-4,marginBottom:9}}>{base.identity.gene_name}</div>}
        <p>Карточка отвечает на главный вопрос: где потеря функции {symbol} сильнее всего снижает жизнеспособность опухолевых клеток и какие молекулярные признаки могут быть с этим связаны.</p>
        {aliases.length>0&&<p style={{marginTop:8,fontSize:12}}>Синонимы: {aliases.slice(0,8).join(", ")}{aliases.length>8?" …":""}</p>}
      </div>
      <div className={styles.heroMeta}>
        <div><span>HGNC</span><b>{base.identity?.hgnc_id||"—"}</b></div>
        <div><span>Ensembl</span><b>{base.identity?.ensembl_gene_id||"—"}</b></div>
        <div><span>UniProt</span><b>{base.identity?.uniprot_id||"не подключено"}</b></div>
        <div><span>Функциональный класс</span><b>{annotations.protein_classes?.[0]||annotations.mcl_domains?.[0]||"пока не определён"}</b></div>
      </div>
    </section>

    <section className={styles.metrics}>
      <div className={styles.metric}><span>Зависимые CRISPR-модели</span><strong>{summary.dependent_models_n??0}/{summary.models_n??0}</strong><small>{pct(summary.dependent_fraction,1)} при Gene Effect ≤ −0,5</small></div>
      <div className={styles.metric}><span>Тип зависимости</span><strong className={styles.metricText}>{summary.dependency_type_ru||"—"}</strong></div>
      <div className={styles.metric}><span>Где выражена сильнее всего</span><strong className={styles.metricText}>{summary.best_cancer_name||"—"}</strong><small>{summary.best_cancer_organ_ru||""}</small></div>
      <div className={styles.metric}><span>Специфичность</span><strong className={styles.metricText}>{summary.specificity_label_ru||"—"}</strong><small>{Number.isFinite(Number(summary.specificity_score))?`оценка ${(Number(summary.specificity_score)*100).toFixed(0)} п.п.`:""}</small></div>
      <div className={styles.metric}><span>Глобальная медиана Gene Effect</span><strong>{formatNumber(summary.global_median_gene_effect)}</strong><small>по всем доступным моделям</small></div>
    </section>

    <section className={styles.section}>
      <div className={styles.sectionHeader}>
        <div><div className="eyebrow">ГЕН → ОПУХОЛИ</div><h2>Где опухолевые клетки зависят от {symbol}?</h2><p>Полный профиль строится по всем CRISPR-моделям Атласа, а не только по заранее заданным сравнительным контекстам.</p></div>
        <span className={styles.badge}>{landscape.cancers?.length||0} опухолевых групп</span>
      </div>
      <GeneDependencyLandscape symbol={symbol} payload={landscape}/>
    </section>

    <section className={styles.section}>
      <div className={styles.sectionHeader}>
        <div><div className="eyebrow">ГЕН → КЛЕТОЧНЫЕ МОДЕЛИ</div><h2>Какие модели сильнее всего зависят от {symbol}?</h2><p>Показаны наиболее отрицательные Gene Effect. Это помогает перейти от опухолевой группы к конкретной экспериментальной модели и её генетическому фону.</p></div>
        <div className={styles.layerState}><span className={models.available_layers?.gene_effect?styles.on:""}>CRISPR</span><span className={models.available_layers?.expression?styles.on:""}>РНК</span><span className={models.available_layers?.copy_number?styles.on:""}>число копий</span></div>
      </div>
      {models.available&&topModels.length?<div className={styles.tableWrap}><table className={styles.table}>
        <thead><tr><th>Модель</th><th>Опухоль</th><th>Gene Effect</th><th>РНК</th><th>Число копий</th></tr></thead>
        <tbody>{topModels.map((row:any)=><tr key={row.model_id}>
          <td><Link className={styles.modelLink} href={`/models/${encodeURIComponent(row.model_id)}`}>{row.cell_line_name||row.model_id}</Link>{row.cell_line_name&&<span className={styles.modelName}>{row.model_id}</span>}</td>
          <td>{row.oncotree_subtype||row.oncotree_primary_disease||row.oncotree_lineage||"—"}</td>
          <td className={Number(row.gene_effect)<0?styles.negative:""}><b>{formatNumber(row.gene_effect)}</b></td>
          <td>{formatNumber(row.expression)}</td><td>{formatNumber(row.copy_number)}</td>
        </tr>)}</tbody>
      </table></div>:<div className={styles.empty}>{models.note||`Для ${symbol} данные по моделям пока недоступны.`}</div>}
      {models.total>topModels.length&&<div className={styles.sectionNote}>Показаны 30 моделей с наиболее сильной зависимостью из {models.total}. Полный список доступен через API; следующим этапом добавим интерактивное раскрытие без перегрузки страницы.</div>}
    </section>

    <section className={styles.section}>
      <div className={styles.sectionHeader}>
        <div><div className="eyebrow">ЧТО МОЖЕТ ОБЪЯСНЯТЬ ЗАВИСИМОСТЬ?</div><h2>Связана ли зависимость с экспрессией РНК или числом копий?</h2><p>Знак ρ автоматически переводится в биологический смысл с учётом того, что более отрицательный Gene Effect означает более сильную зависимость. Отдельно показываем силу связи и её статистическую уверенность.</p></div>
        <span className={styles.badge}>{correlations?.statistics_engine||"нет данных"}</span>
      </div>
      {correlations?.available?<GeneInsightCharts points={correlations.points||[]} relationships={correlations.relationships||[]}/>:<div className={styles.empty}>Недостаточно совместных данных для оценки связи.</div>}
    </section>

    <section className={styles.section}>
      <div className={styles.sectionHeader}>
        <div><div className="eyebrow">ФУНКЦИЯ ГЕНА</div><h2>Что делает {symbol} в клетке?</h2><p>Функциональные категории нужны для биологической интерпретации зависимости. Все проекции сохраняют связь с исходными аннотациями.</p></div>
        <span className={styles.badge}>{annotations.reference_available?"справочный слой подключён":"частичное покрытие"}</span>
      </div>
      <div className={styles.annotationBox}>
        <div className="eyebrow">ФУНКЦИОНАЛЬНЫЕ КЛАСТЕРЫ MCL</div>
        {annotations.mcl_domains?.length?<div className={styles.annotationList} style={{marginTop:8}}>{annotations.mcl_domains.map((x)=><span className={styles.annotation} key={x}><b>Кластер</b>{x}</span>)}{annotations.subdomains?.map((x)=><span className={styles.annotation} key={x}><b>Подфункция</b>{x}</span>)}</div>:<div className={styles.annotationEmpty}>Функциональная разметка пока неполная.</div>}
        <div className={styles.annotationColumns}>
          <div><div className="eyebrow">КЛАСС БЕЛКА</div>{annotations.protein_classes?.length?<div className={styles.annotationList}>{annotations.protein_classes.map(x=><span className={styles.annotation} key={x}>{x}</span>)}</div>:<div className={styles.annotationEmpty}>—</div>}</div>
          <div><div className="eyebrow">КОМПАРТМЕНТ</div>{annotations.compartments?.length?<div className={styles.annotationList}>{annotations.compartments.map(x=><span className={styles.annotation} key={x}>{x}</span>)}</div>:<div className={styles.annotationEmpty}>—</div>}</div>
          <div><div className="eyebrow">ПРИЗНАКИ ОПУХОЛЕВОЙ БИОЛОГИИ</div>{annotations.hallmarks?.length?<div className={styles.annotationList}>{annotations.hallmarks.map(x=><span className={styles.annotation} key={x}>{x}</span>)}</div>:<div className={styles.annotationEmpty}>—</div>}</div>
        </div>
        <details className={styles.details}><summary>Показать происхождение аннотаций и формальные термины</summary><div className={styles.annotationList}>{[...domainDetails,...classDetails,...compartmentDetails,...hallmarkDetails].map((item:any,index:number)=><span className={styles.annotation} key={`${item.annotation_type}-${item.source_id}-${index}`}><b>{item.source||"источник"}</b>{item.source_term_name||item.annotation_label_ru}</span>)}</div>{formal.length>0&&<div className={styles.annotationList} style={{marginTop:10}}>{formal.map((item:any,index:number)=><span className={styles.annotation} key={`${item.source}-${item.term_id}-${index}`}><b>{item.source}</b>{item.term_name||item.term_id}</span>)}</div>}</details>
      </div>
    </section>

    <section className={styles.section}>
      <div className={styles.sectionHeader}>
        <div><div className="eyebrow">МУТАЦИОННЫЙ ФОН</div><h2>Какие генетические изменения могут быть связаны с зависимостью?</h2><p>Этот анализ сравнивает Gene Effect между моделями с функциональными мутациями и без них. Он предназначен для генерации гипотез, а не для доказательства причинности.</p></div>
      </div>
      {mutationAssociations?.available&&mutationAssociations.items?.length?<div className={styles.tableWrap}><table className={styles.table}><thead><tr><th>Мутированный ген</th><th>mut n</th><th>WT n</th><th>GE mut</th><th>GE WT</th><th>Δ GE</th><th>FDR</th></tr></thead><tbody>{mutationAssociations.items.slice(0,20).map((row:any)=><tr key={row.mutation_gene}><td><Link className={styles.geneLink} href={`/genes/${encodeURIComponent(row.mutation_gene)}`}>{row.mutation_gene}</Link></td><td>{row.mutated_n}</td><td>{row.wildtype_n}</td><td>{formatNumber(row.mutated_median_gene_effect)}</td><td>{formatNumber(row.wildtype_median_gene_effect)}</td><td className={Number(row.delta_median_gene_effect)<0?styles.negative:""}>{formatNumber(row.delta_median_gene_effect)}</td><td>{formatNumber(row.q_value,4)}</td></tr>)}</tbody></table></div>:<div className={styles.empty}>{mutationAssociations?.deferred?"Мутационный скрининг вынесен из первичной загрузки, чтобы карточка открывалась быстро. Его подключим как расчёт по запросу пользователя.":mutationAssociations?.note||"Подходящие мутационные группы пока не найдены."}</div>}
    </section>

    <section className={styles.section}>
      <details className={styles.technicalDetails}>
        <summary><div><div className="eyebrow">СТАТИСТИКА И КОНТРОЛЬ КАЧЕСТВА</div><h2>Показать подробные сравнительные анализы</h2><p>Δ Gene Effect, FDR, Cliff’s δ, размеры групп и технические предупреждения скрыты здесь, потому что они нужны для проверки вывода, а не для первичного чтения карточки.</p></div><span className={styles.badge}>{contexts.total} сравнений</span></summary>
        <div className={styles.contextGrid}>{contexts.items.map((row:any)=>{const flags=[row.broad_dependency_warning?"широкая зависимость":null,row.low_sample_size?"малая выборка":null].filter(Boolean);return <article className={styles.contextCard} key={row.comparison_id}><div className={styles.contextCardTop}><span>{row.cancer_name||row.cancer_id}</span><Link className={styles.geneLink} href={`/comparisons/${encodeURIComponent(row.comparison_id)}`}>открыть →</Link></div><h3>{row.comparison_label||row.comparison_id}</h3><div className={styles.contextMetrics}><div><span>Δ Gene Effect</span><b className={Number(row.delta_gene_effect)<0?styles.negative:""}>{formatNumber(row.delta_gene_effect)}</b></div><div><span>Target / comparator</span><b>{formatNumber(row.context_median_gene_effect)} / {formatNumber(row.comparator_median_gene_effect)}</b></div><div><span>FDR</span><b>{formatNumber(row.q_value,4)}</b></div><div><span>Cliff’s δ</span><b>{formatNumber(row.cliffs_delta)}</b></div><div><span>Target n</span><b>{formatNumber(row.context_models_n,0)}</b></div><div><span>Comparator n</span><b>{formatNumber(row.comparator_models_n,0)}</b></div></div><div className={styles.flags}>{flags.length?flags.map(x=><span className={styles.flag} key={x}>{x}</span>):<span className={styles.clean}>без отмеченных флагов</span>}</div></article>})}</div>
      </details>
    </section>

    <section className={styles.guardrails}>
      <div><b>CRISPR-зависимость ≠ лекарственная мишень</b>Потеря гена после CRISPR-выключения не гарантирует аналогичный эффект фармакологического ингибирования.</div>
      <div><b>Корреляция ≠ причинность</b>Связь с РНК, числом копий или мутацией помогает формировать гипотезы, но не устанавливает механизм.</div>
      <div><b>Клеточные модели ≠ пациенты</b>DepMap описывает экспериментальные модели; клиническую значимость нужно подтверждать отдельно.</div>
    </section>
  </>;
}
