import Link from "next/link";
import { apiGet, formatNumber } from "../../../lib/api";
import GeneInsightCharts from "./GeneInsightCharts";
import styles from "./GeneWorkbench.module.css";

export const dynamic = "force-dynamic";

type DescriptivePayload={available:boolean;items:Record<string,any>[];guardrail?:string};
type CorrelationPayload={available:boolean;relationships:Record<string,any>[];points:Record<string,any>[];statistics_engine?:string};
type MutationPayload={available:boolean;profiled_models_n?:number;items:Record<string,any>[];statistics_engine?:string;note?:string};
type BasePayload = {
  identity:Record<string,any>;
  summary:Record<string,any>;
  stability:Record<string,any>;
  comparisons?:Record<string,any>[];
  pathways?:Record<string,any>[];
  model_insights?:{
    descriptive_contexts?:DescriptivePayload;
    correlations?:CorrelationPayload;
    mutation_associations?:MutationPayload;
  };
};
type ContextPayload = { gene_symbol:string; total:number; items:Record<string,any>[] };
type ModelsPayload = {
  gene_symbol:string;
  available:boolean;
  available_layers:Record<string,boolean>;
  total:number;
  page:number;
  page_size:number;
  pages?:number;
  items:Record<string,any>[];
  guardrails?:Record<string,string>;
  note?:string;
};
type AnnotationPayload = {
  gene_symbol:string;
  status:string;
  mcl_domains:string[];
  subdomains:string[];
  protein_classes:string[];
  compartments:string[];
  hallmarks:string[];
  mcl_domain_details?:Record<string,any>[];
  protein_class_details?:Record<string,any>[];
  compartment_details?:Record<string,any>[];
  hallmark_details?:Record<string,any>[];
  formal_annotations:Record<string,any>[];
  coverage?:{annotated_genes_n:number;gene_universe_n:number;coverage_fraction:number|null;status:string;reason:string};
  reference_coverage?:{resolved_genes_n:number;gene_universe_n:number;coverage_fraction:number|null;available:boolean};
  reference_available?:boolean;
  taxonomy_version?:string;
  note?:string;
};
type AtlasContext = { id:string; short_ru?:string; cancer_ru?:string; molecular_ru?:string };
type AtlasPayload = { contexts?:AtlasContext[]; organs?:{contexts?:AtlasContext[]}[] } | AtlasContext[];
type SearchParams = Record<string,string|string[]|undefined>;

function one(value:string|string[]|undefined){ return Array.isArray(value)?value[0]:value; }
function yes(value:unknown){ return String(value).toLowerCase()==="true"; }
function atlasContexts(value:AtlasPayload):AtlasContext[]{
  if(Array.isArray(value)) return value;
  if(value.contexts?.length) return value.contexts;
  return (value.organs || []).flatMap((organ)=>organ.contexts || []);
}
function sourceDetails(items:Record<string,any>[]|undefined){ return (items||[]).slice(0,24); }
function pct(value:unknown){const n=Number(value);return Number.isFinite(n)?`${Math.round(n*100)}%`:"—";}

export default async function GenePage({params,searchParams}:{params:Promise<{gene:string}>;searchParams:Promise<SearchParams>}){
  const [{gene},sp]=await Promise.all([params,searchParams]);
  const symbol=decodeURIComponent(gene).trim().toUpperCase();
  const modelQuery=new URLSearchParams({page_size:"100",sort_by:"gene_effect",sort_order:"asc"});
  const cancerId=one(sp.cancer_id); if(cancerId) modelQuery.set("cancer_id",cancerId);
  const geneEffectMax=one(sp.gene_effect_max); if(geneEffectMax) modelQuery.set("gene_effect_max",geneEffectMax);

  const [base,contexts,models,annotations,atlasRaw]=await Promise.all([
    apiGet<BasePayload>(`/api/genes/${encodeURIComponent(symbol)}`),
    apiGet<ContextPayload>(`/api/genes/${encodeURIComponent(symbol)}/contexts`),
    apiGet<ModelsPayload>(`/api/genes/${encodeURIComponent(symbol)}/models?${modelQuery.toString()}`),
    apiGet<AnnotationPayload>(`/api/genes/${encodeURIComponent(symbol)}/annotations`),
    apiGet<AtlasPayload>("/api/atlas"),
  ]);
  const atlas=atlasContexts(atlasRaw);
  const stable=yes(base.stability?.present_all_thresholds);
  const significant=Number(base.summary?.significant_comparisons_n||0);
  const comparisonN=Number(base.summary?.comparisons_n||0);
  const bestDelta=base.summary?.best_delta_gene_effect;
  const bestModel=base.summary?.best_model_gene_effect;
  const bestQ=base.summary?.best_q_value;
  const formal=(annotations.formal_annotations||[]).slice(0,50);
  const domainDetails=sourceDetails(annotations.mcl_domain_details);
  const classDetails=sourceDetails(annotations.protein_class_details);
  const compartmentDetails=sourceDetails(annotations.compartment_details);
  const hallmarkDetails=sourceDetails(annotations.hallmark_details);
  const aliases=(base.identity?.aliases || []) as string[];
  const descriptive=base.model_insights?.descriptive_contexts;
  const correlations=base.model_insights?.correlations;
  const mutationAssociations=base.model_insights?.mutation_associations;

  return <>
    <div className="breadcrumbs"><Link href="/genes">Гены и мишени</Link><span>›</span><strong>{symbol}</strong></div>

    <section className={styles.hero}>
      <div>
        <div className="eyebrow" style={{color:"#b9d5ee"}}>GENE WORKBENCH · TARGET → CANCER</div>
        <h1>{symbol}</h1>
        {base.identity?.gene_name&&<div style={{fontSize:18,fontWeight:700,color:"#e5eff9",marginTop:-4,marginBottom:9}}>{base.identity.gene_name}</div>}
        <p>Карточка связывает фундаментальную биологию гена с описательной зависимостью клеточных линий, валидированными опухолевыми сравнениями и multi-omics контекстом.</p>
        {aliases.length>0&&<p style={{marginTop:8,fontSize:12.5}}>Синонимы: {aliases.slice(0,8).join(", ")}{aliases.length>8?" …":""}</p>}
      </div>
      <div className={styles.heroMeta}>
        <div><span>HGNC</span><b>{base.identity?.hgnc_id || "—"}</b></div>
        <div><span>Ensembl</span><b>{base.identity?.ensembl_gene_id || "—"}</b></div>
        <div><span>UniProt</span><b>{base.identity?.uniprot_id || "не подключено"}</b></div>
        <div><span>Статус MCL</span><b>{stable?"устойчивый 50/100/200":"произвольный / контекстный ген"}</b></div>
      </div>
    </section>

    <section className={styles.metrics}>
      <div className={styles.metric}><span>Лучший Δ Gene Effect</span><strong>{formatNumber(bestDelta)}</strong></div>
      <div className={styles.metric}><span>Сильнейший model-level Gene Effect</span><strong>{formatNumber(bestModel)}</strong></div>
      <div className={styles.metric}><span>Минимальный q-value</span><strong>{formatNumber(bestQ,4)}</strong></div>
      <div className={styles.metric}><span>Значимые сравнения</span><strong>{significant}/{comparisonN}</strong></div>
      <div className={styles.metric}><span>Клеточные модели с доступными слоями</span><strong>{models.total}</strong></div>
    </section>

    <section className={styles.section}>
      <div className={styles.sectionHeader}>
        <div><div className="eyebrow">ФУНКЦИОНАЛЬНАЯ ПРИНАДЛЕЖНОСТЬ</div><h2>Что это за ген с точки зрения клеточной биологии?</h2><p>Категории multi-label и выводятся из формальных source terms; они помогают навигации, но не заменяют первичную онтологическую аннотацию.</p></div>
        <span className={styles.badge}>{annotations.reference_available?"reference snapshot подключён":"частичное покрытие"}</span>
      </div>
      <div className={styles.annotationBox}>
        <div className="eyebrow">MCL FUNCTIONAL DOMAINS</div>
        {annotations.mcl_domains?.length ? <>
          <div className={styles.annotationList} style={{marginTop:8}}>{annotations.mcl_domains.map((domain)=><span className={styles.annotation} key={domain}><b>Домен</b>{domain}</span>)}</div>
          {annotations.subdomains?.length>0 && <div className={styles.annotationList} style={{marginTop:10}}>{annotations.subdomains.map((sub)=><span className={styles.annotation} key={sub}><b>Подфункция</b>{sub}</span>)}</div>}
        </> : <div className={styles.annotationEmpty}>Для {symbol} доменная разметка пока отсутствует. Это означает неполное покрытие, а не отсутствие фундаментальной функции.</div>}
        <div style={{marginTop:18}}><div className="eyebrow">МОЛЕКУЛЯРНЫЙ КЛАСС БЕЛКА</div>{annotations.protein_classes?.length?<div className={styles.annotationList} style={{marginTop:8}}>{annotations.protein_classes.map((item)=><span className={styles.annotation} key={item}><b>Класс</b>{item}</span>)}</div>:<div className={styles.annotationEmpty}>Класс пока не определён.</div>}</div>
        <div style={{marginTop:18}}><div className="eyebrow">КЛЕТОЧНЫЙ КОМПАРТМЕНТ</div>{annotations.compartments?.length?<div className={styles.annotationList} style={{marginTop:8}}>{annotations.compartments.map((item)=><span className={styles.annotation} key={item}><b>Локализация</b>{item}</span>)}</div>:<div className={styles.annotationEmpty}>Компартмент пока не определён из подключённых GO:CC-аннотаций.</div>}</div>
        <div style={{marginTop:18}}><div className="eyebrow">HALLMARKS OF CANCER</div>{annotations.hallmarks?.length?<div className={styles.annotationList} style={{marginTop:8}}>{annotations.hallmarks.map((item)=><span className={styles.annotation} key={item}><b>Hallmark</b>{item}</span>)}</div>:<div className={styles.annotationEmpty}>Hallmark-проекция пока не сформирована.</div>}</div>
        {(domainDetails.length+classDetails.length+compartmentDetails.length+hallmarkDetails.length)>0&&<div style={{marginTop:20}}><div className="eyebrow">ПРОИСХОЖДЕНИЕ КАТЕГОРИЙ</div><div className={styles.annotationList} style={{marginTop:8}}>{[...domainDetails,...classDetails,...compartmentDetails,...hallmarkDetails].slice(0,36).map((item:any,index:number)=><span className={styles.annotation} key={`${item.annotation_type}-${item.annotation_id}-${item.source_id}-${index}`}><b>{item.source || "source"} · {item.source_id || "—"}</b>{item.source_term_name || item.annotation_label_ru}</span>)}</div></div>}
        <div style={{marginTop:20}}><div className="eyebrow">ФОРМАЛЬНЫЕ ТЕРМИНЫ</div>{formal.length ? <div className={styles.annotationList} style={{marginTop:8}}>{formal.map((item:any,index:number)=><span className={styles.annotation} key={`${item.source}-${item.term_id}-${index}`}><b>{item.source} · {item.term_id || "—"}</b>{item.term_name || item.source_term_name || item.term_id}</span>)}</div> : <div className={styles.annotationEmpty}>Формальные GO/pathway-термины для {symbol} в локальном snapshot пока отсутствуют.</div>}</div>
        {annotations.coverage&&<p className={styles.annotationEmpty} style={{marginBottom:0}}>Доменная разметка покрывает {annotations.coverage.annotated_genes_n} из {annotations.coverage.gene_universe_n} генов текущего universe.</p>}
      </div>
    </section>

    <section className={styles.section}>
      <div className={styles.sectionHeader}>
        <div><div className="eyebrow">ОПИСАТЕЛЬНО · БЕЗ ТРЕБОВАНИЯ COMPARATOR</div><h2>Как выглядит зависимость от {symbol} в разных модельных контекстах?</h2><p>Этот блок показывает наблюдаемый Gene Effect в целевых MCL-группах даже там, где корректного comparator нет. Он нужен для поиска сигналов, но не доказывает контекстную селективность.</p></div>
        <span className={styles.badge}>descriptive layer</span>
      </div>
      {descriptive?.available&&descriptive.items?.length?<div className={styles.tableWrap}><table className={styles.table}>
        <thead><tr><th>Контекст</th><th>Target median GE</th><th>IQR</th><th>Target n</th><th>GE &lt; −0.5</th><th>GE &lt; −1</th><th>Comparator</th></tr></thead>
        <tbody>{descriptive.items.map((row:any)=><tr key={row.cancer_id}>
          <td><Link className={styles.geneLink} href={`/atlas/${encodeURIComponent(row.cancer_id)}`}>{row.cancer_name||row.cancer_id}</Link><span className={styles.modelName}>{row.context_definition||row.cancer_id}</span></td>
          <td className={Number(row.context?.median)<0?styles.negative:""}>{formatNumber(row.context?.median)}</td>
          <td>{formatNumber(row.context?.q1)} … {formatNumber(row.context?.q3)}</td>
          <td>{row.context?.n??0}</td>
          <td>{pct(row.context?.fraction_lt_minus_0_5)}</td>
          <td>{pct(row.context?.fraction_lt_minus_1)}</td>
          <td>{row.has_comparator?`${formatNumber(row.comparator?.median)} · n=${row.comparator?.n??0}`:"нет корректной группы"}</td>
        </tr>)}</tbody>
      </table></div>:<div className={styles.empty}>Описательный model-level Gene Effect для {symbol} пока недоступен.</div>}
    </section>

    <section className={styles.section}>
      <div className={styles.sectionHeader}>
        <div><div className="eyebrow">СТАТИСТИЧЕСКОЕ СРАВНЕНИЕ</div><h2>Где зависимость контекстно селективна?</h2><p>Здесь остаются только заранее определённые target vs comparator сравнения с Δ Gene Effect, Cliff's δ и FDR.</p></div>
        <span className={styles.badge}>{contexts.total} сравнений</span>
      </div>
      {contexts.items.length ? <div className={styles.contextGrid}>{contexts.items.map((row:any)=>{
        const delta=Number(row.delta_gene_effect);
        const flags=[row.broad_dependency_warning?"broad dependency":null,row.low_sample_size?"low sample":null].filter(Boolean);
        return <article className={styles.contextCard} key={row.comparison_id}>
          <div className={styles.contextCardTop}><span>{row.cancer_name || row.cancer_id}</span><Link className={styles.geneLink} href={`/comparisons/${encodeURIComponent(row.comparison_id)}`}>открыть сравнение →</Link></div>
          <h3>{row.comparison_label || row.comparison_id}</h3>
          <div className={styles.contextMetrics}>
            <div><span>Δ Gene Effect</span><b className={Number.isFinite(delta)&&delta<0?styles.negative:""}>{formatNumber(row.delta_gene_effect)}</b></div>
            <div><span>Target / comparator</span><b>{formatNumber(row.context_median_gene_effect)} / {formatNumber(row.comparator_median_gene_effect)}</b></div>
            <div><span>q-value</span><b>{formatNumber(row.q_value,4)}</b></div>
            <div><span>Cliff's δ</span><b>{formatNumber(row.cliffs_delta)}</b></div>
            <div><span>Target n</span><b>{formatNumber(row.context_models_n,0)}</b></div>
            <div><span>Comparator n</span><b>{formatNumber(row.comparator_models_n,0)}</b></div>
          </div>
          <div className={styles.flags}>{flags.length?flags.map((flag)=><span className={styles.flag} key={flag}>{flag}</span>):<span className={styles.clean}>без отмеченных технических флагов</span>}</div>
        </article>;
      })}</div> : <div className={styles.empty}>Для {symbol} нет рассчитанных target-vs-comparator сравнений. Это не отменяет описательные сигналы выше.</div>}
    </section>

    <section className={styles.section}>
      <div className={styles.sectionHeader}>
        <div><div className="eyebrow">GENE → CELL MODELS</div><h2>Какие клеточные линии сильнее всего зависят от {symbol}?</h2><p>Модели сортируются по Gene Effect по возрастанию.</p></div>
        <div className={styles.layerState}><span className={models.available_layers?.gene_effect?styles.on:""}>CRISPR {models.available_layers?.gene_effect?"доступен":"нет"}</span><span className={models.available_layers?.expression?styles.on:""}>RNA {models.available_layers?.expression?"доступна":"нет"}</span><span className={models.available_layers?.copy_number?styles.on:""}>CN {models.available_layers?.copy_number?"доступен":"нет"}</span></div>
      </div>
      <form action={`/genes/${encodeURIComponent(symbol)}`} method="get" className={styles.modelControls}>
        <label><span>Опухолевый контекст</span><select name="cancer_id" defaultValue={cancerId||""}><option value="">Все модели MCL</option>{atlas.map((ctx)=><option key={ctx.id} value={ctx.id}>{ctx.short_ru || `${ctx.cancer_ru||ctx.id} · ${ctx.molecular_ru||""}`}</option>)}</select></label>
        <label><span>Gene Effect ≤</span><input type="number" step="0.05" name="gene_effect_max" defaultValue={geneEffectMax||""} placeholder="например -0.5"/></label>
        <button type="submit">Фильтровать модели</button>
        {(cancerId||geneEffectMax)&&<Link className={styles.geneLink} href={`/genes/${encodeURIComponent(symbol)}`}>Сбросить</Link>}
      </form>
      {models.available && models.items.length ? <div className={styles.tableWrap}><table className={styles.table}>
        <thead><tr><th>Модель</th><th>Происхождение</th><th>Контексты MCL</th><th>Gene Effect</th><th>RNA · log2(TPM+1)</th><th>Relative CN</th></tr></thead>
        <tbody>{models.items.map((row:any)=>{const memberships=(row.memberships||[]) as any[];return <tr key={row.model_id}><td><Link className={styles.modelLink} href={`/models/${encodeURIComponent(row.model_id)}`}>{row.cell_line_name || row.model_id}</Link>{row.cell_line_name&&<span className={styles.modelName}>{row.model_id}</span>}</td><td>{row.oncotree_subtype || row.oncotree_primary_disease || row.oncotree_lineage || "—"}</td><td>{memberships.length?memberships.map((m:any)=>`${m.cancer_id} · ${m.assigned_group}`).join("; "):"—"}</td><td className={Number(row.gene_effect)<0?styles.negative:""}>{formatNumber(row.gene_effect)}</td><td>{formatNumber(row.expression)}</td><td>{formatNumber(row.copy_number)}</td></tr>;})}</tbody>
      </table></div> : <div className={styles.empty}>{models.note || `Для ${symbol} model-level multi-omics данные пока не подключены.`}</div>}
    </section>

    <section className={styles.section}>
      <div className={styles.sectionHeader}>
        <div><div className="eyebrow">MULTI-OMICS RELATIONSHIPS</div><h2>Связана ли зависимость с RNA или copy number?</h2><p>Spearman ρ рассчитывается по клеточным моделям с обоими доступными измерениями. Это исследовательская корреляция, а не причинный вывод.</p></div>
        <span className={styles.badge}>{correlations?.statistics_engine||"нет данных"}</span>
      </div>
      {correlations?.available?<GeneInsightCharts points={correlations.points||[]} relationships={correlations.relationships||[]}/>:<div className={styles.empty}>Недостаточно multi-omics данных для расчёта корреляций.</div>}
    </section>

    <section className={styles.section}>
      <div className={styles.sectionHeader}>
        <div><div className="eyebrow">МУТАЦИОННЫЙ ФОН</div><h2>Какие мутации ассоциированы с зависимостью от {symbol}?</h2><p>Для каждого часто встречающегося функционального варианта сравнивается Gene Effect в mutation-positive и mutation-negative моделях. Более отрицательный Δ median GE означает более сильную зависимость в мутантной группе.</p></div>
        <span className={styles.badge}>{mutationAssociations?.profiled_models_n??0} профилированных моделей</span>
      </div>
      {mutationAssociations?.available&&mutationAssociations.items?.length?<div className={styles.tableWrap}><table className={styles.table}>
        <thead><tr><th>Мутированный ген</th><th>mut n</th><th>WT n</th><th>median GE mut</th><th>median GE WT</th><th>Δ median GE</th><th>Cliff's δ</th><th>p</th><th>q/FDR</th></tr></thead>
        <tbody>{mutationAssociations.items.slice(0,30).map((row:any)=><tr key={row.mutation_gene}><td><Link className={styles.geneLink} href={`/genes/${encodeURIComponent(row.mutation_gene)}`}>{row.mutation_gene}</Link></td><td>{row.mutated_n}</td><td>{row.wildtype_n}</td><td>{formatNumber(row.mutated_median_gene_effect)}</td><td>{formatNumber(row.wildtype_median_gene_effect)}</td><td className={Number(row.delta_median_gene_effect)<0?styles.negative:""}>{formatNumber(row.delta_median_gene_effect)}</td><td>{formatNumber(row.cliffs_delta)}</td><td>{formatNumber(row.p_value,4)}</td><td>{formatNumber(row.q_value,4)}</td></tr>)}</tbody>
      </table></div>:<div className={styles.empty}>{mutationAssociations?.note||"Статистически/технически пригодных мутационных групп пока не найдено."}</div>}
    </section>

    <section className={styles.guardrails}>
      <div><b>Descriptive ≠ selective</b>Низкий Gene Effect в группе моделей — повод исследовать сигнал, но контекстную селективность подтверждает только корректное target-vs-comparator сравнение.</div>
      <div><b>Correlation ≠ causation</b>Связь Gene Effect с RNA, CN или мутацией не доказывает причинный механизм.</div>
      <div><b>CRISPR knockout ≠ ингибитор</b>Gene Effect описывает генетическое выключение и не гарантирует эффект малой молекулы.</div>
      <div><b>Модели ≠ пациенты</b>DepMap характеризует клеточные линии и не подменяет пациентские когорты или клиническую частоту.</div>
    </section>
  </>;
}
