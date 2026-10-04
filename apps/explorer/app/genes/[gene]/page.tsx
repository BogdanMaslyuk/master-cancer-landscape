import Link from "next/link";
import { apiGet, formatNumber } from "../../../lib/api";
import styles from "./GeneWorkbench.module.css";

export const dynamic = "force-dynamic";

type BasePayload = {
  identity:Record<string,any>;
  summary:Record<string,any>;
  stability:Record<string,any>;
  comparisons?:Record<string,any>[];
  pathways?:Record<string,any>[];
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
  formal_annotations:Record<string,any>[];
  coverage?:{annotated_genes_n:number;gene_universe_n:number;coverage_fraction:number|null;status:string;reason:string};
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
  const formal=(annotations.formal_annotations||[]).slice(0,40);
  const domainDetails=(annotations.mcl_domain_details||[]).slice(0,20);

  return <>
    <div className="breadcrumbs"><Link href="/genes">Гены и мишени</Link><span>›</span><strong>{symbol}</strong></div>

    <section className={styles.hero}>
      <div>
        <div className="eyebrow" style={{color:"#b9d5ee"}}>GENE WORKBENCH · TARGET → CANCER</div>
        <h1>{symbol}</h1>
        <p>Карточка показывает обратный маршрут: функциональная принадлежность → опухолевые сравнения → отдельные клеточные модели → Gene Effect → RNA → copy number.</p>
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
        <div><div className="eyebrow">ФУНКЦИОНАЛЬНАЯ ПРИНАДЛЕЖНОСТЬ</div><h2>К каким биологическим процессам относится {symbol}?</h2><p>MCL Functional Domains — человекочитаемый multi-label слой над формальными терминами. Домен назначается только тогда, когда существующий GO/Reactome/KEGG/CORUM-термин удовлетворяет прозрачному правилу; исходный термин сохраняется как provenance.</p></div>
        <span className={styles.badge}>{annotations.mcl_domains?.length ? `${annotations.mcl_domains.length} доменов` : "нет доменной разметки"}</span>
      </div>
      <div className={styles.annotationBox}>
        {annotations.mcl_domains?.length ? <>
          <div className={styles.annotationList}>{annotations.mcl_domains.map((domain)=><span className={styles.annotation} key={domain}><b>MCL domain</b>{domain}</span>)}</div>
          {annotations.subdomains?.length>0 && <div className={styles.annotationList} style={{marginTop:10}}>{annotations.subdomains.map((sub)=><span className={styles.annotation} key={sub}><b>Подфункция</b>{sub}</span>)}</div>}
        </> : <div className={styles.annotationEmpty}>Для {symbol} доменная разметка пока отсутствует. Это означает неполное покрытие текущего аннотационного слоя, а не отсутствие фундаментальной функции.</div>}

        {domainDetails.length>0 && <div style={{marginTop:16}}>
          <div className="eyebrow">ПРОИСХОЖДЕНИЕ ДОМЕННОЙ РАЗМЕТКИ</div>
          <div className={styles.annotationList} style={{marginTop:8}}>{domainDetails.map((item:any,index:number)=><span className={styles.annotation} key={`${item.annotation_id}-${item.source_id}-${index}`}><b>{item.source || "source"} · {item.source_id || "—"}</b>{item.source_term_name || item.annotation_label_ru}</span>)}</div>
        </div>}

        <div style={{marginTop:18}}>
          <div className="eyebrow">ФОРМАЛЬНЫЕ ТЕРМИНЫ</div>
          {formal.length ? <div className={styles.annotationList} style={{marginTop:8}}>{formal.map((item:any,index:number)=><span className={styles.annotation} key={`${item.source}-${item.term_id}-${index}`}><b>{item.source} · {item.term_id || "—"}</b>{item.term_name || item.source_term_name || item.term_id}</span>)}</div> : <div className={styles.annotationEmpty}>В текущем MCL enrichment-layer для {symbol} формальных терминов не найдено. Полное gene-to-ontology покрытие ещё не загружено.</div>}
        </div>

        {annotations.coverage&&<p className={styles.annotationEmpty} style={{marginBottom:0}}>Текущее функциональное покрытие: {annotations.coverage.annotated_genes_n} из {annotations.coverage.gene_universe_n} генов. Непокрытые гены нельзя интерпретировать как «не относящиеся» к функциям.</p>}
        {annotations.note&&<p className={styles.annotationEmpty} style={{marginBottom:0}}>{annotations.note}</p>}
      </div>
    </section>

    <section className={styles.section}>
      <div className={styles.sectionHeader}>
        <div><div className="eyebrow">ГДЕ ГЕН ВАЖЕН В MCL</div><h2>Опухолевые контексты</h2><p>Сравнения отсортированы по Δ Gene Effect: более отрицательное значение означает более сильную зависимость целевой группы относительно comparator. Отдельно сохраняются статистика и предупреждения.</p></div>
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
      })}</div> : <div className={styles.empty}>Для {symbol} нет рассчитанных опухолевых сравнений в текущем наборе MCL. Ген всё равно может иметь model-level данные DepMap ниже.</div>}
    </section>

    <section className={styles.section}>
      <div className={styles.sectionHeader}>
        <div><div className="eyebrow">GENE → CELL MODELS</div><h2>Какие клеточные линии сильнее всего зависят от {symbol}?</h2><p>По умолчанию модели сортируются по Gene Effect по возрастанию. Это позволяет быстро увидеть линии, где CRISPR-выключение гена связано с наиболее сильной потерей приспособленности.</p></div>
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
        <tbody>{models.items.map((row:any)=>{
          const memberships=(row.memberships||[]) as any[];
          return <tr key={row.model_id}>
            <td><Link className={styles.modelLink} href={`/models/${encodeURIComponent(row.model_id)}`}>{row.cell_line_name || row.model_id}</Link>{row.cell_line_name&&<span className={styles.modelName}>{row.model_id}</span>}</td>
            <td>{row.oncotree_subtype || row.oncotree_primary_disease || row.oncotree_lineage || "—"}</td>
            <td>{memberships.length?memberships.map((m:any)=>`${m.cancer_id} · ${m.assigned_group}`).join("; "):"—"}</td>
            <td className={Number(row.gene_effect)<0?styles.negative:""}>{formatNumber(row.gene_effect)}</td>
            <td>{formatNumber(row.expression)}</td>
            <td>{formatNumber(row.copy_number)}</td>
          </tr>;
        })}</tbody>
      </table></div> : <div className={styles.empty}>{models.note || `Для ${symbol} model-level multi-omics данные пока не подключены.`} {models.available && !models.items.length ? "По текущим фильтрам моделей не найдено." : ""}</div>}
    </section>

    <section className={styles.guardrails}>
      <div><b>Домен ≠ новое биологическое доказательство</b>MCL Functional Domain — навигационная проекция формального термина; исходный GO/Reactome/KEGG/CORUM-термин остаётся первичным доказательным слоем.</div>
      <div><b>CRISPR knockout ≠ ингибитор</b>Gene Effect описывает генетическое выключение функции и не доказывает воспроизводимость эффекта малой молекулой.</div>
      <div><b>RNA ≠ активный белок</b>Экспрессия помогает интерпретации модели, но сама по себе не подтверждает уровень, локализацию или активность белка.</div>
      <div><b>Модели ≠ пациентская частота</b>Данные DepMap характеризуют экспериментальные клеточные линии и не подменяют реальные пациентские когорты.</div>
    </section>
  </>;
}
