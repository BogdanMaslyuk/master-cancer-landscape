import Link from "next/link";
import styles from "./ModelDependencyPanel.module.css";

type DependencyItem = {
  rank:number;
  rank_percentile?:number|null;
  gene:string;
  gene_name?:string|null;
  gene_effect:number;
  dependency_type:string;
  dependency_type_ru:string;
  pan_dependency_fraction?:number|null;
  cancer_dependency_fraction?:number|null;
  cancer_models_n?:number|null;
  domains?:{id:string;label_ru:string}[];
};

type FacetItem = {id:string;label_ru:string;genes_n:number};
type Payload = {
  available?:boolean;
  status?:string;
  note_ru?:string;
  build_command?:string;
  genes_measured_n?:number;
  strong_dependencies_n?:number;
  filtered_total?:number;
  offset?:number;
  limit?:number;
  items?:DependencyItem[];
  type_counts?:Record<string,number>;
  facets?:{dependency_types?:FacetItem[];domains?:FacetItem[]};
  interpretation?:Record<string,string>;
};

type Filters = {q?:string;type?:string;domain?:string;offset?:number};

function pct(value:number|null|undefined){
  return value===null||value===undefined ? "—" : `${Math.round(value*100)}%`;
}

function effectWidth(value:number){
  const strength=Math.max(0,Math.min(1.5,-value));
  return `${Math.round((strength/1.5)*100)}%`;
}

function pageHref(modelId:string, filters:Filters, offset:number){
  const q=new URLSearchParams();
  if(filters.q)q.set("dep_q",filters.q);
  if(filters.type)q.set("dep_type",filters.type);
  if(filters.domain)q.set("dep_domain",filters.domain);
  if(offset>0)q.set("dep_offset",String(offset));
  return `/models/${encodeURIComponent(modelId)}${q.toString()?`?${q.toString()}`:""}`;
}

export default function ModelDependencyPanel({payload,modelId,filters}:{payload:Payload;modelId:string;filters:Filters}){
  if(!payload.available){
    return <div className={styles.build}>
      <b>Полный CRISPR-профиль этой модели ещё не проиндексирован</b>
      <p>{payload.note_ru || "Нужно перестроить локальный индекс Gene Effect для полного CRISPR-атласа."}</p>
      {payload.build_command && <code>{payload.build_command}</code>}
    </div>;
  }

  const items=payload.items||[];
  const typeCounts=payload.type_counts||{};
  const total=payload.genes_measured_n||0;
  const filtered=payload.filtered_total||0;
  const offset=payload.offset||0;
  const limit=payload.limit||200;
  const hasPrev=offset>0;
  const hasNext=offset+items.length<filtered;

  return <div className={styles.panel}>
    <div className={styles.summary}>
      <div><strong>{total.toLocaleString("ru-RU")}</strong><span>генов с измеренным CRISPR Gene Effect</span></div>
      <div><strong>{(payload.strong_dependencies_n||0).toLocaleString("ru-RU")}</strong><span>сильных зависимостей при Gene Effect ≤ −0.5</span></div>
      <div><strong>{(typeCounts.broad_core||0).toLocaleString("ru-RU")}</strong><span>широких базовых зависимостей в опухолевых моделях DepMap</span></div>
      <div><strong>{((typeCounts.selective||0)+(typeCounts.model_selective||0)+(typeCounts.cancer_enriched||0)).toLocaleString("ru-RU")}</strong><span>селективных или обогащённых зависимостей для дальнейшего изучения</span></div>
    </div>

    <div className={styles.guardrail}>
      <span>i</span>
      <div><b>Что означает «специфичность» здесь</b><p>{payload.interpretation?.classification_ru || "Классификация основана на распространённости зависимости среди опухолевых моделей DepMap и не является оценкой безопасности для нормальных тканей."}</p></div>
    </div>

    <form className={styles.toolbar} method="get">
      <label className={styles.field}><span>Ген</span><input name="dep_q" defaultValue={filters.q||""} placeholder="KRAS, EGFR, PARP1…"/></label>
      <label className={styles.field}><span>Тип зависимости</span><select name="dep_type" defaultValue={filters.type||""}><option value="">Все типы</option>{(payload.facets?.dependency_types||[]).map(x=><option value={x.id} key={x.id}>{x.label_ru} · {x.genes_n}</option>)}</select></label>
      <label className={styles.field}><span>Функциональный кластер</span><select name="dep_domain" defaultValue={filters.domain||""}><option value="">Все кластеры</option>{(payload.facets?.domains||[]).map(x=><option value={x.id} key={x.id}>{x.label_ru} · {x.genes_n}</option>)}</select></label>
      <button className={styles.apply} type="submit">Применить</button>
    </form>

    {items.length ? <>
      <div className={styles.tableWrap}>
        <table className={styles.table}>
          <thead><tr><th>Ранг</th><th>Ген</th><th>Влияние на жизнеспособность</th><th>Тип зависимости</th><th>Функциональный кластер</th><th>Насколько распространена</th></tr></thead>
          <tbody>{items.map(item=><tr key={item.gene}>
            <td className={styles.rank}>#{item.rank}</td>
            <td><div className={styles.gene}><Link href={`/genes/${encodeURIComponent(item.gene)}`}>{item.gene}</Link>{item.gene_name&&<small>{item.gene_name}</small>}</div></td>
            <td><div className={styles.effect}><strong>{item.gene_effect.toFixed(3)}</strong><div className={styles.track}><span className={styles.fill} style={{width:effectWidth(item.gene_effect)}}/></div><small>{item.rank_percentile!==null&&item.rank_percentile!==undefined?`сильнее ${item.rank_percentile.toFixed(1)}% измеренных генов этой модели`:"—"}</small></div></td>
            <td><span className={`${styles.tag} ${styles[item.dependency_type as keyof typeof styles]||""}`}>{item.dependency_type_ru}</span></td>
            <td><div className={styles.clusters}>{(item.domains||[]).length?(item.domains||[]).slice(0,3).map(d=><span className={styles.cluster} key={d.id}>{d.label_ru}</span>):<span className={styles.muted}>не классифицирован</span>}</div></td>
            <td><div className={styles.prevalence}><b>{pct(item.pan_dependency_fraction)} моделей DepMap</b><span>имеют Gene Effect ≤ −0.5{item.cancer_dependency_fraction!==null&&item.cancer_dependency_fraction!==undefined&&item.cancer_models_n?`; в этой опухоли — ${pct(item.cancer_dependency_fraction)} (n=${item.cancer_models_n})`:""}</span></div></td>
          </tr>)}</tbody>
        </table>
      </div>

      <div className={styles.pagination}>
        <span>Показано {offset+1}–{offset+items.length} из {filtered.toLocaleString("ru-RU")} генов после фильтров.</span>
        <div className={styles.pageLinks}>{hasPrev&&<Link href={pageHref(modelId,filters,Math.max(0,offset-limit))}>← Предыдущие</Link>}{hasNext&&<Link href={pageHref(modelId,filters,offset+limit)}>Следующие →</Link>}</div>
      </div>
    </> : <div className={styles.empty}>По выбранным фильтрам гены не найдены.</div>}

    <div className={styles.guardrail}>
      <span>!</span>
      <div><b>CRISPR-нокаут ≠ лекарственное ингибирование</b><p>{payload.interpretation?.pharmacology_ru || "Даже сильная зависимость требует отдельной проверки лекарственной доступности мишени, селективности и токсичности."}</p></div>
    </div>
  </div>;
}
