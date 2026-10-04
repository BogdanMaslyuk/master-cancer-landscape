import Link from "next/link";
import { apiGet, formatNumber } from "../../../lib/api";
import styles from "./GeneCancerMatrix.module.css";

export const dynamic = "force-dynamic";

type Params = Record<string,string|string[]|undefined>;
type Facet = {id:string;label_ru:string;genes_n:number};
type Domain = Facet & {subdomains?:Facet[]};
type Facets = {domains:Domain[];protein_classes?:Facet[];compartments?:Facet[];hallmarks?:Facet[]};
type AtlasContext = {id:string;short_ru?:string;cancer_ru?:string;molecular_ru?:string};
type AtlasPayload = {contexts?:AtlasContext[];organs?:{contexts?:AtlasContext[]}[]} | AtlasContext[];
type Comparison = {id:string;label:string;cancer_id:string;context_models_n?:number;comparator_models_n?:number;qc_status?:string};
type Cell = {
  delta_gene_effect?:number|null;
  context_median_gene_effect?:number|null;
  comparator_median_gene_effect?:number|null;
  q_value?:number|null;
  cliffs_delta?:number|null;
  fdr_0_05?:boolean;
  broad_dependency_warning?:boolean;
  low_sample_size?:boolean;
};
type MatrixRow = {gene_symbol:string;gene_name?:string|null;present_all_thresholds?:boolean;cells:Record<string,Cell|null>};
type MatrixPayload = {
  genes_total_after_filters:number;
  genes_returned_n:number;
  comparisons_n:number;
  comparisons:Comparison[];
  rows:MatrixRow[];
};

function one(value:string|string[]|undefined){return Array.isArray(value)?value[0]:value;}
function on(value:string|undefined,defaultValue=false){if(value===undefined)return defaultValue;return value==="true"||value==="1"||value==="on";}
function contexts(value:AtlasPayload):AtlasContext[]{if(Array.isArray(value))return value;if(value.contexts?.length)return value.contexts;return (value.organs||[]).flatMap(x=>x.contexts||[]);}
function cellClass(cell:Cell|null|undefined){
  if(!cell || cell.delta_gene_effect===null || cell.delta_gene_effect===undefined) return `${styles.cell} ${styles.noData}`;
  const delta=Number(cell.delta_gene_effect);
  const magnitude=delta<=-0.5?styles.strongNegative:delta<=-0.2?styles.mediumNegative:delta<0?styles.weakNegative:styles.positive;
  const significant=cell.fdr_0_05?styles.significant:"";
  const warning=cell.broad_dependency_warning||cell.low_sample_size?styles.warning:"";
  return `${styles.cell} ${magnitude} ${significant} ${warning}`;
}

function queryFrom(sp:Params){
  const q=new URLSearchParams();
  for(const key of ["q","domain","subdomain","protein_class","compartment","hallmark","cancer_id","gene_effect_max","delta_gene_effect_max","q_value_max","cliffs_delta_abs_min","limit","sort_by","sort_order"]){const value=one(sp[key]);if(value)q.set(key,value);}
  if(on(one(sp.stable_only)))q.set("stable_only","true");
  if(on(one(sp.exclude_broad),true))q.set("exclude_broad","true"); else q.set("exclude_broad","false");
  if(on(one(sp.exclude_low_sample),true))q.set("exclude_low_sample","true"); else q.set("exclude_low_sample","false");
  if(!q.has("limit"))q.set("limit","60");
  return q;
}

export default async function MatrixPage({searchParams}:{searchParams:Promise<Params>}){
  const sp=await searchParams;
  const query=queryFrom(sp);
  const [matrix,facets,atlasRaw]=await Promise.all([
    apiGet<MatrixPayload>(`/api/gene-matrix?${query.toString()}`),
    apiGet<Facets>("/api/genes/facets"),
    apiGet<AtlasPayload>("/api/atlas"),
  ]);
  const atlas=contexts(atlasRaw);
  const currentDomain=one(sp.domain)||"";
  const selectedDomain=facets.domains.find(x=>x.id===currentDomain);
  const subdomains=selectedDomain?.subdomains?.length?selectedDomain.subdomains:facets.domains.flatMap(x=>x.subdomains||[]);

  return <>
    <section className={styles.header}>
      <div><div className="eyebrow">GENE × CANCER CONTEXT</div><h1>Матрица функциональных зависимостей</h1><p>Строка — ген, столбец — выполненное MCL-сравнение. Основное число в клетке — Δ Gene Effect = median(context) − median(comparator). Цвет показывает направление и величину эффекта; статистическая значимость отображается отдельно точкой.</p></div>
      <Link className={styles.back} href="/genes">← Вернуться к каталогу генов</Link>
    </section>

    <div className={styles.layout}>
      <aside className={styles.filters}>
        <form action="/genes/matrix" method="get">
          <div className={styles.group}><div className={styles.groupTitle}>Гены</div><label className={styles.field}><span>Символ / название</span><input name="q" defaultValue={one(sp.q)||""} placeholder="AHR, KRAS…"/></label></div>
          <div className={styles.group}><div className={styles.groupTitle}>Биология</div>
            <label className={styles.field}><span>Functional Domain</span><select name="domain" defaultValue={currentDomain}><option value="">Все домены</option>{facets.domains.map(x=><option key={x.id} value={x.id}>{x.label_ru} · {x.genes_n}</option>)}</select></label>
            <label className={styles.field}><span>Подфункция</span><select name="subdomain" defaultValue={one(sp.subdomain)||""}><option value="">Все</option>{subdomains.map(x=><option key={x.id} value={x.id}>{x.label_ru} · {x.genes_n}</option>)}</select></label>
            <label className={styles.field}><span>Класс белка</span><select name="protein_class" defaultValue={one(sp.protein_class)||""}><option value="">Все</option>{(facets.protein_classes||[]).map(x=><option key={x.id} value={x.id}>{x.label_ru} · {x.genes_n}</option>)}</select></label>
            <label className={styles.field}><span>Компартмент</span><select name="compartment" defaultValue={one(sp.compartment)||""}><option value="">Все</option>{(facets.compartments||[]).map(x=><option key={x.id} value={x.id}>{x.label_ru} · {x.genes_n}</option>)}</select></label>
            <label className={styles.field}><span>Hallmark</span><select name="hallmark" defaultValue={one(sp.hallmark)||""}><option value="">Все</option>{(facets.hallmarks||[]).map(x=><option key={x.id} value={x.id}>{x.label_ru} · {x.genes_n}</option>)}</select></label>
          </div>
          <div className={styles.group}><div className={styles.groupTitle}>Контекст</div><label className={styles.field}><span>Опухоль / молекулярный контекст</span><select name="cancer_id" defaultValue={one(sp.cancer_id)||""}><option value="">Все доступные сравнения</option>{atlas.map(x=><option key={x.id} value={x.id}>{x.short_ru||`${x.cancer_ru||x.id} · ${x.molecular_ru||""}`}</option>)}</select></label></div>
          <div className={styles.group}><div className={styles.groupTitle}>CRISPR / статистика</div>
            <label className={styles.field}><span>min model Gene Effect ≤</span><input name="gene_effect_max" type="number" step="0.05" defaultValue={one(sp.gene_effect_max)||""} placeholder="-0.5"/></label>
            <label className={styles.field}><span>Δ Gene Effect ≤</span><input name="delta_gene_effect_max" type="number" step="0.05" defaultValue={one(sp.delta_gene_effect_max)||""} placeholder="-0.2"/></label>
            <label className={styles.field}><span>q-value ≤</span><input name="q_value_max" type="number" min="0" max="1" step="0.01" defaultValue={one(sp.q_value_max)||""} placeholder="0.05"/></label>
            <label className={styles.field}><span>|Cliff's δ| ≥</span><input name="cliffs_delta_abs_min" type="number" min="0" max="1" step="0.05" defaultValue={one(sp.cliffs_delta_abs_min)||""} placeholder="0.3"/></label>
          </div>
          <div className={styles.group}><div className={styles.groupTitle}>Надёжность</div>
            <label className={styles.check}><input name="stable_only" type="checkbox" value="true" defaultChecked={on(one(sp.stable_only))}/><span>Только устойчивые Top-50/100/200</span></label>
            <label className={styles.check}><input name="exclude_broad" type="checkbox" value="true" defaultChecked={on(one(sp.exclude_broad),true)}/><span>Исключить broad dependency</span></label>
            <label className={styles.check}><input name="exclude_low_sample" type="checkbox" value="true" defaultChecked={on(one(sp.exclude_low_sample),true)}/><span>Исключить low sample</span></label>
          </div>
          <div className={styles.group}><div className={styles.groupTitle}>Отображение</div><label className={styles.field}><span>Генов</span><select name="limit" defaultValue={one(sp.limit)||"60"}><option value="30">30</option><option value="60">60</option><option value="100">100</option><option value="120">120</option></select></label></div>
          <button className={styles.apply} type="submit">Обновить матрицу</button>
        </form>
      </aside>

      <section className={styles.content}>
        <div className={styles.summary}>
          <div><strong>{matrix.genes_returned_n} из {matrix.genes_total_after_filters} генов</strong><span>{matrix.comparisons_n} выполненных сравнений. • рядом с ΔGE означает FDR &lt; 0.05; нижняя янтарная линия — broad/low-sample предупреждение.</span></div>
          <div className={styles.legend}><span>Δ ≤ −0.5</span><span>−0.5…−0.2</span><span>−0.2…0</span><span>Δ ≥ 0</span></div>
        </div>

        {matrix.rows.length && matrix.comparisons.length ? <div className={styles.matrixWrap}><table className={styles.matrix}>
          <thead><tr><th className={styles.geneHead}>Ген</th>{matrix.comparisons.map(c=><th className={styles.comparisonHead} key={c.id}><b>{c.label}</b><small>target n={c.context_models_n??"—"} · comparator n={c.comparator_models_n??"—"}<br/>QC: {c.qc_status||"—"}</small></th>)}</tr></thead>
          <tbody>{matrix.rows.map(row=><tr key={row.gene_symbol}><td className={styles.geneCell}><Link className={styles.geneLink} href={`/genes/${encodeURIComponent(row.gene_symbol)}`}>{row.gene_symbol}</Link>{row.gene_name&&<span className={styles.geneName}>{row.gene_name}</span>}</td>{matrix.comparisons.map(c=>{const cell=row.cells[c.id];return <td className={cellClass(cell)} key={c.id}>{cell?<><span className={styles.cellValue}>{formatNumber(cell.delta_gene_effect)}</span><span className={styles.cellMeta}>q {formatNumber(cell.q_value,3)} · δ {formatNumber(cell.cliffs_delta)}</span></>:"—"}</td>;})}</tr>)}</tbody>
        </table></div> : <div className={styles.empty}>Для выбранных условий матрица пуста. Ослабьте биологический или статистический фильтр.</div>}

        <div className={styles.noteGrid}>
          <div><b>Цвет ≠ значимость</b>Фон кодирует только направление и величину Δ Gene Effect. Статистическая поддержка показана отдельно.</div>
          <div><b>Межбелковое сравнение осторожно</b>Матрица нужна для навигации и поиска паттернов; абсолютные значения не являются универсальным рейтингом «лучшей мишени».</div>
          <div><b>Покрытие ограничено MCL</b>Столбцы появляются только для реально выполненных genome-wide сравнений. Отсутствующий опухолевый контекст нельзя интерпретировать как отрицательный результат.</div>
        </div>
      </section>
    </div>
  </>;
}
