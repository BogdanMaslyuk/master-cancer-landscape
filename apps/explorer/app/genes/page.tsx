import Link from "next/link";
import GeneSearch from "../../components/GeneSearch";
import { apiGet, formatNumber } from "../../lib/api";
import styles from "./GeneExplorer.module.css";

export const dynamic = "force-dynamic";

type SearchPayload = {
  page:number;
  page_size:number;
  total:number;
  pages:number;
  sort_by:string;
  sort_order:string;
  items:Record<string,any>[];
};

type Comparison = { id:string; label:string; cancer_id:string };
type AtlasContext = { id:string; cancer_ru?:string; molecular_ru?:string; short_ru?:string };
type AtlasPayload = { contexts?:AtlasContext[] } | AtlasContext[];
type Params = Record<string, string | string[] | undefined>;

function one(value:string|string[]|undefined){ return Array.isArray(value) ? value[0] : value; }
function on(value:string|undefined){ return value === "true" || value === "1" || value === "on"; }

function apiQuery(searchParams:Params){
  const out=new URLSearchParams();
  const names=["q","cancer_id","comparison_id","gene_effect_max","delta_gene_effect_max","q_value_max","cliffs_delta_abs_min","page","page_size","sort_by","sort_order"];
  for(const name of names){ const value=one(searchParams[name]); if(value) out.set(name,value); }
  for(const name of ["stable_only","exclude_broad","exclude_low_sample"]){ if(on(one(searchParams[name]))) out.set(name,"true"); }
  if(!out.has("page_size")) out.set("page_size","50");
  return out;
}

function pageHref(searchParams:Params,page:number){
  const q=apiQuery(searchParams);
  q.set("page",String(page));
  return `/genes?${q.toString()}`;
}

function activeFilterLabels(searchParams:Params,contexts:AtlasContext[],comparisons:Comparison[]){
  const labels:string[]=[];
  const q=one(searchParams.q); if(q) labels.push(`Поиск: ${q}`);
  const cancer=one(searchParams.cancer_id); if(cancer){ const found=contexts.find(x=>x.id===cancer); labels.push(found?.short_ru || found?.cancer_ru || cancer); }
  const cmp=one(searchParams.comparison_id); if(cmp){ const found=comparisons.find(x=>x.id===cmp); labels.push(found?.label || cmp); }
  const ge=one(searchParams.gene_effect_max); if(ge) labels.push(`min Gene Effect ≤ ${ge}`);
  const delta=one(searchParams.delta_gene_effect_max); if(delta) labels.push(`ΔGE ≤ ${delta}`);
  const qv=one(searchParams.q_value_max); if(qv) labels.push(`FDR ≤ ${qv}`);
  const cliff=one(searchParams.cliffs_delta_abs_min); if(cliff) labels.push(`|Cliff's δ| ≥ ${cliff}`);
  if(on(one(searchParams.stable_only))) labels.push("устойчивые 50/100/200");
  if(on(one(searchParams.exclude_broad))) labels.push("без broad dependency");
  if(on(one(searchParams.exclude_low_sample))) labels.push("без малых выборок");
  return labels;
}

export default async function GenesPage({searchParams}:{searchParams:Promise<Params>}){
  const sp=await searchParams;
  const query=apiQuery(sp);
  const [result, comparisons, atlasRaw]=await Promise.all([
    apiGet<SearchPayload>(`/api/genes/search?${query.toString()}`),
    apiGet<Comparison[]>("/api/comparisons"),
    apiGet<AtlasPayload>("/api/atlas"),
  ]);
  const contexts=Array.isArray(atlasRaw) ? atlasRaw : (atlasRaw.contexts || []);
  const active=activeFilterLabels(sp,contexts,comparisons);
  const currentCancer=one(sp.cancer_id) || "";
  const availableComparisons=currentCancer ? comparisons.filter(x=>x.cancer_id===currentCancer) : comparisons;

  return <>
    <section className={styles.pageHeader}>
      <div>
        <div className="eyebrow">UNIVERSAL GENE EXPLORER</div>
        <h1>Гены и мишени</h1>
        <p>Начните с любого гена и проследите обратный маршрут: в каких опухолевых контекстах возникает CRISPR-сигнал, какие клеточные модели сильнее всего зависят от гена и какие multi-omics данные поддерживают интерпретацию.</p>
      </div>
      <div className={styles.searchSlot}>
        <GeneSearch />
        <div className={styles.searchHint}>Можно открыть AHR, KRAS, DNM1L и любой другой gene symbol из текущего gene universe — не только устойчивые Top-N кандидаты.</div>
      </div>
    </section>

    <div className={styles.layout}>
      <aside className={styles.filters}>
        <div className={styles.filterHead}><strong>Фильтры</strong><Link className={styles.reset} href="/genes">Сбросить всё</Link></div>
        <form action="/genes" method="get">
          <div className={styles.group}>
            <div className={styles.groupTitle}>Поиск</div>
            <label className={styles.field}><span>Gene symbol содержит</span><input name="q" defaultValue={one(sp.q) || ""} placeholder="например AHR" /></label>
          </div>

          <div className={styles.group}>
            <div className={styles.groupTitle}>Опухолевый контекст</div>
            <label className={styles.field}><span>Контекст MCL</span><select name="cancer_id" defaultValue={currentCancer}><option value="">Все контексты</option>{contexts.map((ctx)=><option key={ctx.id} value={ctx.id}>{ctx.short_ru || `${ctx.cancer_ru || ctx.id} · ${ctx.molecular_ru || ""}`}</option>)}</select></label>
            <label className={styles.field}><span>Функциональное сравнение</span><select name="comparison_id" defaultValue={one(sp.comparison_id) || ""}><option value="">Все сравнения</option>{availableComparisons.map((cmp)=><option key={cmp.id} value={cmp.id}>{cmp.label}</option>)}</select></label>
          </div>

          <div className={styles.group}>
            <div className={styles.groupTitle}>CRISPR и статистика</div>
            <label className={styles.field}><span>Сильнейший model-level Gene Effect ≤</span><input name="gene_effect_max" type="number" step="0.05" defaultValue={one(sp.gene_effect_max) || ""} placeholder="-0.5" /></label>
            <label className={styles.field}><span>Δ Gene Effect ≤</span><input name="delta_gene_effect_max" type="number" step="0.05" defaultValue={one(sp.delta_gene_effect_max) || ""} placeholder="-0.2" /></label>
            <label className={styles.field}><span>q-value / FDR ≤</span><input name="q_value_max" type="number" min="0" max="1" step="0.01" defaultValue={one(sp.q_value_max) || ""} placeholder="0.05" /></label>
            <label className={styles.field}><span>|Cliff's δ| ≥</span><input name="cliffs_delta_abs_min" type="number" min="0" max="1" step="0.05" defaultValue={one(sp.cliffs_delta_abs_min) || ""} placeholder="0.3" /></label>
          </div>

          <div className={styles.group}>
            <div className={styles.groupTitle}>Надёжность</div>
            <label className={styles.check}><input type="checkbox" name="stable_only" value="true" defaultChecked={on(one(sp.stable_only))}/><span>Только устойчивые при Top-50/100/200</span></label>
            <label className={styles.check}><input type="checkbox" name="exclude_broad" value="true" defaultChecked={on(one(sp.exclude_broad))}/><span>Исключить broad dependency</span></label>
            <label className={styles.check}><input type="checkbox" name="exclude_low_sample" value="true" defaultChecked={on(one(sp.exclude_low_sample))}/><span>Исключить low sample</span></label>
          </div>

          <div className={styles.group}>
            <div className={styles.groupTitle}>Сортировка</div>
            <label className={styles.field}><span>Показатель</span><select name="sort_by" defaultValue={one(sp.sort_by) || "best_delta_gene_effect"}><option value="best_delta_gene_effect">Контекстная селективность ΔGE</option><option value="best_model_gene_effect">Сильнейший Gene Effect модели</option><option value="best_q_value">Минимальный q-value</option><option value="best_cliffs_delta">Cliff's δ</option><option value="significant_comparisons_n">Число значимых сравнений</option><option value="gene_symbol">Алфавит</option></select></label>
            <label className={styles.field}><span>Порядок</span><select name="sort_order" defaultValue={one(sp.sort_order) || "asc"}><option value="asc">По возрастанию</option><option value="desc">По убыванию</option></select></label>
            <input type="hidden" name="page_size" value={one(sp.page_size) || "50"}/>
          </div>
          <button className={styles.apply} type="submit">Применить фильтры</button>
        </form>
      </aside>

      <section className={styles.results}>
        <div className={styles.summaryBar}>
          <div><strong>{new Intl.NumberFormat("ru-RU").format(result.total)} генов</strong><span>Показаны данные, которые реально присутствуют в текущих MCL/DepMap индексах.</span></div>
          {active.length>0 && <div className={styles.activeFilters}>{active.map((label)=><span className={styles.chip} key={label}>{label}</span>)}</div>}
        </div>

        <div className={styles.tableWrap}>
          {result.items.length ? <table className={styles.table}>
            <thead><tr><th>Ген</th><th>Лучший контекст</th><th>Δ Gene Effect</th><th>min model GE</th><th>min q-value</th><th>Cliff's δ</th><th>Сравнения</th><th>Ограничения</th></tr></thead>
            <tbody>{result.items.map((gene:any)=>{
              const stable=String(gene.present_all_thresholds).toLowerCase()==="true";
              const delta=Number(gene.best_delta_gene_effect);
              const warnings=[gene.broad_dependency_any?"broad":null,gene.low_sample_any?"low n":null].filter(Boolean);
              return <tr key={gene.gene_symbol}>
                <td><Link className={styles.geneLink} href={`/genes/${encodeURIComponent(gene.gene_symbol)}`}>{stable && <span className={styles.stableDot}/>} {gene.gene_symbol}</Link></td>
                <td className={styles.contextCell}>{gene.best_context_label || "—"}</td>
                <td className={Number.isFinite(delta)&&delta<0?styles.negative:styles.metricStrong}>{formatNumber(gene.best_delta_gene_effect)}</td>
                <td className={styles.metricStrong}>{formatNumber(gene.best_model_gene_effect)}</td>
                <td>{formatNumber(gene.best_q_value,4)}</td>
                <td>{formatNumber(gene.best_cliffs_delta)}</td>
                <td>{gene.significant_comparisons_n ?? 0} знач. / {gene.comparisons_n ?? 0}</td>
                <td>{warnings.length ? <span className={styles.warning}>{warnings.join(" · ")}</span> : <span className={styles.clean}>нет флагов</span>}</td>
              </tr>;
            })}</tbody>
          </table> : <div className={styles.empty}><strong>По выбранным условиям гены не найдены</strong><span>Ослабьте один из порогов или сбросьте фильтры. Пустой результат не означает отсутствие биологической роли.</span></div>}
        </div>

        <div className={styles.pagination}>
          <span>Страница {result.page} из {result.pages}</span>
          <div className={styles.pageButtons}>
            <Link className={`${styles.pageButton} ${result.page<=1?styles.disabled:""}`} href={pageHref(sp,Math.max(1,result.page-1))}>← Назад</Link>
            <Link className={`${styles.pageButton} ${result.page>=result.pages?styles.disabled:""}`} href={pageHref(sp,Math.min(result.pages,result.page+1))}>Дальше →</Link>
          </div>
        </div>

        <div className={styles.guardrail}>
          <div><b>CRISPR ≠ лекарство</b>Сильная зависимость после knockout не доказывает, что фармакологическое ингибирование воспроизведёт эффект.</div>
          <div><b>Клеточные линии ≠ пациенты</b>Контекстная селективность относится к экспериментальным моделям и не является частотой признака в опухолях пациентов.</div>
          <div><b>Нет данных ≠ нет эффекта</b>RNA, CN или статистика могут отсутствовать из-за неполного покрытия конкретного слоя.</div>
        </div>
      </section>
    </div>
  </>;
}
