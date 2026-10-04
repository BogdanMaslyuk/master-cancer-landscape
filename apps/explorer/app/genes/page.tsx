import Link from "next/link";
import GeneSearch from "../../components/GeneSearch";
import { apiGet } from "../../lib/api";
import styles from "./GeneExplorer.module.css";
import catalog from "./GeneCatalogV2.module.css";

export const dynamic = "force-dynamic";

type SearchPayload = {
  page:number;
  page_size:number;
  total:number;
  pages:number;
  sort_by:string;
  sort_order:string;
  functional_coverage?:{annotated_genes_n:number;gene_universe_n:number;coverage_fraction:number|null;status:string;reason:string};
  reference_coverage?:{resolved_genes_n:number;gene_universe_n:number;coverage_fraction:number|null;available:boolean};
  items:Record<string,any>[];
};
type Comparison = { id:string; label:string; cancer_id:string };
type AtlasContext = { id:string; cancer_ru?:string; molecular_ru?:string; short_ru?:string };
type AtlasPayload = { contexts?:AtlasContext[]; organs?:{contexts?:AtlasContext[]}[] } | AtlasContext[];
type Facet = { id:string; label_ru:string; genes_n:number };
type FacetDomain = { id:string; label_ru:string; label_en?:string; genes_n:number; subdomains:Facet[] };
type FacetPayload = {
  domains:FacetDomain[];
  protein_classes?:Facet[];
  compartments?:Facet[];
  hallmarks?:Facet[];
  sources:{id:string;genes_n:number}[];
  coverage:{annotated_genes_n:number;gene_universe_n:number;coverage_fraction:number|null;status:string;reason:string};
  reference_coverage?:{resolved_genes_n:number;gene_universe_n:number;coverage_fraction:number|null;available:boolean};
};
type Params = Record<string,string|string[]|undefined>;

const dependencyTypes = [
  ["broad_baseline","Широкая базовая"],
  ["broad_tumor","Широкая опухолевая"],
  ["cancer_enriched","Обогащена в типе опухоли"],
  ["selective","Селективная"],
  ["intermediate","Промежуточная"],
  ["weak","Слабая / невыраженная"],
] as const;

function one(value:string|string[]|undefined){return Array.isArray(value)?value[0]:value;}
function on(value:string|undefined){return value==="true"||value==="1"||value==="on";}
function atlasContexts(value:AtlasPayload):AtlasContext[]{if(Array.isArray(value))return value;if(value.contexts?.length)return value.contexts;return (value.organs||[]).flatMap(x=>x.contexts||[]);}
function jsonList(value:unknown):string[]{if(Array.isArray(value))return value.map(String);try{const parsed=JSON.parse(String(value||"[]"));return Array.isArray(parsed)?parsed.map(String):[];}catch{return [];}}
function facetLabel(items:Facet[]|undefined,id:string|undefined){return items?.find(x=>x.id===id)?.label_ru||id;}
function pct(value:unknown){const n=Number(value);return Number.isFinite(n)?`${Math.round(n*100)}%`:"—";}
function dependencyClass(type:string|undefined){return ({broad_baseline:catalog.broadBaseline,broad_tumor:catalog.broadTumor,cancer_enriched:catalog.cancerEnriched,selective:catalog.selective,intermediate:catalog.intermediate,weak:catalog.weak} as Record<string,string>)[String(type||"")]||catalog.intermediate;}

function apiQuery(sp:Params){
  const out=new URLSearchParams();
  const names=["q","domain","subdomain","pathway","annotation_source","protein_class","compartment","hallmark","cancer_id","comparison_id","gene_effect_max","delta_gene_effect_max","q_value_max","cliffs_delta_abs_min","dependency_type","dependency_fraction_min","specificity_score_min","page","page_size","sort_by","sort_order"];
  for(const name of names){const value=one(sp[name]);if(value)out.set(name,value);}
  for(const name of ["stable_only","exclude_broad","exclude_low_sample"]){if(on(one(sp[name])))out.set(name,"true");}
  if(!out.has("page_size"))out.set("page_size","50");
  return out;
}
function pageHref(sp:Params,page:number){const q=apiQuery(sp);q.set("page",String(page));return `/genes?${q.toString()}`;}

function activeFilterLabels(sp:Params,contexts:AtlasContext[],comparisons:Comparison[],facets:FacetPayload){
  const labels:string[]=[];
  const q=one(sp.q);if(q)labels.push(`Поиск: ${q}`);
  const domain=one(sp.domain);if(domain)labels.push(facets.domains.find(x=>x.id===domain)?.label_ru||domain);
  const sub=one(sp.subdomain);if(sub)labels.push(facets.domains.flatMap(x=>x.subdomains||[]).find(x=>x.id===sub)?.label_ru||sub);
  const pc=one(sp.protein_class);if(pc)labels.push(`Класс: ${facetLabel(facets.protein_classes,pc)}`);
  const compartment=one(sp.compartment);if(compartment)labels.push(`Компартмент: ${facetLabel(facets.compartments,compartment)}`);
  const hallmark=one(sp.hallmark);if(hallmark)labels.push(`Hallmark: ${facetLabel(facets.hallmarks,hallmark)}`);
  const dep=one(sp.dependency_type);if(dep)labels.push(dependencyTypes.find(x=>x[0]===dep)?.[1]||dep);
  const frac=one(sp.dependency_fraction_min);if(frac)labels.push(`Зависимых моделей ≥ ${Math.round(Number(frac)*100)}%`);
  const spec=one(sp.specificity_score_min);if(spec)labels.push(`Специфичность ≥ ${Math.round(Number(spec)*100)} п.п.`);
  const cancer=one(sp.cancer_id);if(cancer){const found=contexts.find(x=>x.id===cancer);labels.push(found?.short_ru||found?.cancer_ru||cancer);}
  const cmp=one(sp.comparison_id);if(cmp)labels.push(comparisons.find(x=>x.id===cmp)?.label||cmp);
  return labels;
}

export default async function GenesPage({searchParams}:{searchParams:Promise<Params>}){
  const sp=await searchParams;
  const query=apiQuery(sp);
  const [result,comparisons,atlasRaw,facets]=await Promise.all([
    apiGet<SearchPayload>(`/api/genes/search?${query.toString()}`),
    apiGet<Comparison[]>("/api/comparisons"),
    apiGet<AtlasPayload>("/api/atlas"),
    apiGet<FacetPayload>("/api/genes/facets"),
  ]);
  const contexts=atlasContexts(atlasRaw);
  const active=activeFilterLabels(sp,contexts,comparisons,facets);
  const currentCancer=one(sp.cancer_id)||"";
  const currentDomain=one(sp.domain)||"";
  const selectedDomain=facets.domains.find(x=>x.id===currentDomain);
  const availableSubdomains=selectedDomain?.subdomains?.length?selectedDomain.subdomains:facets.domains.flatMap(x=>x.subdomains||[]);
  const availableComparisons=currentCancer?comparisons.filter(x=>x.cancer_id===currentCancer):comparisons;
  const domainLabels=new Map(facets.domains.map(x=>[x.id,x.label_ru]));
  const coverage=facets.coverage;
  const referenceCoverage=facets.reference_coverage;

  return <>
    <section className={styles.pageHeader}>
      <div>
        <div className="eyebrow">ГЕН → ФУНКЦИЯ → ЗАВИСИМОСТЬ → ОПУХОЛЬ</div>
        <h1>Гены и мишени</h1>
        <p>Каталог показывает не статистику отдельного сравнения, а биологический профиль гена: что он делает, насколько часто опухолевые модели зависят от него и в каком типе опухоли эта зависимость наиболее выражена.</p>
      </div>
      <div className={styles.searchSlot}><GeneSearch/><div className={styles.searchHint}>Поиск работает по символу, названию и синонимам, если локальный reference snapshot построен.</div></div>
    </section>

    <div className={styles.layout}>
      <aside className={styles.filters}>
        <div className={styles.filterHead}><strong>Фильтры</strong><Link className={styles.reset} href="/genes">Сбросить всё</Link></div>
        <form action="/genes" method="get">
          <div className={styles.group}><div className={styles.groupTitle}>Поиск</div><label className={styles.field}><span>Ген / название / синоним</span><input name="q" defaultValue={one(sp.q)||""} placeholder="KRAS, EGFR, PARP1…"/></label></div>

          <div className={styles.group}><div className={styles.groupTitle}>Функция</div>
            <label className={styles.field}><span>Функциональный домен</span><select name="domain" defaultValue={currentDomain}><option value="">Все домены</option>{facets.domains.map(x=><option key={x.id} value={x.id}>{x.label_ru} · {x.genes_n}</option>)}</select></label>
            <label className={styles.field}><span>Подфункция</span><select name="subdomain" defaultValue={one(sp.subdomain)||""}><option value="">Все подфункции</option>{availableSubdomains.map(x=><option key={x.id} value={x.id}>{x.label_ru} · {x.genes_n}</option>)}</select></label>
            <label className={styles.field}><span>Класс белка</span><select name="protein_class" defaultValue={one(sp.protein_class)||""}><option value="">Все классы</option>{(facets.protein_classes||[]).map(x=><option key={x.id} value={x.id}>{x.label_ru} · {x.genes_n}</option>)}</select></label>
            <label className={styles.field}><span>Компартмент</span><select name="compartment" defaultValue={one(sp.compartment)||""}><option value="">Все компартменты</option>{(facets.compartments||[]).map(x=><option key={x.id} value={x.id}>{x.label_ru} · {x.genes_n}</option>)}</select></label>
            <label className={styles.field}><span>Hallmark of Cancer</span><select name="hallmark" defaultValue={one(sp.hallmark)||""}><option value="">Все Hallmarks</option>{(facets.hallmarks||[]).map(x=><option key={x.id} value={x.id}>{x.label_ru} · {x.genes_n}</option>)}</select></label>
          </div>

          <div className={styles.group}><div className={styles.groupTitle}>CRISPR-зависимость</div>
            <label className={styles.field}><span>Тип зависимости</span><select name="dependency_type" defaultValue={one(sp.dependency_type)||""}><option value="">Все типы</option>{dependencyTypes.map(([id,label])=><option key={id} value={id}>{label}</option>)}</select></label>
            <label className={styles.field}><span>Зависимых моделей не менее</span><select name="dependency_fraction_min" defaultValue={one(sp.dependency_fraction_min)||""}><option value="">Без порога</option><option value="0.05">5%</option><option value="0.10">10%</option><option value="0.25">25%</option><option value="0.50">50%</option><option value="0.75">75%</option></select></label>
            <label className={styles.field}><span>Специфичность не менее</span><select name="specificity_score_min" defaultValue={one(sp.specificity_score_min)||""}><option value="">Без порога</option><option value="0.15">средняя</option><option value="0.30">высокая</option><option value="0.50">очень высокая</option></select></label>
            <div className={catalog.typeHint}>Рабочий порог каталога: Gene Effect ≤ −0,5. Он используется только для описательной навигации по моделям DepMap.</div>
          </div>

          <details className={catalog.advanced}>
            <summary>Расширенный анализ</summary>
            <div className={catalog.advancedBody}>
              <label className={styles.field}><span>Контекст MCL</span><select name="cancer_id" defaultValue={currentCancer}><option value="">Все контексты</option>{contexts.map(x=><option key={x.id} value={x.id}>{x.short_ru||x.cancer_ru||x.id}</option>)}</select></label>
              <label className={styles.field}><span>Сравнение</span><select name="comparison_id" defaultValue={one(sp.comparison_id)||""}><option value="">Все сравнения</option>{availableComparisons.map(x=><option key={x.id} value={x.id}>{x.label}</option>)}</select></label>
              <label className={styles.field}><span>Сильнейший Gene Effect ≤</span><input name="gene_effect_max" type="number" step="0.05" defaultValue={one(sp.gene_effect_max)||""} placeholder="-0.5"/></label>
              <label className={styles.field}><span>Δ Gene Effect ≤</span><input name="delta_gene_effect_max" type="number" step="0.05" defaultValue={one(sp.delta_gene_effect_max)||""} placeholder="-0.2"/></label>
              <label className={styles.field}><span>FDR ≤</span><input name="q_value_max" type="number" min="0" max="1" step="0.01" defaultValue={one(sp.q_value_max)||""} placeholder="0.05"/></label>
              <label className={styles.field}><span>|Cliff's δ| ≥</span><input name="cliffs_delta_abs_min" type="number" min="0" max="1" step="0.05" defaultValue={one(sp.cliffs_delta_abs_min)||""} placeholder="0.3"/></label>
              <label className={styles.check}><input type="checkbox" name="stable_only" value="true" defaultChecked={on(one(sp.stable_only))}/><span>Только устойчивые Top-50/100/200</span></label>
              <label className={styles.check}><input type="checkbox" name="exclude_broad" value="true" defaultChecked={on(one(sp.exclude_broad))}/><span>Исключить broad dependency</span></label>
              <label className={styles.check}><input type="checkbox" name="exclude_low_sample" value="true" defaultChecked={on(one(sp.exclude_low_sample))}/><span>Исключить малые выборки</span></label>
            </div>
          </details>

          <div className={styles.group}><div className={styles.groupTitle}>Сортировка</div>
            <label className={styles.field}><span>Показатель</span><select name="sort_by" defaultValue={one(sp.sort_by)||"specificity_score"}><option value="specificity_score">Специфичность</option><option value="dependency_fraction">Распространённость зависимости</option><option value="best_cancer_dependency_fraction">Зависимость в лучшем типе опухоли</option><option value="global_median_gene_effect">Медианный Gene Effect</option><option value="functional_annotations_n">Функциональная изученность</option><option value="gene_symbol">Алфавит</option></select></label>
            <label className={styles.field}><span>Порядок</span><select name="sort_order" defaultValue={one(sp.sort_order)||"desc"}><option value="desc">Сначала больше / выше</option><option value="asc">Сначала меньше / ниже</option></select></label>
            <input type="hidden" name="page_size" value={one(sp.page_size)||"50"}/>
          </div>
          <button className={styles.apply} type="submit">Применить фильтры</button>
        </form>
      </aside>

      <section className={styles.results}>
        <div className={styles.summaryBar}>
          <div className={catalog.tableIntroText}><strong>{new Intl.NumberFormat("ru-RU").format(result.total)} генов</strong><span>Функционально размечено {new Intl.NumberFormat("ru-RU").format(coverage.annotated_genes_n)} из {new Intl.NumberFormat("ru-RU").format(coverage.gene_universe_n)}. {referenceCoverage?.available?`Справочник генов: ${new Intl.NumberFormat("ru-RU").format(referenceCoverage.resolved_genes_n)}.`:"Полный справочник генов ещё не построен."}</span></div>
          {active.length>0&&<div className={styles.activeFilters}>{active.map(label=><span className={styles.chip} key={label}>{label}</span>)}</div>}
        </div>

        <div className={catalog.methodNote}><b>Как читать:</b><span>«Зависимая модель» здесь означает Gene Effect ≤ −0,5. «Специфичность» — насколько доля зависимых моделей в наиболее чувствительном типе опухоли выше общей доли по доступным CRISPR-моделям. Это описательный показатель, а не доказательство лекарственной мишени.</span></div>

        <div className={styles.tableWrap}>
          {result.items.length?<table className={styles.table}>
            <thead><tr><th>Ген</th><th>Функция</th><th>Тип зависимости</th><th>Где особенно важен</th><th>Зависимые модели</th><th>Специфичность</th><th></th></tr></thead>
            <tbody>{result.items.map((gene:any)=>{
              const stable=String(gene.present_all_thresholds).toLowerCase()==="true";
              const domains=jsonList(gene.mcl_domains_json).map(id=>domainLabels.get(id)||id).slice(0,3);
              const fraction=Number(gene.dependency_fraction);
              const bestFraction=Number(gene.best_cancer_dependency_fraction);
              const score=Number(gene.specificity_score);
              return <tr key={gene.gene_symbol}>
                <td><Link className={styles.geneLink} href={`/genes/${encodeURIComponent(gene.gene_symbol)}`}>{stable&&<span className={styles.stableDot}/>} {gene.gene_symbol}</Link>{gene.gene_name&&<div style={{fontSize:10.5,color:"#7b8999",marginTop:3,maxWidth:190}}>{gene.gene_name}</div>}</td>
                <td className={catalog.functionCell}>{domains.length?domains.join(" · "):"Функция пока не размечена"}</td>
                <td><span className={`${catalog.dependencyTag} ${dependencyClass(gene.dependency_type)}`}>{gene.dependency_type_ru||"Не классифицирована"}</span></td>
                <td><div className={catalog.whereCell}><b>{gene.best_cancer_name||"—"}</b><span>{gene.best_cancer_organ_ru||""}{Number.isFinite(bestFraction)?` · зависимы ${pct(bestFraction)}`:""}</span></div></td>
                <td className={catalog.countCell}><div className={catalog.countTop}><b>{gene.dependency_models_n??0}</b><span>из {gene.dependency_models_total_n??"—"}</span></div><div className={catalog.bar}><i style={{width:Number.isFinite(fraction)?`${Math.max(0,Math.min(100,fraction*100))}%`:"0%"}}/></div><span className={catalog.percent}>{pct(fraction)} моделей</span></td>
                <td><div className={catalog.specificity}><b>{gene.specificity_label_ru||"—"}</b><span>{Number.isFinite(score)?`+${Math.round(score*100)} п.п. к общей доле`:"недостаточно данных"}</span></div></td>
                <td><Link className={catalog.openLink} href={`/genes/${encodeURIComponent(gene.gene_symbol)}`} aria-label={`Открыть ${gene.gene_symbol}`}>→</Link></td>
              </tr>;
            })}</tbody>
          </table>:<div className={styles.empty}><strong>По выбранным условиям гены не найдены</strong><span>Ослабьте фильтры. Пустой результат не означает отсутствие биологической роли.</span></div>}
        </div>

        <div className={styles.pagination}><span>Страница {result.page} из {result.pages}</span><div className={styles.pageButtons}><Link className={`${styles.pageButton} ${result.page<=1?styles.disabled:""}`} href={pageHref(sp,Math.max(1,result.page-1))}>← Назад</Link><Link className={`${styles.pageButton} ${result.page>=result.pages?styles.disabled:""}`} href={pageHref(sp,Math.min(result.pages,result.page+1))}>Дальше →</Link></div></div>

        <div className={styles.guardrail}>
          <div><b>Широкая зависимость ≠ хорошая мишень</b>Если ген нужен большинству опухолевых моделей, его выключение может отражать базовую клеточную функцию и потенциальный риск токсичности.</div>
          <div><b>Селективность — повод исследовать</b>Высокая специфичность помогает найти уязвимость, но требует проверки генетического контекста, воспроизводимости и нормальных тканей.</div>
          <div><b>CRISPR ≠ лекарство</b>Нокаут гена не равен частичному фармакологическому ингибированию белка.</div>
          <div><b>Модели ≠ пациенты</b>Все показатели на странице относятся к экспериментальным моделям DepMap, а не к частоте признака у пациентов.</div>
        </div>
      </section>
    </div>
  </>;
}
