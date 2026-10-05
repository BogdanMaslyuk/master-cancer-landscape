import Link from "next/link";
import { apiGet } from "../../lib/api";
import styles from "../pharmacology.module.css";

export const dynamic = "force-dynamic";

type Params=Record<string,string|string[]|undefined>;
type Payload=Record<string,any>;
function one(value:string|string[]|undefined){return Array.isArray(value)?value[0]:value;}
function n(value:any){const x=Number(value);return Number.isFinite(x)?x.toLocaleString("ru-RU"):"—";}
function pct(value:any){if(value===null||value===undefined||value==="")return "—";const x=Number(value);return Number.isFinite(x)?`${Math.round(x*100)}%`:"—";}
const statusRu:Record<string,string>={
  priority_for_in_vitro:"Приоритет для in vitro",
  supported_hypothesis:"Поддержанная",
  exploratory_hypothesis:"Исследовательская",
  insufficient_evidence:"Недостаточно данных",
};
const axisRu:Record<string,string>={
  strong:"сильная",supportive:"поддерживает",weak:"слабая",sparse:"мало данных",
  conflicting:"противоречие",uncertain:"неопределённо",not_assessable:"не сопоставимо",not_available:"нет данных",
  concordant_enrichment:"согласованное обогащение",partial_enrichment:"частичное обогащение",
  weak_enrichment:"слабое обогащение",not_enriched:"не обогащено",limited_sample:"малая группа",
  not_assessed:"не оценено",
};
function pageHref(sp:Params,offset:number){const q=new URLSearchParams();for(const k of ["q","status","target_gene","cancer","limit"]){const v=one(sp[k]);if(v)q.set(k,v);}q.set("offset",String(Math.max(0,offset)));return `/hypotheses?${q.toString()}`;}

export default async function HypothesesPage({searchParams}:{searchParams:Promise<Params>}){
  const sp=await searchParams;
  const limit=Math.max(1,Math.min(500,Number(one(sp.limit)||"100")||100));
  const offset=Math.max(0,Number(one(sp.offset)||"0")||0);
  const query=new URLSearchParams({limit:String(limit),offset:String(offset)});
  for(const k of ["q","status","target_gene","cancer"]){const v=one(sp[k]);if(v)query.set(k,v);}
  const [catalog,summary]=await Promise.all([
    apiGet<Payload>(`/api/hypotheses?${query.toString()}`),
    apiGet<Payload>("/api/hypotheses/summary"),
  ]);
  const items=(catalog.items||[]) as any[];
  const counts=summary.status_counts||{};
  const end=Math.min(catalog.total||0,offset+items.length);
  return <>
    <section className={styles.pageHeader}>
      <div>
        <div className="eyebrow">ВЕЩЕСТВО × БЕЛКОВАЯ МИШЕНЬ × ОПУХОЛЬ × КЛЕТОЧНАЯ МОДЕЛЬ</div>
        <h1>Исследовательские гипотезы</h1>
        <p>Центральный слой выбора следующего эксперимента. MCL не выдаёт искусственный единый балл: фенотип, CRISPR-зависимость, согласованность механизма, контекст и полнота данных показаны отдельно.</p>
      </div>
      <div className={styles.heroStats}>
        <div className={styles.heroStat}><strong>{n(counts.priority_for_in_vitro)}</strong><span>приоритет для in vitro</span></div>
        <div className={styles.heroStat}><strong>{n(counts.supported_hypothesis)}</strong><span>поддержанных гипотез</span></div>
        <div className={styles.heroStat}><strong>{n(summary.hypotheses_n)}</strong><span>всего сочетаний</span></div>
      </div>
    </section>

    <form action="/hypotheses" method="get" className={styles.toolbar}>
      <label className={styles.field}><span>Поиск</span><input name="q" defaultValue={one(sp.q)||""} placeholder="vemurafenib, BRAF, melanoma…"/></label>
      <label className={styles.field}><span>Статус</span><select name="status" defaultValue={one(sp.status)||""}><option value="">все</option><option value="priority_for_in_vitro">приоритет для in vitro</option><option value="supported_hypothesis">поддержанная</option><option value="exploratory_hypothesis">исследовательская</option><option value="insufficient_evidence">недостаточно данных</option></select></label>
      <label className={styles.field}><span>Мишень / ген</span><input name="target_gene" defaultValue={one(sp.target_gene)||""} placeholder="BRAF"/></label>
      <label className={styles.field}><span>Опухоль / орган</span><input name="cancer" defaultValue={one(sp.cancer)||""} placeholder="melanoma"/></label>
      <input type="hidden" name="limit" value={String(limit)}/>
      <button className={styles.button} type="submit">Применить</button>
      <Link className={styles.reset} href="/hypotheses">Сбросить</Link>
    </form>

    {!catalog.available?<div className={styles.empty}><b>Слой гипотез ещё не построен.</b><br/><code>{catalog.build_command||summary.build_command||"Запустите build_candidate_hypotheses.py"}</code></div>:
    <>
      <div className={styles.note}>{catalog.interpretation_ru||summary.contract_note_ru}</div>
      <div className={styles.tableWrap}><table className={styles.table}>
        <thead><tr><th>Кандидат</th><th>Белковая мишень</th><th>Опухолевый контекст</th><th>Фенотип</th><th>CRISPR</th><th>Механизм</th><th>Контекст</th><th>Статус</th></tr></thead>
        <tbody>{items.map((row:any)=><tr key={row.hypothesis_id}>
          <td><Link className={styles.primary} href={`/hypotheses/${encodeURIComponent(row.hypothesis_id)}`}>{row.preferred_name||row.compound_id}</Link><span className={styles.sub}>{row.compound_id}</span><span className={styles.sub}>{n(row.joint_support_models_n)} совместно поддерживающих моделей</span></td>
          <td><Link className={styles.primary} href={`/targets/${encodeURIComponent(row.target_gene)}`}>{row.protein_preferred_name||row.target_gene}</Link><span className={styles.sub}>ген <Link href={`/genes/${encodeURIComponent(row.target_gene)}`}>{row.target_gene}</Link>{row.uniprot_primary_accession?` · UniProt ${row.uniprot_primary_accession}`:""}</span></td>
          <td><b>{row.mcl_cancer_name||"—"}</b><span className={styles.sub}>{row.mcl_organ_ru||""} · {n(row.models_n)} моделей</span></td>
          <td><b>{axisRu[row.phenotype_axis]||row.phenotype_axis||"—"}</b><span className={styles.sub}>{pct(row.sensitive_fraction)} чувствительных</span></td>
          <td><b>{axisRu[row.dependency_axis]||row.dependency_axis||"—"}</b><span className={styles.sub}>{pct(row.dependency_fraction_in_cancer)} с GE ≤ −0,5</span></td>
          <td><b>{axisRu[row.mechanism_axis]||row.mechanism_axis||"—"}</b><span className={styles.sub}>{row.concordance_label||"профиль не оценён"}</span></td>
          <td><b>{axisRu[row.specificity_axis]||row.specificity_axis||"—"}</b><span className={styles.sub}>молекулярная генетика: ещё не включена</span></td>
          <td><Link className={styles.chip} href={`/hypotheses/${encodeURIComponent(row.hypothesis_id)}`}>{statusRu[row.priority_status]||row.priority_status_ru||row.priority_status}</Link></td>
        </tr>)}</tbody>
      </table></div>
      <div className={styles.pagination}>
        <Link className={offset<=0?styles.disabled:""} href={pageHref(sp,Math.max(0,offset-limit))}>← Назад</Link>
        <span className={styles.sub}>{catalog.total?`${offset+1}–${end} из ${n(catalog.total)}`:"0 записей"}</span>
        <Link className={end>=catalog.total?styles.disabled:""} href={pageHref(sp,offset+limit)}>Дальше →</Link>
      </div>
    </>}
  </>;
}
