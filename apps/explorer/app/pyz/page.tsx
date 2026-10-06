import Link from "next/link";
import { apiGet, formatNumber } from "../../lib/api";
import styles from "./pyz.module.css";

export const dynamic = "force-dynamic";

type Params = Record<string,string|string[]|undefined>;
type Payload = Record<string,any>;
function one(v:string|string[]|undefined){return Array.isArray(v)?v[0]:v;}
function riskRu(value:any){return value==="high_predicted_risk"?"высокий предиктивный риск":"требует проверки";}

export default async function PyzCatalogPage({searchParams}:{searchParams:Promise<Params>}){
  const sp=await searchParams;
  const q=one(sp.q)||"";
  const target=one(sp.target_gene)||"";
  const shortlist=one(sp.shortlist_only)==="true"||one(sp.shortlist_only)==="1";
  const query=new URLSearchParams({limit:"100",offset:"0"});
  if(q)query.set("q",q);
  if(target)query.set("target_gene",target);
  if(shortlist)query.set("shortlist_only","true");
  const [summary,catalog]=await Promise.all([
    apiGet<Payload>("/api/pyz/summary"),
    apiGet<Payload>(`/api/pyz?${query.toString()}`),
  ]);
  const admet=summary.admet||{};
  const items=(catalog.items||[]) as any[];

  return <>
    <section className={styles.pageHeader}>
      <div>
        <div className="eyebrow">СОБСТВЕННАЯ ХИМИЧЕСКАЯ БИБЛИОТЕКА · PYZ-001…PYZ-065</div>
        <h1>Молекулы PYZ</h1>
        <p>Единая точка входа от нашей структуры к экспериментальным лигандам мишени, опухолевому контексту, клеточным моделям, лабораторной выполнимости и ADMET. Клеточные культуры здесь означают, где рационально проверять гипотезу, а не где PYZ уже доказанно активна.</p>
      </div>
      <div className={styles.heroStats}>
        <div className={styles.heroStat}><strong>{formatNumber(summary.compounds_n,0)}</strong><span>собственных молекул</span></div>
        <div className={styles.heroStat}><strong>{formatNumber(summary.targets_n,0)}</strong><span>core-мишеней</span></div>
        <div className={styles.heroStat}><strong>{formatNumber(summary.shortlist_rows_n,0)}</strong><span>пар Tanimoto ≥ 0,35</span></div>
        <div className={styles.heroStat}><strong>{(summary.shortlist_targets||[]).join(" · ")||"—"}</strong><span>мишени текущего shortlist</span></div>
      </div>
    </section>

    <div className={styles.notice}><b>Научная граница:</b> {summary.interpretation_ru}</div>

    <section className={styles.section}>
      <div className={styles.sectionHead}>
        <div><div className="eyebrow">ADMET · СНИМОК SYL-RPT-2026-022 v1.6</div><h2>Разработческие ограничения серии</h2></div>
        <Link className={styles.secondaryButton} href="/pyz/matrix">Открыть матрицу PYZ × мишень →</Link>
      </div>
      <div className={styles.admetBar}>
        <div className={styles.admetMetric}><b>{formatNumber(admet.ames_positive_n,0)}/{formatNumber(admet.compounds_n,0)}</b><span>Ames+ предиктивный сигнал</span></div>
        <div className={styles.admetMetric}><b>{formatNumber(admet.dili_positive_n,0)}/{formatNumber(admet.compounds_n,0)}</b><span>DILI+ предиктивный сигнал</span></div>
        <div className={styles.admetMetric}><b>{formatNumber(admet.outside_applicability_domain_n,0)}/{formatNumber(admet.compounds_n,0)}</b><span>вне области применимости части моделей</span></div>
        <div className={styles.admetMetric}><b>{formatNumber(admet.high_integrated_risk_n,0)}</b><span>в красной очереди экспериментальной токсикологии</span></div>
      </div>
      <div className={styles.warning} style={{marginTop:10}}>{admet.interpretation_ru||"ADMET-слой пока не материализован."}</div>
    </section>

    <section className={styles.section}>
      <div className={styles.sectionHead}><div><div className="eyebrow">65 МОЛЕКУЛ · 845 ПАР С 13 МИШЕНЯМИ</div><h2>Каталог собственной серии</h2></div><div className={styles.sectionCopy}>Сначала показаны соединения, имеющие хотя бы одну пару Tanimoto ≥ 0,35.</div></div>
      <form action="/pyz" method="get" className={styles.toolbar}>
        <label className={styles.field}><span>Поиск</span><input name="q" defaultValue={q} placeholder="PYZ-035, антипиринат, SMILES…"/></label>
        <label className={styles.field}><span>Мишень shortlist</span><select name="target_gene" defaultValue={target}><option value="">все</option>{(summary.target_genes||[]).map((gene:string)=><option key={gene} value={gene}>{gene}</option>)}</select></label>
        <label className={styles.field}><span>Показать</span><select name="shortlist_only" defaultValue={shortlist?"true":"false"}><option value="false">все 65</option><option value="true">только молекулы с парой ≥0,35</option></select></label>
        <button className={styles.button} type="submit">Применить</button><Link className={styles.reset} href="/pyz">Сбросить</Link>
      </form>

      {items.length===0?<div className={styles.empty}>По выбранным условиям молекулы не найдены.</div>:
      <div className={styles.catalog}>{items.map((row:any)=>
        <Link href={`/pyz/${encodeURIComponent(row.own_compound_id)}`} className={styles.compoundCard} key={row.own_compound_id}>
          <div className={styles.cardTop}><div><div className={styles.compoundId}>{row.own_compound_id}</div><div className={styles.sectionCopy}>{row.chemotype_ru} · {row.ring_variant}</div></div><span className={row.integrated_risk==="high_predicted_risk"?styles.badgeRisk:styles.badgeSoft}>{riskRu(row.integrated_risk)}</span></div>
          <div className={styles.badges}><span className={styles.badge}>тир {row.tier}</span>{(row.shortlist_target_genes||[]).map((gene:string)=><span className={styles.badgePriority} key={gene}>{gene}</span>)}</div>
          <div className={styles.cardMetric}><span>Лучший target-сосед</span><b>{row.best_target_gene||"—"}</b></div>
          <div className={styles.cardMetric}><span>Лучший Tanimoto</span><b>{row.best_tanimoto!==null&&row.best_tanimoto!==undefined?Number(row.best_tanimoto).toFixed(3):"—"}</b></div>
          <div className={styles.cardMetric}><span>QED</span><b>{row.qed!==null&&row.qed!==undefined?Number(row.qed).toFixed(2):"—"}</b></div>
          <div className={styles.cardFooter}>Открыть молекулу →</div>
        </Link>
      )}</div>}
    </section>
  </>;
}
