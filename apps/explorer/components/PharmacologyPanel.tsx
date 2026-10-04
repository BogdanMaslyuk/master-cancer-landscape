import Link from "next/link";
import styles from "./PharmacologyPanel.module.css";

type Pharmacology = Record<string, any>;

function fmt(value:any,digits=3){const n=Number(value);return Number.isFinite(n)?n.toFixed(digits):"—";}
function resultText(row:any){
  if(row.endpoint && row.value!==null && row.value!==undefined && row.value!=="") return `${row.endpoint}: ${fmt(row.value)}${row.unit?` ${row.unit}`:""}`;
  for(const key of ["ic50","ec50","gi50","auc","viability"]){
    const value=row[key]; if(value!==null&&value!==undefined&&value!=="") return `${key.toUpperCase()}: ${fmt(value)}${row.unit?` ${row.unit}`:""}`;
  }
  return "результат сохранён в исходной записи";
}
function supportClass(level:string){
  if(level==="strong_dependency") return styles.strong;
  if(level==="dependency") return styles.support;
  return styles.weak;
}
function supportLabel(level:string){
  return ({strong_dependency:"сильная CRISPR-зависимость",dependency:"CRISPR-зависимость",weak_or_none:"слабая / нет зависимости",not_available:"CRISPR нет"} as Record<string,string>)[level]||level||"—";
}

export default function PharmacologyPanel({payload}:{payload:Pharmacology}){
  if(!payload?.available){
    return <div className={styles.empty}>
      <b>Фармакологический слой для этой модели пока не заполнен.</b><br/>
      MCL уже готов хранить результаты PRISM, GDSC, CTRP/CTD², ChEMBL и PubChem как отдельные экспериментальные наблюдения и связывать их с известными мишенями и CRISPR-профилем.
      {payload?.build_command&&<code>{payload.build_command}</code>}
    </div>;
  }
  const items=(payload.items||[]) as any[];
  return <>
    <div className={styles.summary}>
      <div className={styles.metric}><strong>{Number(payload.compounds_n||0).toLocaleString("ru-RU")}</strong><span>веществ тестировалось</span></div>
      <div className={styles.metric}><strong>{Number(payload.observations_n||0).toLocaleString("ru-RU")}</strong><span>экспериментальных наблюдений</span></div>
      <div className={styles.metric}><strong>{(payload.sources||[]).length}</strong><span>источников данных</span></div>
    </div>
    {(payload.sources||[]).length>0&&<div className={styles.sources}>{payload.sources.map((s:string)=><span className={styles.source} key={s}>{s}</span>)}</div>}
    {items.length?<div className={styles.tableWrap}><table className={styles.table}>
      <thead><tr><th>Вещество</th><th>Источник / анализ</th><th>Результат</th><th>Условия</th><th>Мишени и согласованность с CRISPR</th></tr></thead>
      <tbody>{items.map((row:any,index:number)=><tr key={row.observation_id||`${row.compound_id}-${index}`}>
        <td><span className={styles.compound}>{row.preferred_name||row.compound_id}</span><span className={styles.sub}>{row.compound_id}{row.chembl_id?` · ${row.chembl_id}`:""}</span></td>
        <td><b>{row.source||"—"}</b><span className={styles.sub}>{row.assay_type||row.source_assay_id||"—"}</span></td>
        <td><b>{resultText(row)}</b><span className={styles.sub}>{row.endpoint||"endpoint не указан"}</span></td>
        <td>{row.dose!==null&&row.dose!==undefined?`доза ${fmt(row.dose)} ${row.dose_unit||""}`:"—"}<span className={styles.sub}>{row.exposure_time_h?`${fmt(row.exposure_time_h,1)} ч`:""}</span></td>
        <td>{(row.target_hypotheses||[]).length?<div className={styles.targets}>{row.target_hypotheses.map((t:any,j:number)=><Link href={`/genes/${encodeURIComponent(t.target_gene||"")}`} className={`${styles.target} ${supportClass(t.crispr_support_level)}`} key={`${t.target_gene}-${j}`}><b>{t.target_gene||"?"}</b><span>{t.action||t.evidence_type||"мишень"}</span><span>{supportLabel(t.crispr_support_level)}{t.gene_effect!==null&&t.gene_effect!==undefined?` · GE ${fmt(t.gene_effect)}`:""}</span></Link>)}</div>:<span className={styles.sub}>мишень не аннотирована в текущем слое</span>}</td>
      </tr>)}</tbody>
    </table></div>:<div className={styles.empty}>Для этой модели нет нормализованных фармакологических наблюдений.</div>}
    <p className={styles.note}>Важно: чувствительность клетки к веществу, аннотация мишени и CRISPR-зависимость — разные типы доказательств. Их совпадение поддерживает механизм, но само по себе его не доказывает. IC50, AUC, GI50 и viability не объединяются в одну общую шкалу.</p>
  </>;
}
