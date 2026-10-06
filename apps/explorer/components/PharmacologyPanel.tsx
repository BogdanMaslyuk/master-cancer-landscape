import Link from "next/link";
import styles from "./PharmacologyPanel.module.css";

type Pharmacology = Record<string, any>;

function fmt(value:any,digits=3){if(value===null||value===undefined||value==="")return "—";const n=Number(value);return Number.isFinite(n)?n.toFixed(digits):String(value);}
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
function rankedTargets(row:any){
  const order:Record<string,number>={strong_dependency:0,dependency:1,weak_or_none:2,not_available:3};
  return [...(row.target_hypotheses||[])].sort((a:any,b:any)=>{
    const byClass=(order[a.crispr_support_level]??9)-(order[b.crispr_support_level]??9);
    if(byClass!==0)return byClass;
    const av=Number(a.gene_effect),bv=Number(b.gene_effect);
    if(Number.isFinite(av)&&Number.isFinite(bv))return av-bv;
    if(Number.isFinite(av))return -1;if(Number.isFinite(bv))return 1;return String(a.target_gene||"").localeCompare(String(b.target_gene||""));
  }).slice(0,6);
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
      <thead><tr><th>Вещество</th><th>Источник / анализ</th><th>Результат</th><th>Условия</th><th>Кандидатные белковые мишени в этой модели</th></tr></thead>
      <tbody>{items.map((row:any,index:number)=>{const targets=rankedTargets(row);return <tr key={row.observation_id||`${row.compound_id}-${index}`}>
        <td><Link href={`/compounds/${encodeURIComponent(row.compound_id||"")}`} className={styles.compound}>{row.preferred_name||row.compound_id}</Link><span className={styles.sub}>{row.compound_id}{row.chembl_id?` · ${row.chembl_id}`:""}</span>{row.canonical_smiles&&<span className={styles.sub}>SMILES доступен</span>}</td>
        <td><b>{row.source||"—"}</b><span className={styles.sub}>{row.assay_type||row.source_assay_id||"—"}</span></td>
        <td><b>{resultText(row)}</b><span className={styles.sub}>{row.endpoint||"endpoint не указан"}</span></td>
        <td>{row.dose!==null&&row.dose!==undefined?`доза ${fmt(row.dose)} ${row.dose_unit||""}`:"—"}<span className={styles.sub}>{row.exposure_time_h?`${fmt(row.exposure_time_h,1)} ч`:""}</span></td>
        <td>{targets.length?<div className={styles.targets}>{targets.map((t:any,j:number)=>{const protein=t.protein_preferred_name||null;return <Link href={`/targets/${encodeURIComponent(t.target_gene||"")}`} className={`${styles.target} ${supportClass(t.crispr_support_level)}`} key={`${t.target_gene}-${j}`}><b>{protein||`мишень гена ${t.target_gene||"?"}`}</b><span>ген {t.target_gene||"?"}{t.uniprot_primary_accession?` · UniProt ${t.uniprot_primary_accession}`:""}</span><span>{supportLabel(t.crispr_support_level)}{t.gene_effect!==null&&t.gene_effect!==undefined?` · GE ${fmt(t.gene_effect)}`:""}</span></Link>})}</div>:<span className={styles.sub}>мишень не аннотирована в текущем слое</span>}</td>
      </tr>})}</tbody>
    </table></div>:<div className={styles.empty}>Для этой модели нет нормализованных фармакологических наблюдений.</div>}
    <p className={styles.note}>Фармакологический источник обычно задаёт мишень через символ гена; MCL отдельно добавляет название соответствующего reviewed-белка из UniProt, если сопоставление однозначно. CRISPR относится к кодирующему гену. Эти слои поддерживают механистическую гипотезу, но не доказывают связывание или конкретную изоформу.</p>
  </>;
}
