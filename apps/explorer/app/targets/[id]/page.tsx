import Link from "next/link";
import { apiGet } from "../../../lib/api";
import styles from "../../pharmacology.module.css";

export const dynamic = "force-dynamic";

type Payload = Record<string,any>;
function n(value:any){const x=Number(value);return Number.isFinite(x)?x.toLocaleString("ru-RU"):"—";}
function pct(value:any){if(value===null||value===undefined||value==="")return "—";const x=Number(value);return Number.isFinite(x)?`${Math.round(x*100)}%`:"—";}
function fmt(value:any,digits=3){if(value===null||value===undefined||value==="")return "—";const x=Number(value);return Number.isFinite(x)?x.toFixed(digits):String(value);}

export default async function TargetPage({params}:{params:Promise<{id:string}>}){
  const {id}=await params;
  const targetId=decodeURIComponent(id).toUpperCase();
  const data=await apiGet<Payload>(`/api/targets/${encodeURIComponent(targetId)}?limit=120`);
  const identity=data.identity||{};
  const compounds=(data.compounds||[]) as any[];
  const models=(data.model_examples||[]) as any[];
  const concordance=(data.target_concordance?.items||[]) as any[];
  return <>
    <div className="breadcrumbs"><Link href="/targets">Белковые мишени</Link><span>›</span><strong>{identity.target_gene||targetId}</strong></div>

    <section className={styles.detailHero}>
      <div className={styles.detailCard}>
        <div className="eyebrow">БЕЛКОВАЯ МИШЕНЬ · GENE-MAPPED TARGET</div>
        <h1>{identity.target_gene||targetId}</h1>
        <span className={styles.idPill}>{data.resolution||"gene_mapped_protein_target"}</span>
        <div className={styles.section}>
          <div className={styles.callout}><b>Кодирующий ген:</b> <Link href={`/genes/${encodeURIComponent(data.coding_gene||targetId)}`}>{data.coding_gene||targetId} → открыть генетическую карточку</Link><br/>{data.resolution_note_ru}</div>
        </div>
        {(identity.actions||[]).length>0&&<div className={styles.section}><div className={styles.chips}>{identity.actions.slice(0,8).map((x:string)=><span className={styles.chip} key={x}>{x}</span>)}</div></div>}
      </div>
      <div className={styles.detailCard}>
        <div className="eyebrow">ПОКРЫТИЕ</div>
        <div className={styles.metricGrid}>
          <div className={styles.metric}><strong>{n(identity.compounds_n)}</strong><span>связанных веществ</span></div>
          <div className={styles.metric}><strong>{n(identity.models_n)}</strong><span>моделей с фармакологией</span></div>
          <div className={styles.metric}><strong>{pct(identity.dependency_fraction)}</strong><span>CRISPR-зависимых среди покрытых</span></div>
          <div className={styles.metric}><strong>{n(identity.target_evidence_n)}</strong><span>строк доказательств</span></div>
        </div>
      </div>
    </section>

    <section className={styles.section}>
      <div className={styles.sectionHead}><div><div className="eyebrow">МИШЕНЬ → ВЕЩЕСТВА</div><h2>Какие соединения аннотированы для этой мишени?</h2></div></div>
      {compounds.length?<div className={styles.tableWrap}><table className={styles.table}>
        <thead><tr><th>Вещество</th><th>SMILES</th><th>Тип действия</th><th>Модели</th><th>CRISPR-поддержка в покрытых моделях</th></tr></thead>
        <tbody>{compounds.map((row:any)=><tr key={row.compound_id}>
          <td><Link className={styles.primary} href={`/compounds/${encodeURIComponent(row.compound_id)}`}>{row.preferred_name||row.compound_id}</Link><span className={styles.sub}>{row.compound_id}</span></td>
          <td>{row.canonical_smiles?<div className={styles.smiles}>{String(row.canonical_smiles).slice(0,86)}{String(row.canonical_smiles).length>86?"…":""}</div>:<span className={styles.sub}>не указан</span>}</td>
          <td>{(row.actions||[]).length?<div className={styles.chips}>{row.actions.slice(0,3).map((x:string)=><span className={styles.chip} key={x}>{x}</span>)}</div>:<span className={styles.sub}>не уточнено</span>}</td>
          <td><b>{n(row.models_n)}</b></td>
          <td><b>{n(row.dependent_models_n)}</b><span className={styles.sub}>GE ≤ −0,5 · сильных {n(row.strong_dependency_models_n)}</span></td>
        </tr>)}</tbody>
      </table></div>:<div className={styles.empty}>Связанные вещества не найдены.</div>}
    </section>

    <section className={styles.section}>
      <div className={styles.sectionHead}><div><div className="eyebrow">МИШЕНЬ × КЛЕТОЧНАЯ МОДЕЛЬ</div><h2>Примеры CRISPR-зависимости среди фармакологически покрытых моделей</h2></div></div>
      {models.length?<div className={styles.tableWrap}><table className={styles.table}>
        <thead><tr><th>Модель</th><th>Опухолевый контекст</th><th>Gene Effect</th><th>Категория</th></tr></thead>
        <tbody>{models.map((row:any)=><tr key={row.model_id}>
          <td><Link className={styles.primary} href={`/models/${encodeURIComponent(row.model_id)}`}>{row.cell_line_name||row.model_id}</Link><span className={styles.sub}>{row.model_id}</span></td>
          <td>{row.mcl_cancer_name||"—"}<span className={styles.sub}>{row.mcl_organ_ru||""}</span></td>
          <td><b>{fmt(row.gene_effect)}</b></td><td>{row.crispr_support_level||"—"}</td>
        </tr>)}</tbody>
      </table></div>:<div className={styles.empty}>Нет доступных пересечений фармакологического покрытия с CRISPR для этой мишени.</div>}
      <div className={styles.note}>Эти модели не доказывают действие конкретного вещества через данную мишень. Они показывают только функциональную CRISPR-зависимость от кодирующего гена среди моделей, где тестировались связанные вещества.</div>
    </section>

    <section className={styles.section}>
      <div className={styles.sectionHead}><div><div className="eyebrow">ПРОВЕРКА МЕХАНИЗМА</div><h2>Вещество × мишень: согласованность с CRISPR</h2></div></div>
      {concordance.length?<div className={styles.tableWrap}><table className={styles.table}>
        <thead><tr><th>Вещество</th><th>Статус</th><th>Spearman ρ</th><th>Моделей</th><th>Интерпретация</th></tr></thead>
        <tbody>{concordance.map((row:any,index:number)=><tr key={`${row.compound_id}-${index}`}>
          <td><Link className={styles.primary} href={`/compounds/${encodeURIComponent(row.compound_id||"")}`}>{row.preferred_name||row.compound_id||"вещество"}</Link></td>
          <td>{row.concordance_label||"—"}</td><td>{fmt(row.spearman_rho)}</td><td>{n(row.models_n)}</td><td>{row.interpretation_ru||"—"}</td>
        </tr>)}</tbody>
      </table></div>:<div className={styles.callout}><b>Слой профильной согласованности пока не построен.</b> Следующий этап сможет проверять, связана ли чувствительность к конкретному веществу с CRISPR-зависимостью этой мишени по панели клеточных моделей.</div>}
    </section>
  </>;
}
