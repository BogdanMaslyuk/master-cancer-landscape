import Link from "next/link";
import { apiGet } from "../../../lib/api";
import styles from "../../pharmacology.module.css";

export const dynamic = "force-dynamic";

type Payload = Record<string,any>;
function n(value:any){const x=Number(value);return Number.isFinite(x)?x.toLocaleString("ru-RU"):"—";}
function fmt(value:any,digits=3){if(value===null||value===undefined||value==="")return "—";const x=Number(value);return Number.isFinite(x)?x.toFixed(digits):String(value);}

export default async function CompoundPage({params}:{params:Promise<{id:string}>}){
  const {id}=await params;
  const compoundId=decodeURIComponent(id);
  const data=await apiGet<Payload>(`/api/compounds/${encodeURIComponent(compoundId)}?limit=80`);
  const identity=data.identity||{};
  const targets=(data.targets||[]) as any[];
  const responses=(data.response_examples||[]) as any[];
  const concordance=(data.target_concordance?.items||[]) as any[];
  return <>
    <div className="breadcrumbs"><Link href="/compounds">Вещества</Link><span>›</span><strong>{identity.preferred_name||compoundId}</strong></div>

    <section className={styles.detailHero}>
      <div className={styles.detailCard}>
        <div className="eyebrow">ВЕЩЕСТВО · ФАРМАКОЛОГИЧЕСКИЙ ПРОФИЛЬ</div>
        <h1>{identity.preferred_name||compoundId}</h1>
        <span className={styles.idPill}>{compoundId}</span>
        <div className={styles.section}>
          <div className={styles.definition}>
            <div><b>SMILES</b><code>{identity.canonical_smiles||"не указан"}</code></div>
            <div><b>InChIKey</b><span>{identity.inchikey||"—"}</span></div>
            <div><b>ChEMBL</b><span>{identity.chembl_id||"—"}</span></div>
            <div><b>Broad / PRISM</b><span>{identity.broad_id||"—"}</span></div>
            <div><b>PubChem CID</b><span>{identity.pubchem_cid||"—"}</span></div>
            <div><b>Источники ответа</b><span>{(identity.sources||[]).join(" · ")||"—"}</span></div>
          </div>
        </div>
      </div>
      <div className={styles.detailCard}>
        <div className="eyebrow">ПОКРЫТИЕ</div>
        <div className={styles.metricGrid}>
          <div className={styles.metric}><strong>{n(identity.models_n||data.response_summary?.models_n)}</strong><span>клеточных моделей</span></div>
          <div className={styles.metric}><strong>{n(identity.observations_n||data.response_summary?.observations_n)}</strong><span>наблюдений</span></div>
          <div className={styles.metric}><strong>{n(identity.targets_n||targets.length)}</strong><span>аннотированных мишеней</span></div>
          <div className={styles.metric}><strong>{n(identity.target_evidence_n)}</strong><span>строк доказательств</span></div>
        </div>
      </div>
    </section>

    <section className={styles.section}>
      <div className={styles.sectionHead}><div><div className="eyebrow">ВЕЩЕСТВО → БЕЛКОВАЯ МИШЕНЬ → ГЕН</div><h2>Аннотированные мишени</h2></div></div>
      {targets.length?<div className={styles.tableWrap}><table className={styles.table}>
        <thead><tr><th>Белковая мишень</th><th>Кодирующий ген</th><th>Тип действия / механизм</th><th>Доказательства</th><th>Источник</th></tr></thead>
        <tbody>{targets.map((row:any)=>{const protein=row.protein_preferred_name||null;return <tr key={row.target_gene}>
          <td><Link className={styles.primary} href={`/targets/${encodeURIComponent(row.target_gene)}`}>{protein||`Белковая мишень, связанная с ${row.target_gene}`}</Link><span className={styles.sub}>{row.uniprot_primary_accession?`UniProt ${row.uniprot_primary_accession}`:row.protein_mapping_status||"gene-mapped"}</span></td>
          <td><Link className={styles.primary} href={`/genes/${encodeURIComponent(row.target_gene)}`}>{row.target_gene} →</Link><span className={styles.sub}>исходная target_gene аннотация</span></td>
          <td>{(row.actions||[]).length?<div className={styles.chips}>{row.actions.slice(0,4).map((x:string)=><span className={styles.chip} key={x}>{x}</span>)}</div>:<span className={styles.sub}>тип действия не уточнён</span>}</td>
          <td><b>{n(row.evidence_rows_n)}</b><span className={styles.sub}>{(row.evidence_types||[]).join(" · ")||"аннотация мишени"}</span></td>
          <td>{(row.sources||[]).join(" · ")||"—"}</td>
        </tr>})}</tbody>
      </table></div>:<div className={styles.empty}>В текущем источнике для вещества нет аннотированной мишени.</div>}
      <div className={styles.note}>Белковое название — отдельное справочное сопоставление MCL через UniProtKB/Swiss-Prot. Исходная фармакологическая аннотация может оставаться только на уровне гена.</div>
    </section>

    <section className={styles.section}>
      <div className={styles.sectionHead}><div><div className="eyebrow">ВЕЩЕСТВО × КЛЕТОЧНАЯ МОДЕЛЬ</div><h2>Примеры экспериментальных ответов</h2></div></div>
      {responses.length?<div className={styles.tableWrap}><table className={styles.table}>
        <thead><tr><th>Модель</th><th>Опухолевый контекст</th><th>Источник</th><th>Результат</th><th>Условия</th></tr></thead>
        <tbody>{responses.map((row:any,index:number)=><tr key={`${row.model_id}-${index}`}>
          <td><Link className={styles.primary} href={`/models/${encodeURIComponent(row.model_id)}`}>{row.cell_line_name||row.model_id}</Link><span className={styles.sub}>{row.model_id}</span></td>
          <td>{row.mcl_cancer_name||"—"}<span className={styles.sub}>{row.mcl_organ_ru||""}</span></td>
          <td>{row.source||"—"}</td>
          <td><b>{row.endpoint||"показатель"}: {fmt(row.value)}</b><span className={styles.sub}>{row.unit||""}</span></td>
          <td>{row.dose!==null&&row.dose!==undefined?`${fmt(row.dose)} ${row.dose_unit||""}`:"—"}<span className={styles.sub}>{row.exposure_time_h?`${fmt(row.exposure_time_h,1)} ч`:""}</span></td>
        </tr>)}</tbody>
      </table></div>:<div className={styles.empty}>Нет нормализованных ответов клеточных моделей.</div>}
      <div className={styles.note}>Для PRISM LFC более отрицательное значение обычно соответствует более сильному снижению сигнала относительно контроля. В будущих источниках IC50, AUC, GI50 и другие показатели будут показаны раздельно и не станут автоматически сравниваться между собой.</div>
    </section>

    <section className={styles.section}>
      <div className={styles.sectionHead}><div><div className="eyebrow">ФАРМАКОЛОГИЯ ↔ CRISPR</div><h2>Согласованность заявленных мишеней</h2></div></div>
      {concordance.length?<div className={styles.tableWrap}><table className={styles.table}>
        <thead><tr><th>Белковая мишень / ген</th><th>Статус</th><th>Spearman ρ</th><th>Моделей</th><th>Интерпретация</th></tr></thead>
        <tbody>{concordance.map((row:any,index:number)=><tr key={`${row.target_gene}-${index}`}>
          <td><Link className={styles.primary} href={`/targets/${encodeURIComponent(row.target_gene)}`}>{row.protein_preferred_name||row.target_gene}</Link><span className={styles.sub}>ген {row.target_gene}{row.uniprot_primary_accession?` · UniProt ${row.uniprot_primary_accession}`:""}</span></td>
          <td>{row.concordance_label||"—"}</td><td>{fmt(row.spearman_rho)}</td><td>{n(row.models_n)}</td><td>{row.interpretation_ru||"—"}</td>
        </tr>)}</tbody>
      </table></div>:<div className={styles.callout}><b>Слой согласованности ещё не построен.</b> Он отдельно проверяет, совпадает ли профиль чувствительности клеток к веществу с CRISPR-зависимостью кодирующего гена по панели моделей.</div>}
    </section>

    <div className={styles.callout}><b>Граница интерпретации.</b> Аннотация target_gene, справочное белковое сопоставление и ответ конкретной клеточной линии — разные типы данных. Даже их согласованность не заменяет прямое подтверждение связывания и причинного механизма.</div>
  </>;
}
