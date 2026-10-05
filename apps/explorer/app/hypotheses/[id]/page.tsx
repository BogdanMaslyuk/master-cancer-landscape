import Link from "next/link";
import { apiGet } from "../../../lib/api";
import styles from "../../pharmacology.module.css";

export const dynamic = "force-dynamic";

type Payload=Record<string,any>;
function n(value:any){const x=Number(value);return Number.isFinite(x)?x.toLocaleString("ru-RU"):"—";}
function pct(value:any){if(value===null||value===undefined||value==="")return "—";const x=Number(value);return Number.isFinite(x)?`${Math.round(x*100)}%`:"—";}
function fmt(value:any,digits=3){if(value===null||value===undefined||value==="")return "—";const x=Number(value);return Number.isFinite(x)?x.toFixed(digits):String(value);}
const statusRu:Record<string,string>={priority_for_in_vitro:"Приоритет для in vitro проверки",supported_hypothesis:"Поддержанная гипотеза",exploratory_hypothesis:"Исследовательская гипотеза",insufficient_evidence:"Недостаточно данных"};
const axisRu:Record<string,string>={strong:"сильная",supportive:"поддерживает",weak:"слабая",sparse:"мало данных",conflicting:"противоречит",uncertain:"неопределённо",not_assessable:"не сопоставимо",not_available:"нет данных",concordant_enrichment:"согласованное обогащение",partial_enrichment:"частичное обогащение",weak_enrichment:"слабое обогащение",not_enriched:"не обогащено",limited_sample:"малая группа",not_assessed:"не оценено"};

function ModelTable({rows,title}:{rows:any[];title:string}){
  return <section className={styles.section}><div className={styles.sectionHead}><div><h2>{title}</h2></div></div>
    {rows.length?<div className={styles.tableWrap}><table className={styles.table}><thead><tr><th>Модель</th><th>Опухоль</th><th>Ответ на вещество</th><th>Относительная чувствительность</th><th>Gene Effect</th><th>Почему выбрана</th></tr></thead><tbody>
      {rows.map((row:any)=><tr key={`${row.model_role}-${row.model_id}`}>
        <td><Link className={styles.primary} href={`/models/${encodeURIComponent(row.model_id)}`}>{row.cell_line_name||row.model_id}</Link><span className={styles.sub}>{row.model_id}</span></td>
        <td>{row.mcl_cancer_name||"—"}<span className={styles.sub}>{row.oncotree_subtype||""}</span></td>
        <td><b>LFC {fmt(row.response_value)}</b></td>
        <td><b>{pct(row.relative_sensitivity)}</b><span className={styles.sub}>внутри профиля этого вещества</span></td>
        <td><b>{fmt(row.gene_effect)}</b></td>
        <td>{row.role_reason_ru||"—"}</td>
      </tr>)}
    </tbody></table></div>:<div className={styles.empty}>Подходящие модели этого типа не найдены в текущем пересечении данных.</div>}
  </section>;
}

export default async function HypothesisPage({params}:{params:Promise<{id:string}>}){
  const {id}=await params;
  const hypothesisId=decodeURIComponent(id);
  const data=await apiGet<Payload>(`/api/hypotheses/${encodeURIComponent(hypothesisId)}`);
  const h=data.hypothesis||{};
  const roles=data.models_by_role||{};
  const reasons=(h.priority_reasons||[]) as string[];
  const gaps=(h.evidence_gaps||[]) as string[];
  const lab=(data.laboratory_route||[]) as any[];
  const falsification=(data.falsification_criteria||[]) as string[];
  return <>
    <div className="breadcrumbs"><Link href="/hypotheses">Исследовательские гипотезы</Link><span>›</span><strong>{h.preferred_name||hypothesisId}</strong></div>

    <section className={styles.detailHero}>
      <div className={styles.detailCard}>
        <div className="eyebrow">ВЕЩЕСТВО → МИШЕНЬ → ОПУХОЛЬ → ЭКСПЕРИМЕНТ</div>
        <h1>{h.preferred_name||h.compound_id}</h1>
        <span className={styles.idPill}>{statusRu[h.priority_status]||h.priority_status_ru||h.priority_status||"—"}</span>
        <div className={styles.section}>
          <div className={styles.definition}>
            <div><b>Белковая мишень</b><span><Link href={`/targets/${encodeURIComponent(h.target_gene||"")}`}>{h.protein_preferred_name||h.target_gene||"—"}</Link></span></div>
            <div><b>Кодирующий ген</b><span><Link href={`/genes/${encodeURIComponent(h.target_gene||"")}`}>{h.target_gene||"—"}</Link></span></div>
            <div><b>UniProt</b><span>{h.uniprot_primary_accession||"не разрешён однозначно"}</span></div>
            <div><b>Опухолевый контекст</b><span>{h.mcl_cancer_name||"—"}<br/>{h.mcl_organ_ru||""}</span></div>
          </div>
        </div>
      </div>
      <div className={styles.detailCard}>
        <div className="eyebrow">ЦЕЛОСТНОСТЬ ДОКАЗАТЕЛЬСТВ</div>
        <div className={styles.metricGrid}>
          <div className={styles.metric}><strong>{axisRu[h.phenotype_axis]||h.phenotype_axis||"—"}</strong><span>фармакологический фенотип</span></div>
          <div className={styles.metric}><strong>{axisRu[h.dependency_axis]||h.dependency_axis||"—"}</strong><span>CRISPR-зависимость</span></div>
          <div className={styles.metric}><strong>{axisRu[h.mechanism_axis]||h.mechanism_axis||"—"}</strong><span>профильная согласованность</span></div>
          <div className={styles.metric}><strong>{axisRu[h.specificity_axis]||h.specificity_axis||"—"}</strong><span>обогащение в контексте</span></div>
        </div>
      </div>
    </section>

    <section className={styles.section}>
      <div className={styles.sectionHead}><div><div className="eyebrow">ПОЧЕМУ ЭТА ГИПОТЕЗА ВЫДЕЛЕНА</div><h2>Раздельные оси доказательств</h2></div></div>
      <div className={styles.definition}>
        <div><b>Фенотип</b><span>{n(h.sensitive_models_n)} / {n(h.models_n)} моделей в чувствительной четверти ({pct(h.sensitive_fraction)}). Медианный LFC: {fmt(h.median_response)}.</span></div>
        <div><b>CRISPR</b><span>{n(h.dependency_models_n)} / {n(h.crispr_models_n)} моделей с Gene Effect ≤ −0,5 ({pct(h.dependency_fraction_in_cancer)}). Медианный GE: {fmt(h.median_gene_effect)}.</span></div>
        <div><b>Совместная поддержка</b><span>{n(h.joint_support_models_n)} моделей одновременно чувствительны к веществу и зависимы от кодирующего гена.</span></div>
        <div><b>Профиль вещество ↔ мишень</b><span>{h.concordance_label||"не оценён"}{h.spearman_rho!==null&&h.spearman_rho!==undefined?` · Spearman ρ ${fmt(h.spearman_rho)}`:""}</span></div>
        <div><b>Обогащение чувствительности</b><span>{h.sensitivity_enrichment===null||h.sensitivity_enrichment===undefined?"—":`${h.sensitivity_enrichment>=0?"+":""}${Math.round(Number(h.sensitivity_enrichment)*100)} п.п.`} относительно общей панели вещества.</span></div>
        <div><b>Обогащение зависимости</b><span>{h.dependency_enrichment===null||h.dependency_enrichment===undefined?"—":`${h.dependency_enrichment>=0?"+":""}${Math.round(Number(h.dependency_enrichment)*100)} п.п.`} относительно глобальной CRISPR-зависимости от гена.</span></div>
      </div>
      {reasons.length>0&&<div className={styles.callout} style={{marginTop:12}}><b>Что поддерживает гипотезу:</b><ul>{reasons.map((x:string)=><li key={x}>{x}</li>)}</ul></div>}
    </section>

    <ModelTable rows={(roles.positive||[]) as any[]} title="Положительные модели для проверки"/>
    <ModelTable rows={(roles.negative_same_cancer||[]) as any[]} title="Отрицательные контроли того же опухолевого контекста"/>

    {((roles.discordant_sensitive_without_dependency||[]).length>0||(roles.discordant_dependency_without_sensitivity||[]).length>0)&&<section className={styles.section}>
      <div className={styles.sectionHead}><div><div className="eyebrow">НЕ ИГНОРИРОВАТЬ ПРОТИВОРЕЧИЯ</div><h2>Модели, которые могут опровергнуть простую версию механизма</h2></div></div>
      <div className={styles.callout}>Чувствительность без CRISPR-зависимости может указывать на альтернативную мишень или полифармакологию. CRISPR-зависимость без чувствительности к веществу может указывать на недостаточное фармакологическое воздействие, контекстную резистентность или несоответствие механизма.</div>
      <ModelTable rows={(roles.discordant_sensitive_without_dependency||[]) as any[]} title="Чувствительны без ожидаемой CRISPR-зависимости"/>
      <ModelTable rows={(roles.discordant_dependency_without_sensitivity||[]) as any[]} title="CRISPR-зависимы, но не чувствительны к веществу"/>
    </section>}

    <section className={styles.section}>
      <div className={styles.sectionHead}><div><div className="eyebrow">ЧТО ДЕЛАТЬ В ЛАБОРАТОРИИ</div><h2>Маршрут экспериментальной проверки</h2></div></div>
      <div className={styles.definition}>{lab.map((row:any)=><div key={row.stage}><b>{row.stage}</b><span>{row.goal}</span></div>)}</div>
    </section>

    <section className={styles.section}>
      <div className={styles.sectionHead}><div><div className="eyebrow">КРИТЕРИИ ОПРОВЕРЖЕНИЯ</div><h2>Что заставит отказаться от гипотезы или пересмотреть механизм</h2></div></div>
      <div className={styles.callout}><ul>{falsification.map((x:string)=><li key={x}>{x}</li>)}</ul></div>
    </section>

    <section className={styles.section}>
      <div className={styles.sectionHead}><div><div className="eyebrow">ЧЕГО ПОКА НЕ ХВАТАЕТ</div><h2>Границы текущего статуса</h2></div></div>
      <div className={styles.callout}><b>Статус ещё не учитывает молекулярную генетику модели и нормальные ткани.</b><ul>{gaps.map((x:string)=><li key={x}>{x}</li>)}</ul></div>
    </section>
  </>;
}
