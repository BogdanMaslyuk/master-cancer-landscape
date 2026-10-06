import Link from "next/link";
import { apiGet } from "../../../lib/api";
import styles from "../../pharmacology.module.css";

export const dynamic = "force-dynamic";

type Payload=Record<string,any>;
function n(value:any){const x=Number(value);return Number.isFinite(x)?x.toLocaleString("ru-RU"):"—";}
function pct(value:any){if(value===null||value===undefined||value==="")return "—";const x=Number(value);return Number.isFinite(x)?`${Math.round(x*100)}%`:"—";}
function fmt(value:any,digits=3){if(value===null||value===undefined||value==="")return "—";const x=Number(value);return Number.isFinite(x)?x.toFixed(digits):String(value);}
function pp(value:any){if(value===null||value===undefined||value==="")return "—";const x=Number(value);return Number.isFinite(x)?`${x>=0?"+":""}${Math.round(x*100)} п.п.`:"—";}
function mutationText(row:any){
  const parts:string[]=[];
  if(row.has_hotspot)parts.push("hotspot");
  if(row.has_likely_lof)parts.push("LikelyLoF");
  const changes=(row.protein_changes||[]) as string[];
  if(changes.length)parts.push(changes.slice(0,2).join(", "));
  return parts.length?parts.join(" · "):"нет выделенной аннотации";
}
const statusRu:Record<string,string>={priority_for_in_vitro:"Приоритет для in vitro проверки",supported_hypothesis:"Поддержанная гипотеза",exploratory_hypothesis:"Исследовательская гипотеза",insufficient_evidence:"Недостаточно данных"};
const axisRu:Record<string,string>={
  strong:"сильная",supportive:"поддерживает",weak:"слабая",sparse:"мало данных",conflicting:"противоречит",
  uncertain:"неопределённо",not_assessable:"не сопоставимо",not_available:"нет данных",concordant_enrichment:"согласованное обогащение",
  partial_enrichment:"частичное обогащение",weak_enrichment:"слабое обогащение",not_enriched:"не обогащено",limited_sample:"малая группа",
  not_assessed:"не оценено",mutation_context_present_direction_unresolved:"есть мутационный контекст",expression_context_available:"есть экспрессия",
  omics_context_available:"есть молекулярные данные",
};

function ModelTable({rows,title,isV2}:{rows:any[];title:string;isV2:boolean}){
  return <section className={styles.section}><div className={styles.sectionHead}><div><h2>{title}</h2></div></div>
    {rows.length?<div className={styles.tableWrap}><table className={styles.table}><thead><tr><th>Модель</th><th>Опухоль</th><th>Фармакология</th><th>CRISPR</th>{isV2&&<th>Молекулярный контекст мишени</th>}<th>Почему выбрана</th></tr></thead><tbody>
      {rows.map((row:any)=><tr key={`${row.model_role}-${row.model_id}`}>
        <td><Link className={styles.primary} href={`/models/${encodeURIComponent(row.model_id)}`}>{row.cell_line_name||row.model_id}</Link><span className={styles.sub}>{row.model_id}</span></td>
        <td>{row.mcl_cancer_name||"—"}<span className={styles.sub}>{row.oncotree_subtype||""}</span></td>
        <td><b>LFC {fmt(row.response_value)}</b><span className={styles.sub}>относительная чувствительность {pct(row.relative_sensitivity)}</span>{isV2&&<span className={styles.sub}>{row.absolute_active?"проходит LFC ≤ −1":"не проходит LFC ≤ −1"}</span>}</td>
        <td>{isV2?<><b>P(dep) {fmt(row.dependency_probability)}</b><span className={styles.sub}>{row.dependency_call===true?"dependent (>0,5)":row.dependency_call===false?"не dependent (≤0,5)":"нет probability"}</span><span className={styles.sub}>Gene Effect {fmt(row.gene_effect)} — непрерывный фенотип</span></>:<><b>Gene Effect {fmt(row.gene_effect)}</b></>}</td>
        {isV2&&<td><b>RNA {fmt(row.expression_log2_tpm1)}</b><span className={styles.sub}>log₂(TPM+1) · перцентиль в опухоли {pct(row.expression_percentile_in_cancer)}</span><span className={styles.sub}>relative CN {fmt(row.copy_number_relative)} · перцентиль {pct(row.copy_number_percentile_in_cancer)}</span><span className={styles.sub}>{mutationText(row)}</span></td>}
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
  const isV2=data.candidate_version==="v2";
  const reasons=(h.priority_reasons||[]) as string[];
  const gaps=(h.evidence_gaps||[]) as string[];
  const lab=(data.laboratory_route||[]) as any[];
  const falsification=(data.falsification_criteria||[]) as string[];
  return <>
    <div className="breadcrumbs"><Link href="/hypotheses">Исследовательские гипотезы</Link><span>›</span><strong>{h.preferred_name||hypothesisId}</strong></div>

    <section className={styles.detailHero}>
      <div className={styles.detailCard}>
        <div className="eyebrow">ВЕЩЕСТВО → МИШЕНЬ → ОПУХОЛЬ → ЭКСПЕРИМЕНТ · {isV2?"CANDIDATE v2":"CANDIDATE v1"}</div>
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
          <div className={styles.metric}><strong>{axisRu[h.molecular_context_axis]||h.molecular_context_axis||"—"}</strong><span>молекулярный контекст</span></div>
        </div>
      </div>
    </section>

    <section className={styles.section}>
      <div className={styles.sectionHead}><div><div className="eyebrow">ПОЧЕМУ ЭТА ГИПОТЕЗА ВЫДЕЛЕНА</div><h2>Раздельные оси доказательств</h2></div></div>
      <div className={styles.definition}>
        {isV2?<>
          <div><b>Абсолютная активность</b><span>{n(h.active_models_n)} / {n(h.models_n)} моделей имеют PRISM LFC ≤ −1 ({pct(h.active_fraction)}). Это операционный порог первичного single-dose скрининга, а не клиническая концентрация эффективности.</span></div>
          <div><b>Активность + селективность</b><span>{n(h.sensitive_models_n)} / {n(h.models_n)} одновременно проходят LFC ≤ −1 и входят в наиболее чувствительную четверть ({pct(h.sensitive_fraction)}).</span></div>
          <div><b>Бинарная CRISPR-зависимость</b><span>{n(h.dependency_models_n)} / {n(h.dependency_measured_models_n)} моделей имеют Probability of Dependency &gt; 0,5 ({pct(h.dependency_fraction_in_cancer)}). Медианная P(dep): {fmt(h.median_dependency_probability)}.</span></div>
          <div><b>Непрерывный CRISPR-фенотип</b><span>Медианный Chronos Gene Effect: {fmt(h.median_gene_effect)}. Gene Effect показан как сила loss-of-function фенотипа, а не как бинарный вызов зависимости.</span></div>
        </>:<>
          <div><b>Фенотип</b><span>{n(h.sensitive_models_n)} / {n(h.models_n)} моделей в чувствительной четверти ({pct(h.sensitive_fraction)}). Медианный LFC: {fmt(h.median_response)}.</span></div>
          <div><b>CRISPR v1</b><span>{n(h.dependency_models_n)} / {n(h.crispr_models_n)} моделей по описательному порогу Gene Effect ≤ −0,5. Для лабораторного выбора предпочтительнее v2.</span></div>
        </>}
        <div><b>Совместная поддержка</b><span>{n(h.joint_support_models_n)} моделей одновременно удовлетворяют фармакологическому и CRISPR-критериям текущей версии.</span></div>
        <div><b>Профиль вещество ↔ мишень</b><span>{h.concordance_label||"не оценён"}{h.spearman_rho!==null&&h.spearman_rho!==undefined?` · ρ ${fmt(h.spearman_rho)}`:""}{h.primary_rho_method?` · ${h.primary_rho_method==="lineage_fixed_effect_rank_residual"?"ранговая корреляция после центрирования по lineage":"сырая Spearman-корреляция"}`:""}</span></div>
        <div><b>Обогащение активности</b><span>{pp(isV2?h.active_enrichment:h.sensitivity_enrichment)} относительно общей панели данного вещества.</span></div>
        <div><b>Обогащение зависимости</b><span>{pp(h.dependency_enrichment)} относительно глобальной зависимости от гена.</span></div>
      </div>
      {reasons.length>0&&<div className={styles.callout} style={{marginTop:12}}><b>Что поддерживает гипотезу:</b><ul>{reasons.map((x:string)=><li key={x}>{x}</li>)}</ul></div>}
    </section>

    {isV2&&<section className={styles.section}>
      <div className={styles.sectionHead}><div><div className="eyebrow">МОЛЕКУЛЯРНЫЙ КОНТЕКСТ МИШЕНИ</div><h2>Что известно о мишени в моделях этой опухоли</h2></div></div>
      <div className={styles.definition}>
        <div><b>RNA-экспрессия</b><span>Медиана log₂(TPM+1): {fmt(h.median_expression_log2_tpm1)}. Медианный перцентиль внутри опухолевого контекста: {pct(h.median_expression_percentile_in_cancer)}.</span></div>
        <div><b>Relative copy number</b><span>Медиана: {fmt(h.median_copy_number_relative)}. Перцентиль внутри опухоли: {pct(h.median_copy_number_percentile_in_cancer)}. Это относительный линейный CN; MCL не называет его автоматически амплификацией или делецией.</span></div>
        <div><b>Hotspot-аннотации</b><span>{n(h.hotspot_models_n)} моделей. Hotspot не интерпретируется автоматически как gain-of-function.</span></div>
        <div><b>LikelyLoF</b><span>{n(h.likely_lof_models_n)} моделей. Направление связи с чувствительностью зависит от биологии конкретного гена и пока не добавляет приоритет автоматически.</span></div>
      </div>
      <div className={styles.note}>Молекулярный контекст здесь служит для объяснения и подбора моделей. Для утверждения конкретного биомаркера вида «BRAF V600E → чувствительность к BRAF-ингибитору» будет отдельный направленный слой правил с источником доказательств.</div>
    </section>}

    <ModelTable rows={(roles.positive||[]) as any[]} title="Положительные модели для проверки" isV2={isV2}/>
    <ModelTable rows={(roles.negative_same_cancer||[]) as any[]} title="Отрицательные контроли того же опухолевого контекста" isV2={isV2}/>

    {((roles.discordant_sensitive_without_dependency||[]).length>0||(roles.discordant_dependency_without_activity||[]).length>0||(roles.discordant_dependency_without_sensitivity||[]).length>0)&&<section className={styles.section}>
      <div className={styles.sectionHead}><div><div className="eyebrow">НЕ ИГНОРИРОВАТЬ ПРОТИВОРЕЧИЯ</div><h2>Модели, которые могут опровергнуть простую версию механизма</h2></div></div>
      <div className={styles.callout}>Чувствительность без CRISPR-зависимости может указывать на альтернативную мишень или полифармакологию. CRISPR-зависимость без абсолютного фармакологического эффекта может указывать на недостаточное воздействие, резистентность или несоответствие механизма.</div>
      <ModelTable rows={(roles.discordant_sensitive_without_dependency||[]) as any[]} title="Чувствительны без ожидаемой бинарной CRISPR-зависимости" isV2={isV2}/>
      <ModelTable rows={((roles.discordant_dependency_without_activity||roles.discordant_dependency_without_sensitivity||[]) as any[])} title="CRISPR-зависимы, но выраженного фармакологического ответа нет" isV2={isV2}/>
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
      <div className={styles.callout}><b>Даже Candidate v2 остаётся генератором экспериментальных гипотез, а не доказательством эффективности.</b><ul>{gaps.map((x:string)=><li key={x}>{x}</li>)}</ul></div>
    </section>
  </>;
}
