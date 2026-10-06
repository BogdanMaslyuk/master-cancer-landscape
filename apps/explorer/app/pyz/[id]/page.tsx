import Link from "next/link";
import { apiBase, apiGet, formatNumber } from "../../../lib/api";
import styles from "../pyz.module.css";

export const dynamic = "force-dynamic";

type Params=Record<string,string|string[]|undefined>;
function one(v:string|string[]|undefined){return Array.isArray(v)?v[0]:v;}
function f(v:any,d=3){const x=Number(v);return Number.isFinite(x)?x.toFixed(d):"—";}
function roleRu(value:any){return ({positive:"приоритетная для проверки",negative_same_cancer:"отрицательный контроль той же опухоли",negative_panel_comparator:"отрицательный контроль",discordant_sensitive_without_dependency:"чувствительная без CRISPR-зависимости",discordant_dependency_without_activity:"CRISPR-зависима без фармакологического ответа",unclassified:"дополнительная модель"} as Record<string,string>)[String(value)]||String(value||"дополнительная модель");}
function transferRu(value:any){return ({exploratory_low_structural_similarity:"низкое, но заметное сходство",mechanism_not_transferable_by_2d_similarity:"механизм не переносим по 2D",moderate_structural_transferability:"умеренное сходство",high_structural_transferability:"высокое сходство"} as Record<string,string>)[String(value)]||String(value||"—");}
function riskRu(value:any){return value==="high_predicted_risk"?"высокий предиктивный риск":"требует экспериментальной проверки";}

export default async function PyzDetailPage({params,searchParams}:{params:Promise<{id:string}>,searchParams:Promise<Params>}){
  const {id}=await params;
  const sp=await searchParams;
  const focusTarget=one(sp.target_gene)||"";
  const suffix=focusTarget?`?target_gene=${encodeURIComponent(focusTarget)}`:"";
  const data=await apiGet<Record<string,any>>(`/api/pyz/${encodeURIComponent(id)}${suffix}`);
  const c=data.compound||{};
  const evidence=(data.target_evidence||[]) as any[];
  const contexts=(data.relevant_contexts||[]) as any[];
  const priorityTargets=(data.priority_target_genes||[]) as string[];
  const admet=data.admet_summary||{};

  return <>
    <Link href="/pyz" className={styles.backLink}>← Все молекулы PYZ</Link>
    <section className={styles.detailHero}>
      <div className={styles.structureCard}><img src={`${apiBase}/api/pyz/${encodeURIComponent(c.own_compound_id)}/structure.svg`} alt={`2D-структура ${c.own_compound_id}`}/></div>
      <div className={styles.detailCard}>
        <div className={styles.detailTitle}><div><div className="eyebrow">СОБСТВЕННАЯ МОЛЕКУЛА</div><h1>{c.own_compound_id}</h1><div className={styles.sectionCopy}>{c.chemotype_ru} · {c.ring_variant} · тир {c.tier}</div></div><span className={c.integrated_risk==="high_predicted_risk"?styles.badgeRisk:styles.badgeSoft}>{riskRu(c.integrated_risk)}</span></div>
        <div className={styles.badges}>{priorityTargets.length?priorityTargets.map(g=><span key={g} className={styles.badgePriority}>{g} · проверять</span>):<span className={styles.badge}>нет пар ≥0,35</span>}<span className={styles.badge}>докинг: ожидается</span></div>
        <div className={styles.smiles}>{c.standardized_smiles}</div>
        <div className={styles.facts}>
          <div className={styles.fact}><b>{c.qed!==null&&c.qed!==undefined?Number(c.qed).toFixed(2):"—"}</b><span>QED · характеристика drug-likeness, не безопасность</span></div>
          <div className={styles.fact}><b>{priorityTargets.length}</b><span>мишеней выше исследовательского порога 0,35</span></div>
          <div className={styles.fact}><b>{c.series_range||"—"}</b><span>структурная серия</span></div>
          <div className={styles.fact}><b>не измерено</b><span>прямое связывание PYZ–мишень</span></div>
        </div>
        <div className={styles.warning} style={{marginTop:14}}><b>ADMET:</b> {c.admet_note_ru||admet.interpretation_ru} Ни один предиктивный отрицательный класс не заменяет экспериментальный контроль.</div>
      </div>
    </section>

    <section className={styles.section}>
      <div className={styles.sectionHead}><div><div className="eyebrow">TARGET LIGAND SPACE · 13 CORE-МИШЕНЕЙ</div><h2>Какие мишени сейчас поддерживает химическое пространство</h2><div className={styles.sectionCopy}>Для каждой мишени показан лучший экспериментально измеренный лиганд. Ki, Kd и IC50 не считаются взаимозаменяемыми; значение относится к известному лиганду, а не к {c.own_compound_id}.</div></div><Link href={`/pyz/matrix?q=${encodeURIComponent(c.own_compound_id)}`} className={styles.secondaryButton}>Показать в матрице →</Link></div>
      <div className={styles.targetGrid}>{evidence.map((row:any,idx:number)=>{
        const priority=Number(row.best_tanimoto_morgan_r2_2048)>=0.35;
        const focused=focusTarget&&String(row.target_gene).toUpperCase()===focusTarget.toUpperCase();
        return <article className={`${styles.targetCard} ${focused?styles.focus:""}`} key={row.target_gene}>
          <div className={styles.targetHead}><div><h3>{row.target_gene}</h3><div className={styles.sectionCopy}>{transferRu(row.structural_transfer_class)}</div></div><div className={styles.sim}>{f(row.best_tanimoto_morgan_r2_2048)}</div></div>
          <div className={styles.badges}>{priority?<span className={styles.badgePriority}>ортогональная проверка оправдана</span>:<span className={styles.badgeSoft}>не переносить механизм</span>}{row.best_structural_ligand_cellular_label&&<span className={styles.badge}>{row.best_structural_ligand_cellular_label}</span>}</div>
          <div className={styles.targetMeta}>
            <div><b>{formatNumber(row.ligands_compared_n,0)}</b><span>экспериментальных лигандов</span></div>
            <div><b>{row.best_ligand_id||"—"}</b><span>лучший структурный сосед</span></div>
            <div><b>{row.best_lowest_reported_value_nm_across_endpoints!==null&&row.best_lowest_reported_value_nm_across_endpoints!==undefined?`${formatNumber(row.best_lowest_reported_value_nm_across_endpoints,2)} нМ`:"—"}</b><span>{row.best_lowest_value_endpoint_type||"прямой endpoint неизвестен"}</span></div>
          </div>
          <div className={styles.targetMeta}>
            <div><b>{formatNumber(row.hits_ge_0_35_n,0)}</b><span>соседей ≥0,35</span></div>
            <div><b>{row.strict_oncology_best_reference_name||"—"}</b><span>лучший строгий онкоэталон</span></div>
            <div><b>{row.best_cellular_ligand_id||"—"}</b><span>лучший сосед с PRISM/CRISPR</span></div>
          </div>
          {(row.top_neighbors||[]).length>0&&<div className={styles.neighborList}>{(row.top_neighbors||[]).map((n:any)=><div className={styles.neighbor} key={`${row.target_gene}-${n.rank_in_target}-${n.ligand_id}`}><span>#{n.rank_in_target} · {f(n.tanimoto_morgan_r2_2048)}</span><code>{n.ligand_id}</code><span>{n.lowest_reported_value_nm_across_endpoints!==null&&n.lowest_reported_value_nm_across_endpoints!==undefined?`${formatNumber(n.lowest_reported_value_nm_across_endpoints,2)} нМ ${n.lowest_value_endpoint_type||""}`:"—"}</span></div>)}</div>}
          <div className={styles.splitActions} style={{marginTop:12}}><Link href={`/pyz/${encodeURIComponent(c.own_compound_id)}?target_gene=${encodeURIComponent(row.target_gene)}`} className={styles.secondaryButton}>Органы и культуры →</Link><Link href={`/targets/${encodeURIComponent(row.target_gene)}`} className={styles.reset}>Карточка белка</Link></div>
        </article>;
      })}</div>
    </section>

    <section className={styles.section}>
      <div className={styles.sectionHead}><div><div className="eyebrow">PYZ → МИШЕНЬ → ОПУХОЛЬ → КЛЕТОЧНАЯ ЛИНИЯ</div><h2>Где рационально проверять {c.own_compound_id}</h2><div className={styles.sectionCopy}>{data.context_guardrail_ru}</div></div></div>
      {contexts.length===0?<div className={styles.empty}>Для этой молекулы пока нет target-based маршрута с достаточной опорой. Это не означает отсутствие активности: нужен независимый поиск мишени или результаты докинга.</div>:
      <div className={styles.contextGrid}>{contexts.map((ctx:any,index:number)=><article className={styles.contextCard} key={`${ctx.target_gene}-${ctx.mcl_cancer_id||ctx.mcl_cancer_name}-${index}`}>
        <div className={styles.contextHeader}><div><div className={styles.organ}>{ctx.mcl_organ_ru||ctx.mcl_system_ru||"опухолевый контекст"}</div><h3>{ctx.mcl_cancer_name||ctx.mcl_cancer_id||"контекст Candidate v2"}</h3><div className={styles.sectionCopy}>мишень {ctx.target_gene} · Tanimoto {f(ctx.pyz_target_similarity)} · {ctx.interpretation_ru}</div></div><span className={ctx.molecule_relevance_status==="priority_target_based_test_context"?styles.badgePriority:styles.badgeSoft}>{ctx.molecule_relevance_status==="priority_target_based_test_context"?"приоритет проверки":"контекст мишени"}</span></div>
        <div className={styles.contextStats}><span className={styles.badge}>{formatNumber(ctx.positive_models_n,0)} положительных моделей</span><span className={styles.badge}>{formatNumber(ctx.negative_models_n,0)} отрицательных контролей</span><span className={ctx.laboratory_models_n?styles.badgePriority:styles.badgeSoft}>{formatNumber(ctx.laboratory_models_n,0)} есть у нас</span><span className={styles.badge}>{formatNumber(ctx.priority_hypotheses_n,0)} Candidate v2 гипотез</span></div>
        {(ctx.models||[]).length>0&&<div className={styles.cellTableWrap}><table className={styles.cellTable}><thead><tr><th>Клеточная модель</th><th>Роль</th><th>CRISPR P(dep)</th><th>Chronos Gene Effect</th><th>Известные вещества</th><th>Наша лаборатория</th></tr></thead><tbody>{(ctx.models||[]).map((m:any)=><tr key={`${ctx.target_gene}-${m.model_id}`}><td><div className={styles.cellName}>{m.cell_line_name||m.model_id}</div>{m.model_id&&<Link href={`/models/${encodeURIComponent(m.model_id)}`} className={styles.muted}>{m.model_id}</Link>}</td><td>{roleRu(m.role)}</td><td>{f(m.dependency_probability,2)}</td><td>{f(m.gene_effect,2)}</td><td>{(m.known_reference_compounds||[]).slice(0,4).join(" · ")||"—"}</td><td>{m.available_in_laboratory?<span className={styles.labYes}>есть · {(m.laboratory_lines||[]).join(" / ")}</span>:<span className={styles.labNo}>нет в текущем реестре</span>}</td></tr>)}</tbody></table></div>}
      </article>)}</div>}
    </section>

    <section className={styles.section}>
      <div className={styles.notice}><b>Следующая ось доказательств — докинг.</b> Сейчас статус намеренно указан как «ожидается». После получения результатов сюда должны добавиться PDB, redocking-QC, положение в кармане, контакты и ранжирование внутри одной структуры рецептора. Docking score не будет трактоваться как Ki/Kd/IC50.</div>
    </section>
  </>;
}
