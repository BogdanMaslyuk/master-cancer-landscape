import Link from "next/link";
import { apiBase, apiGet, formatNumber } from "../../../lib/api";
import styles from "../molecule-analysis.module.css";

export const dynamic = "force-dynamic";

function n(value: unknown, digits = 2) {
  if (value === null || value === undefined || value === "") return "—";
  const x = Number(value);
  return Number.isFinite(x) ? x.toFixed(digits).replace(/\.0+$/, "") : String(value);
}

function roleRu(role: string) {
  if (role === "genetically_dependent_candidate") return "генетическая зависимость от мишени";
  if (role === "genetically_nondependent_candidate") return "генетическая независимость от мишени";
  if (role === "known_drug_context_only") return "данные по известным препаратам, не по PYZ";
  if (role === "positive") return "положительная модель по известному препарату (устаревший слой)";
  if (role?.startsWith("negative_")) return "отрицательный опухолевый контроль";
  if (role?.startsWith("discordant_")) return "противоречивый контекст";
  return role || "дополнительная модель";
}

export default async function MoleculeAnalysisDetail({ params }: { params: Promise<{ id: string }> }) {
  const { id } = await params;
  const data = await apiGet<Record<string, any>>(`/api/molecules/${encodeURIComponent(id)}`);
  const d = data.descriptors || {};
  const targets = (data.targets || []) as any[];
  const cells = (data.cell_line_recommendations || []) as any[];
  const panel = data.experiment_panel || {};
  const assay = panel.suggested_first_assay || {};

  return (
    <>
      <Link href="/molecule-analysis" className={styles.back}>← Анализировать другую молекулу</Link>

      <section className={styles.detailHero}>
        <div className={styles.structure}>
          <img src={`${apiBase}/media/molecules/${encodeURIComponent(id)}/structure.svg`} alt={`2D-структура ${data.name || id}`} />
        </div>
        <div className={styles.summary}>
          <div className="eyebrow">ЗАГРУЖЕННАЯ МОЛЕКУЛА</div>
          <h1>{data.name || id}</h1>
          <div className={styles.id}>{id}{data.matched_pyz_id ? ` · совпадает с ${data.matched_pyz_id}` : ""}</div>
          <code className={styles.smiles}>{data.standardized_smiles}</code>
          <div className={styles.statGrid}>
            <div><b>{n(d.molecular_weight, 1)}</b><span>MW</span></div>
            <div><b>{n(d.clogp)}</b><span>cLogP</span></div>
            <div><b>{n(d.tpsa, 1)}</b><span>TPSA</span></div>
            <div><b>{n(d.qed)}</b><span>QED</span></div>
            <div><b>{d.lipinski_violations ?? "—"}</b><span>Lipinski violations</span></div>
            <div><b>{d.pains_alerts_n ?? "—"} / {d.brenk_alerts_n ?? "—"}</b><span>PAINS / Brenk</span></div>
          </div>
          <div className={styles.guardrail}><b>Интерпретация:</b> {data.guardrail_ru}</div>
        </div>
      </section>

      <section className={styles.panel}>
        <div className={styles.panelHead}>
          <div><div className="eyebrow">TARGET LIGAND SPACE</div><h2>Предполагаемые мишени</h2></div>
          <div className={styles.muted}>Значение — Tanimoto до лучшего экспериментального лиганда, не вероятность связывания.</div>
        </div>
        {targets.length === 0 ? <div className={styles.empty}>Локальный Target Ligand Space не найден или ещё не построен.</div> :
          <div className={styles.targetTableWrap}><table className={styles.table}><thead><tr><th>#</th><th>Мишень</th><th>Tanimoto</th><th>Поддержка</th><th>Лучший лиганд</th><th>Ki/Kd/IC50 известного лиганда</th><th>Соседи ≥0,35</th></tr></thead><tbody>
            {targets.map((t:any)=><tr key={t.target_gene}><td>{t.rank}</td><td><Link href={`/targets/${encodeURIComponent(t.target_gene)}`}><b>{t.target_gene}</b></Link></td><td><b>{n(t.best_tanimoto,3)}</b></td><td>{t.support_label_ru}</td><td><code>{t.best_ligand_id || "—"}</code></td><td>{t.best_lowest_reported_value_nm !== null && t.best_lowest_reported_value_nm !== undefined ? `${formatNumber(t.best_lowest_reported_value_nm,2)} нМ · ${t.best_endpoint || "endpoint"}` : "—"}</td><td>{t.hits_ge_0_35_n ?? 0}</td></tr>)}
          </tbody></table></div>}
      </section>

      <section className={styles.panel}>
        <div className={styles.panelHead}>
          <div><div className="eyebrow">МИШЕНЬ → ОПУХОЛЬ → КЛЕТОЧНАЯ ЛИНИЯ</div><h2>Генетический контекст для проверки гипотезы</h2></div>
          <div className={styles.status}>{cells.filter((x:any)=>x.available_in_laboratory).length} моделей в лабораторном реестре</div>
        </div>
        {cells.length === 0 ? <div className={styles.empty}>Для текущих target-гипотез нет достаточного клеточного контекста. Это не означает отсутствие активности молекулы.</div> :
          <div className={styles.cards}>{cells.slice(0,16).map((c:any)=><article className={styles.cellCard} key={`${c.target_gene}-${c.model_id}`}>
            <div className={styles.cellTop}><div><span className={styles.kicker}>{c.organ || "опухолевый контекст"}</span><h3>{c.cell_line_name}</h3><small>{c.cancer_name || c.cancer_id || "—"}</small></div><span className={`${styles.priority} ${styles[c.priority_level] || ""}`}>{c.priority_level_ru}</span></div>
            <div className={styles.badges}><span>{c.target_gene} · Tanimoto {n(c.target_similarity,3)}</span><span>{roleRu(c.role)}</span>{c.available_in_laboratory && <span className={styles.lab}>зарегистрирована в лаборатории; наличие требует проверки</span>}</div>
            <div className={styles.miniStats}><div><b>{n(c.dependency_probability)}</b><span>P(dep)</span></div><div><b>{n(c.gene_effect)}</b><span>Gene Effect</span></div></div>
            {Array.isArray(c.reference_drug_observations) && c.reference_drug_observations.length > 0 && (
              <div className={styles.guardrail}>
                <b>Ответы на известные препараты (не на исследуемую PYZ):</b>
                <div>Всего измерений: {c.reference_drug_observations_n ?? c.reference_drug_observations.length}.
                  На карточке показаны первые восемь наблюдений.</div>
                <ul>{c.reference_drug_observations.slice(0,8).map((o:any, i:number)=>
                  <li key={`${o.observation_id || o.hypothesis_id || o.compound_id}-${i}`}>
                    {o.preferred_name || o.compound_id || "Вещество без названия"}:
                    {" "}LFC {n(o.prism_lfc,3)}
                    {o.reference_model_role ? ` · роль Candidate v2: ${o.reference_model_role}` : ""}
                    {o.candidate_priority_status ? ` · статус гипотезы: ${o.candidate_priority_status}` : ""}
                    {o.observation_id ? ` · наблюдение: ${o.observation_id}` : ""}
                    {o.source_release ? ` · PRISM: ${o.source_release}` : ""}
                  </li>
                )}</ul>
                {c.source_layer === "full_processed_depmap_target_context"
                  ? "Нормализованные исходные записи PRISM, связанные с известными веществами и моделями. Эффект через CDK4 и чувствительность к PYZ не доказаны."
                  : "Исторические записи ограниченного набора Candidate v2; полный массив PRISM здесь не представлен."}
              </div>
            )}
            <ul>{(c.reasons || []).slice(0,5).map((reason:string)=><li key={reason}>{reason}</li>)}</ul>
            {Array.isArray(c.warnings) && c.warnings.length > 0 && (
              <div className={styles.guardrail}>{c.warnings.join(" ")}</div>
            )}
            {c.model_id && <Link className={styles.more} href={`/models/${encodeURIComponent(c.model_id)}`}>Открыть модель →</Link>}
          </article>)}</div>}
      </section>

      <section className={styles.panel}>
        <div className={styles.panelHead}><div><div className="eyebrow">РЕКОМЕНДУЕМЫЙ ПЕРВЫЙ ЭКСПЕРИМЕНТ</div><h2>Минимальная информативная панель</h2></div><div className={styles.status}>мишень · {panel.focus_target_gene || "—"}</div></div>
        {panel.status !== "complete" && (
          <div className={styles.guardrail}>
            <b>{panel.status === "partial" ? "Экспериментальная панель неполная." : "Панель не готова."}</b>
            {" "}{panel.interpretation_ru || "Для формирования панели недостаточно подтверждённых данных."}
            {Array.isArray(panel.missing_components) && panel.missing_components.length > 0 && (
              <div>Отсутствующие элементы: {panel.missing_components.map((code: string) => ({
                supported_target_hypothesis: "мишень с достаточной структурной поддержкой",
                positive_tumour_model: "положительная опухолевая модель",
                independent_confirmatory_model: "вторая независимая модель",
                measured_negative_tumour_control: "отрицательный опухолевый контроль с измеренной CRISPR-зависимостью",
                non_tumour_control: "неопухолевый контроль",
                laboratory_positive_tumour_model: "положительная модель в нашей лаборатории",
                laboratory_confirmatory_model: "вторая положительная модель в нашей лаборатории",
                laboratory_negative_tumour_control: "отрицательный опухолевый контроль в нашей лаборатории",
                same_cancer_negative_control: "сопоставимый отрицательный контроль того же типа опухоли",
                experimentally_supported_uploaded_molecule_target: "подтверждение взаимодействия исследуемой молекулы с мишенью",
                compound_specific_pharmacology: "фармакологические измерения для исследуемой молекулы",
                validated_assay_protocol: "проверенный экспериментальный протокол",
                verified_laboratory_availability: "подтверждённая готовность культур в лаборатории",
                full_depmapped_evidence_layer: "полный проверенный слой данных DepMap (вместо усечённого Candidate v2)",
              } as Record<string, string>)[code] || code).join("; ")}.</div>
            )}
          </div>
        )}
        {!panel.available ? <div className={styles.empty}>Для этой структуры пока нельзя сформировать обоснованную мишень-ориентированную панель.</div> : <>
          {panel.panel_type === "exploratory_genetic_context" && (
            <div className={styles.guardrail}>
              <b>Исследовательский генетический контекст (не прогноз активности PYZ).</b>
              {" "}CRISPR-зависимые модели: {(panel.genetically_dependent_candidates || []).slice(0,6).map((x:any)=>x.cell_line_name).join(" · ") || "не найдены"}.
              {" "}CRISPR-независимые кандидаты для сравнения: {(panel.genetically_nondependent_candidates || []).slice(0,6).map((x:any)=>x.cell_line_name).join(" · ") || "не найдены"}.
              {" "}Эти группы не являются проверенными положительным и отрицательным лекарственными контролями.
            </div>
          )}
          <div className={styles.experimentGrid}>
            <div><span>Положительные модели</span><b>{(panel.positive_models || []).map((x:any)=>x.cell_line_name).join(" · ") || "не найдены"}</b></div>
            <div><span>Отрицательный опухолевый контроль</span><b>{panel.negative_tumor_control?.cell_line_name || "не найден"}</b></div>
            <div><span>Неопухолевый контроль</span><b>{panel.non_tumor_control?.lab_name || "не найден"}</b></div>
            <div><span>Первичный assay</span><b>{assay.assay || "—"} · {assay.exposure_hours || "—"} ч</b></div>
          </div>
          <div className={styles.guardrail}>{assay.note_ru}</div>
        </>}
      </section>

      <section className={styles.panel}>
        <div className={styles.panelHead}><div><div className="eyebrow">NORMAL TISSUE LIABILITY</div><h2>Что ещё не закрыто</h2></div></div>
        <div className={styles.guardrail}>{data.normal_tissue_liability?.interpretation_ru}</div>
      </section>
    </>
  );
}
