import Link from "next/link";
import { apiGet } from "../../lib/api";
import styles from "../pharmacology.module.css";

export const dynamic = "force-dynamic";

type Params=Record<string,string|string[]|undefined>;
type Payload=Record<string,any>;
function one(value:string|string[]|undefined){return Array.isArray(value)?value[0]:value;}
function n(value:any){const x=Number(value);return Number.isFinite(x)?x.toLocaleString("ru-RU"):"—";}
function fmt(value:any,digits=2){if(value===null||value===undefined||value==="")return "—";const x=Number(value);return Number.isFinite(x)?x.toFixed(digits):String(value);}
const readinessRu:Record<string,string>={
  ready_with_internal_tumor_control:"Готово: есть положительная модель и опухолевый контроль",
  ready_positive_model:"Есть положительная модель",
  discordant_current_panel:"Панель даёт противоречивый сигнал",
  no_support_in_current_panel:"На текущей панели поддержки нет",
};
const roleRu:Record<string,string>={
  human_tumor:"опухолевая",
  human_non_tumor_control:"неопухолевый контроль",
  human_technical:"техническая",
  other:"прочая",
};
const matchRu:Record<string,string>={
  matched:"сопоставлена",
  not_found:"не найдена в DepMap",
  requires_review_ambiguous:"нужно проверить идентичность",
  not_applicable_nonhuman:"не относится к human DepMap",
};

export default async function LaboratoryPage({searchParams}:{searchParams:Promise<Params>}){
  const sp=await searchParams;
  const q=one(sp.q)||"";
  const readiness=one(sp.readiness)||"";
  const candidateQuery=new URLSearchParams({limit:"150",offset:"0"});
  if(q)candidateQuery.set("q",q);
  if(readiness)candidateQuery.set("readiness",readiness);
  const [summary,lines,candidates]=await Promise.all([
    apiGet<Payload>("/api/laboratory/summary"),
    apiGet<Payload>("/api/laboratory/lines?species=human"),
    apiGet<Payload>(`/api/laboratory/candidates?${candidateQuery.toString()}`),
  ]);
  const panel=summary.panel||{};
  const triage=summary.candidate_triage||{};
  const counts=triage.readiness_counts||{};
  const lineItems=(lines.items||[]) as any[];
  const candidateItems=(candidates.items||[]) as any[];

  return <>
    <section className={styles.pageHeader}>
      <div>
        <div className="eyebrow">ФИЗИЧЕСКИ ДОСТУПНЫЕ КЛЕТОЧНЫЕ МОДЕЛИ</div>
        <h1>Лабораторная панель кафедры</h1>
        <p>MCL теперь отделяет вычислительно интересные гипотезы от реально выполнимых экспериментов. Эта страница показывает только те клеточные линии, которые есть в лаборатории, их сопоставление с DepMap и пересечение с 222 приоритетами Candidate v2.</p>
      </div>
      <div className={styles.heroStats}>
        <div className={styles.heroStat}><strong>{n(panel.human_tumor_lines_n)}</strong><span>человеческих опухолевых линий</span></div>
        <div className={styles.heroStat}><strong>{n(counts.ready_with_internal_tumor_control||0)}</strong><span>гипотез с внутренним опухолевым контролем</span></div>
        <div className={styles.heroStat}><strong>{n(counts.ready_positive_model||0)}</strong><span>ещё гипотез с положительной моделью</span></div>
      </div>
    </section>

    {!summary.available?<div className={styles.empty}><b>Лабораторная панель ещё не построена.</b><br/><code>{summary.build_command||".\\scripts\\build-laboratory-panel.ps1"}</code></div>:<>
      <section className={styles.section}>
        <div className={styles.sectionHead}><div><div className="eyebrow">КРИТИЧЕСКОЕ ОГРАНИЧЕНИЕ</div><h2>BJ5ta — наш доступный человеческий неопухолевый контроль</h2></div></div>
        <div className={styles.callout}><b>Как её трактовать:</b> {panel.control_interpretation_ru||summary.interpretation_ru} Слабый эффект на BJ5ta при сильном эффекте на опухолевой линии будет полезным указанием на лабораторную селективность, но не доказательством безопасности для нормальной ткани соответствующего органа.</div>
      </section>

      <section className={styles.section}>
        <div className={styles.sectionHead}><div><div className="eyebrow">РЕЕСТР</div><h2>Человеческие линии, которые реально есть на кафедре</h2></div></div>
        <div className={styles.note}>{n(panel.laboratory_lines_n)} линий всего · {n(panel.human_lines_n)} человеческих · {n(panel.human_matched_depmap_n)} человеческих сопоставлены с DepMap · {n(panel.human_in_crispr_atlas_n)} входят в наш CRISPR Cancer Atlas. Неоднозначное сопоставление никогда не разрешается автоматически.</div>
        <div className={styles.tableWrap}><table className={styles.table}>
          <thead><tr><th>Линия</th><th>Роль в лаборатории</th><th>Контекст</th><th>DepMap</th><th>Покрытие MCL</th></tr></thead>
          <tbody>{lineItems.map((row:any)=><tr key={row.lab_id}>
            <td><b>{row.lab_name}</b><span className={styles.sub}>{row.highlighted_on_source?"выделена в исходном списке как приоритетная":"есть в коллекции кафедры"}</span></td>
            <td><b>{roleRu[row.laboratory_role]||row.laboratory_role}</b>{row.control_role&&<span className={styles.sub}>единственный общий человеческий неопухолевый контроль</span>}</td>
            <td>{row.mcl_cancer_name||row.disease_group_ru||"—"}<span className={styles.sub}>{row.mcl_organ_ru||row.note_ru||""}</span></td>
            <td>{row.model_id?<Link className={styles.primary} href={`/models/${encodeURIComponent(row.model_id)}`}>{row.model_id}</Link>:<b>—</b>}<span className={styles.sub}>{matchRu[row.match_status]||row.match_status}</span>{row.depmap_cell_line_name&&<span className={styles.sub}>{row.depmap_cell_line_name}</span>}</td>
            <td><b>{row.has_crispr_atlas?"CRISPR ✓":"CRISPR-атлас —"}</b><span className={styles.sub}>{row.model_class==="technical_transformed"?"не использовать как нормальную ткань":row.control_role?"контрольная роль задана отдельно":""}</span></td>
          </tr>)}</tbody>
        </table></div>
      </section>

      <section className={styles.section}>
        <div className={styles.sectionHead}><div><div className="eyebrow">222 → РЕАЛЬНЫЙ ЭКСПЕРИМЕНТ</div><h2>Кандидаты, проверяемые на нашей панели</h2></div></div>
        <form action="/laboratory" method="get" className={styles.toolbar}>
          <label className={styles.field}><span>Поиск</span><input name="q" defaultValue={q} placeholder="вещество, BRAF, breast…"/></label>
          <label className={styles.field}><span>Лабораторная готовность</span><select name="readiness" defaultValue={readiness}><option value="">все 222 приоритета</option><option value="ready_with_internal_tumor_control">есть положительная + отрицательная опухолевая модель</option><option value="ready_positive_model">есть положительная модель</option><option value="discordant_current_panel">противоречивый сигнал</option><option value="no_support_in_current_panel">нет поддержки на панели</option></select></label>
          <button className={styles.button} type="submit">Применить</button><Link className={styles.reset} href="/laboratory">Сбросить</Link>
        </form>
        <div className={styles.note}>{candidates.interpretation_ru||"Лабораторная готовность не является новым биологическим доказательством — это фильтр выполнимости."}</div>
        {!candidates.available?<div className={styles.empty}><b>Пересечение Candidate v2 с лабораторной панелью ещё не построено.</b><br/><code>{candidates.build_command}</code></div>:
        <div className={styles.tableWrap}><table className={styles.table}>
          <thead><tr><th>Кандидат</th><th>Мишень</th><th>Опухоль</th><th>Положительные линии кафедры</th><th>Отрицательные опухолевые контроли</th><th>BJ5ta</th><th>Готовность</th></tr></thead>
          <tbody>{candidateItems.map((row:any)=><tr key={row.hypothesis_id}>
            <td><Link className={styles.primary} href={`/hypotheses/${encodeURIComponent(row.hypothesis_id)}`}>{row.preferred_name||row.compound_id}</Link><span className={styles.sub}>{row.compound_id}</span></td>
            <td><Link className={styles.primary} href={`/targets/${encodeURIComponent(row.target_gene)}`}>{row.protein_preferred_name||row.target_gene}</Link><span className={styles.sub}>ген {row.target_gene}</span></td>
            <td><b>{row.mcl_cancer_name||"—"}</b><span className={styles.sub}>{row.mcl_organ_ru||""}</span></td>
            <td><b>{n(row.laboratory_positive_models_n)}</b><span className={styles.sub}>{row.positive_lab_lines||"нет"}</span></td>
            <td><b>{n(row.laboratory_negative_models_n)}</b><span className={styles.sub}>{row.negative_lab_lines||"нет"}</span></td>
            <td>{row.bj5ta_prism_available?<><b>LFC {fmt(row.bj5ta_prism_lfc)}</b><span className={styles.sub}>из имеющегося PRISM-слоя</span></>:<><b>доступна физически</b><span className={styles.sub}>ответ на это вещество нужно измерить</span></>}</td>
            <td><Link className={styles.chip} href={`/hypotheses/${encodeURIComponent(row.hypothesis_id)}`}>{readinessRu[row.laboratory_readiness]||row.laboratory_readiness}</Link></td>
          </tr>)}</tbody>
        </table></div>}
      </section>
    </>}
  </>;
}
