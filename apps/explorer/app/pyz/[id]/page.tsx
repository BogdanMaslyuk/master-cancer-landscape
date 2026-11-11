import Link from "next/link";
import { notFound } from "next/navigation";
import { apiBase, apiGet } from "../../../lib/api";
import type { EvidenceClaim, PyzDetailV2Response, PyzTarget, PyzContext } from "../../../lib/pyz-profile";
import { CopySmiles, StructurePreview } from "../../../components/pyz/PyzInteractions";
import styles from "./profile-v2.module.css";

export const dynamic = "force-dynamic";

function num(value: number | null | undefined, digits = 3): string {
  return value === null || value === undefined || !Number.isFinite(value)
    ? "нет данных"
    : new Intl.NumberFormat("ru-RU", { maximumFractionDigits: digits }).format(value);
}

function Claim({ value }: { value: EvidenceClaim }) {
  return <li className={styles.claim}><strong>{value.title}</strong>
    <span>{value.meaning}</span><small>Основание: {value.source} · {value.rule_id}</small></li>;
}

function TargetCard({ item, id, focus }: { item: PyzTarget; id: string; focus: boolean }) {
  return <article className={focus ? styles.targetFocus : styles.targetCard} id={"target-" + item.target_gene}>
    <div className={styles.targetHeader}>
      <div><Link className={styles.targetLink} href={"/targets/" + encodeURIComponent(item.target_gene)}>{item.target_gene}</Link>
        <small>{item.uniprot_accession ? "UniProt " + item.uniprot_accession : "UniProt не указан"}</small></div>
      <div className={styles.similarity}><b>{num(item.tanimoto)}</b><span>Танимото</span></div>
      <div className={item.passes_navigation_cutoff ? styles.pillBlue : styles.pillNeutral}>
        {item.passes_navigation_cutoff ? "Структурная гипотеза" : "Недостаточно сходства"}
      </div>
    </div>
    <details open={focus || undefined} className={styles.details}>
      <summary>Сведения об аналогах и проверка мишени</summary>
      <div className={styles.detailsBody}>
        <p className={styles.subtle}>Сравнивалось лигандов: {num(item.ligands_compared_n, 0)}.
          Сходство не доказывает взаимодействия PYZ с этим белком.</p>
        <div className={styles.referenceGrid}>
          <div><span>Известный лиганд</span><b>{item.best_ligand_id || "не указан"}</b></div>
          <div><span>Тип измерения лиганда</span><b>{item.reference_endpoint || "не указан"}</b></div>
          <div><span>Справочное значение лиганда</span><b>{item.reference_value_nm !== null ? num(item.reference_value_nm, 2) + " нМ" : "нет данных"}</b></div>
        </div>
        <p className={styles.caution}>Это значение относится к известному веществу, а не PYZ.
          Оно получено из прежнего агрегата измерений; Ki, Kd и IC50 нельзя считать взаимозаменяемыми
          или использовать для сравнения активности разных мишеней.</p>
        {item.neighbors.length > 0 && <details className={styles.innerDetails}>
          <summary>Показать ближайшие структуры ({item.neighbors.length})</summary>
          <ul className={styles.neighbors}>{item.neighbors.map((x, i) =>
            <li key={String(x.ligand_id || i)}><span>{x.ligand_id || "без ID"}</span><span>Танимото {num(x.tanimoto)}</span></li>
          )}</ul>
        </details>}
        <div className={styles.actions}>
          <Link className={styles.primaryAction} href={"/pyz/" + encodeURIComponent(id) + "?target_gene=" + encodeURIComponent(item.target_gene) + "#models"}>Опухоли и культуры</Link>
          <Link className={styles.secondaryAction} href={"/targets/" + encodeURIComponent(item.target_gene)}>Карточка белка</Link>
        </div>
      </div>
    </details>
  </article>;
}

const roleNames: Record<string, string> = {
  positive: "Модель для проверки гипотезы о мишени",
  negative_same_cancer: "Отрицательный опухолевый контроль",
  negative_panel_comparator: "Сравнительная модель",
  discordant_sensitive_without_dependency: "Известный препарат активен без CRISPR-зависимости",
  discordant_dependency_without_activity: "CRISPR-зависимость без эффекта известного препарата",
};

function ContextCard({ item }: { item: PyzContext }) {
  return <article className={styles.contextCard}>
    <header className={styles.contextHeader}><div><span className={styles.eyebrow}>{item.mcl_organ_ru || "Орган не указан"}</span>
      <h3>{item.mcl_cancer_name || item.mcl_cancer_id || "Контекст мишени"}</h3>
      <p>{item.target_gene || "Белок не указан"} · Танимото {num(item.tanimoto)}</p></div>
      <span className={styles.pillNeutral}>Не измерение PYZ</span></header>
    {item.models.length ? <div className={styles.tableScroll}><table className={styles.table}>
      <thead><tr><th>Клеточная линия</th><th>Роль при проверке</th><th>CRISPR P(dep)</th><th>Gene Effect</th><th>Наличие</th></tr></thead>
      <tbody>{item.models.map((model, i) => <tr key={String(model.model_id || i)}>
        <td>{model.model_id ? <Link className={styles.tableLink} href={"/models/" + encodeURIComponent(model.model_id)}>{model.cell_line_name || model.model_id}</Link> : model.cell_line_name || "нет названия"}</td>
        <td>{roleNames[model.role || ""] || model.role || "не определена"}</td>
        <td>{num(model.dependency_probability, 2)}</td>
        <td>{num(model.gene_effect, 2)}</td>
        <td>{model.available_in_laboratory_registry ? "В реестре; проверить наличие" : "Не отмечена в реестре"}</td>
      </tr>)}</tbody>
    </table></div> : <p className={styles.empty}>Сопоставленных клеточных моделей нет.</p>}
  </article>;
}

export default async function PyzDetailPage({ params, searchParams }: {
  params: Promise<{ id: string }>;
  searchParams: Promise<Record<string, string | string[] | undefined>>;
}) {
  const { id } = await params;
  const query = await searchParams;
  const input = query.target_gene;
  const focus = (Array.isArray(input) ? input[0] : input) || "";
  let data: PyzDetailV2Response;
  try {
    data = await apiGet<PyzDetailV2Response>(
      "/api/pyz/" + encodeURIComponent(id) + (focus ? "?target_gene=" + encodeURIComponent(focus) : "")
    );
  } catch (error) {
    if (error instanceof Error && error.message.startsWith("404:")) notFound();
    throw error;
  }
  const profile = data.profile_v2;
  if (!profile || profile.schema_version !== "mcl-pyz-profile-v2-p0") {
    throw new Error("Для карточки PYZ 2.0 требуется обновить локальный сервер API.");
  }
  const summary = profile.scientific_summary;
  const targets = profile.target_hypotheses.items;
  const contexts = profile.cell_context.items;
  const smiles = profile.identity.standardized_smiles;
  const cutoffHits = targets.filter(x => x.passes_navigation_cutoff).length;

  return <main className={styles.page}>
    <div className={styles.breadcrumbs}>
      <Link href="/pyz">Молекулы PYZ</Link><span>/</span><strong>{profile.compound_id}</strong>
      <Link className={styles.matrixLink} href={"/pyz/matrix?q=" + encodeURIComponent(profile.compound_id)}>Открыть матрицу</Link>
    </div>
    <div className={styles.pageHeader}><div>
      <span className={styles.eyebrow}>MCL · НАУЧНАЯ КАРТОЧКА 2.0</span>
      <h1>{profile.compound_id}</h1>
      <p>{profile.identity.chemotype_ru || "Не указан химический тип"} · {profile.identity.ring_variant || "вариант не указан"} · {profile.identity.series_range || "серия не указана"}</p>
    </div><span className={styles.pillBlue}>Исследовательская стадия</span></div>

    <nav className={styles.jumpNav} aria-label="Разделы карточки">
      <a href="#conclusion">Заключение</a><a href="#structure">Структура</a>
      <a href="#admet">ADMET</a><a href="#targets">Мишени</a>
      <a href="#models">Клеточные линии</a><a href="#sources">Источники</a>
    </nav>

    <section className={styles.heroGrid} id="structure" aria-label="Идентичность молекулы">
      <div className={styles.structure}><StructurePreview compoundId={profile.compound_id}
        src={apiBase + "/media/pyz/" + encodeURIComponent(profile.compound_id) + "/structure.svg"} />
        <p>2D-структура из реестра PYZ</p></div>
      <div className={styles.identity}>
        <span className={styles.eyebrow}>ХИМИЧЕСКАЯ ИДЕНТИЧНОСТЬ</span>
        <div className={styles.facts}>
          <div><span>QED <abbr title="Описательная оценка drug-likeness, не доказательство эффективности или безопасности">ⓘ</abbr></span><strong>{num(profile.properties.qed, 2)}</strong></div>
          <div><span>Серия</span><strong>{profile.identity.series_range || "нет данных"}</strong></div>
          <div><span>Мишеней с Танимото ≥ 0,35</span><strong>{cutoffHits}</strong></div>
          <div><span>Прямые опыты PYZ–белок</span><strong>Не измерено</strong></div>
        </div>
        <h2>Стандартизованный SMILES</h2>
        <code className={styles.smiles}>{smiles || "нет данных"}</code>
        <CopySmiles compoundId={profile.compound_id} smiles={smiles} />
      </div>
    </section>

    <section className={styles.section} id="conclusion">
      <div className={styles.sectionHeading}><div><span className={styles.eyebrow}>ГЛАВНОЕ ДЛЯ ИССЛЕДОВАТЕЛЯ</span>
        <h2>Научное заключение</h2></div>
        <span className={styles.pillNeutral}>Проверяемые правила, не прогноз эффективности</span></div>
      <div className={styles.summaryGrid}>
        <article className={styles.summaryCard}><span className={styles.eyebrow}>01 · ФАКТЫ</span><h3>Что установлено</h3>
          <ul>{summary.established.map(x => <Claim key={x.rule_id} value={x} />)}</ul></article>
        <article className={styles.summaryCard}><span className={styles.eyebrow}>02 · ГИПОТЕЗЫ</span><h3>Что предполагается</h3>
          {summary.hypotheses.length ? <ul>{summary.hypotheses.map(x => <Claim key={x.rule_id} value={x} />)}</ul> :
            <p className={styles.empty}>Нет структурной опоры выше порога. Это не означает отсутствия активности.</p>}</article>
        <article className={styles.summaryCard}><span className={styles.eyebrow}>03 · НЕИЗВЕСТНО</span><h3>Что требует проверки</h3>
          <ul>{summary.gaps.map(x => <Claim key={x.rule_id} value={x} />)}</ul></article>
      </div>
      <div className={styles.nextStep}><div><span className={styles.eyebrow}>04 · СЛЕДУЮЩИЙ ШАГ</span>
        <h3>{summary.next_step.title}</h3><p>{summary.next_step.meaning}</p>
        <small>Правило {summary.next_step.rule_id} · {summary.next_step.source}</small></div>
        <Link className={styles.primaryAction} href="/laboratory">Лабораторные модели</Link></div>
    </section>

    <div className={styles.metrics}>
      <div><b>{targets.length}</b><span>мишеней текущей библиотеки</span></div>
      <div><b>{cutoffHits}</b><span>исследовательских пар ≥0,35</span></div>
      <div><b>{contexts.length}</b><span>опухолевых контекстов</span></div>
      <div><b>{contexts.reduce((sum, item) => sum + item.models.length, 0)}</b><span>связей с клеточными линиями</span></div>
    </div>

    <section className={styles.section} id="admet">
      <div className={styles.sectionHeading}><div><span className={styles.eyebrow}>ПРОГНОЗНЫЙ СЛОЙ</span><h2>ADMET и ограничения</h2></div><span className={styles.pillNeutral}>Не измеренная безопасность</span></div>
      <div className={styles.admetGrid}>
        <div><span>QED</span><b>{num(profile.properties.qed, 2)}</b><small>Описательное свойство структуры</small></div>
        <div><span>Интегральный прогноз</span><b>{profile.admet.integrated_risk === "high_predicted_risk" ? "Повышенный риск по модели" : profile.admet.integrated_risk ? "Требует проверки" : "Нет данных"}</b><small>Прогноз, не эксперимент</small></div>
        <div><span>Ames для данной PYZ</span><b>{profile.admet.individual_ames || "Нет индивидуального значения"}</b><small>Сводка по серии — не личный показатель PYZ</small></div>
        <div><span>DILI для данной PYZ</span><b>{profile.admet.individual_dili || "Нет индивидуального значения"}</b><small>Прогноз по серии нельзя переносить</small></div>
      </div>
      <p className={styles.caution}>{profile.admet.note_ru || "Индивидуальные прогнозные данные отсутствуют."} Никакой прогноз не заменяет проверку в опыте. Источник: {profile.admet.source_report || "не указан"}.</p>
    </section>

    <section className={styles.section} id="targets">
      <div className={styles.sectionHeading}><div><span className={styles.eyebrow}>СТРУКТУРНОЕ СРАВНЕНИЕ · ДЕЙСТВУЮЩИЕ 13 МИШЕНЕЙ</span><h2>Предполагаемые белковые мишени</h2></div>
        <Link className={styles.secondaryAction} href={"/pyz/matrix?q=" + encodeURIComponent(profile.compound_id)}>Показать в матрице</Link></div>
      <p className={styles.description}>Танимото (Morgan r2 / 2048) — мера химического сходства, а не вероятность связывания.
        Порог 0,35 предназначен только для отбора гипотез.</p>
      {targets.length ? <div className={styles.targetList}>{targets.map(item =>
        <TargetCard key={item.target_gene} item={item} id={profile.compound_id} focus={item.target_gene.toUpperCase() === focus.toUpperCase() && !!focus} />
      )}</div> : <div className={styles.empty}>Локальный реестр структурных аналогов не доступен. Это не означает, что соединение неактивно.</div>}
    </section>

    <section className={styles.section} id="models">
      <div className={styles.sectionHeading}><div><span className={styles.eyebrow}>КОНТЕКСТ МИШЕНИ · НЕ ДЕЙСТВИЕ НАШЕЙ МОЛЕКУЛЫ</span><h2>Опухоли и клеточные модели</h2></div></div>
      <p className={styles.description}>Клеточные модели предложены на основании генетической зависимости от предполагаемых мишеней и результатов с известными препаратами. Они не являются доказательством чувствительности к PYZ.</p>
      {contexts.length ? <div className={styles.contextList}>{contexts.map((item, i) =>
        <ContextCard key={String(item.target_gene) + "-" + String(item.mcl_cancer_id || i)} item={item} />
      )}</div> : <div className={styles.empty}>Для этого соединения нет обоснованного модельного маршрута в текущих данных.</div>}
      <p className={styles.caution}>HEK293/HEK293T/HEK293FT — трансформированные исследовательские линии, не опухоль почки; MRC-5 — неопухолевые фибробласты лёгкого, не универсальное доказательство безопасности.</p>
    </section>

    <section className={styles.section} id="sources">
      <div className={styles.sectionHeading}><div><span className={styles.eyebrow}>ПРОИСХОЖДЕНИЕ И ГРАНИЦЫ ДАННЫХ</span><h2>Источники</h2></div></div>
      <dl className={styles.sources}>
        <div><dt>Реестр структуры</dt><dd>{profile.provenance.own_compounds_registry}</dd></div>
        <div><dt>ADMET</dt><dd>{profile.provenance.admet_report || "нет данных"}</dd></div>
        <div><dt>Лиганды и сходство</dt><dd>{profile.provenance.target_library}</dd></div>
        <div><dt>Опухолевый контекст</dt><dd>{profile.provenance.candidate_context}</dd></div>
        <div><dt>BindingDB P0.3.6</dt><dd>Не подключён: научная приёмка не завершена</dd></div>
        <div><dt>Докинг PYZ</dt><dd>Внешние результаты ожидаются</dd></div>
      </dl>
    </section>
    <footer className={styles.footer}>MCL Explorer · PYZ 2.0 · <Link href="/methodology">Методика и термины</Link></footer>
  </main>;
}
