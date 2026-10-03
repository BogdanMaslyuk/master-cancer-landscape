import Link from "next/link";
import { apiGet } from "../../../lib/api";
import styles from "./page.module.css";

type Model = Record<string, any>;

function variantLabel(variants:any[]|undefined){
  if(!variants?.length) return "Определяющие варианты не зарегистрированы";
  return variants.map(v=>[v.gene,v.protein_change].filter(Boolean).join(" ")).filter(Boolean).join(" · ");
}

function yesNo(value:any){ return value ? "да" : "—"; }

function impactRu(value:any){
  const key=String(value||"").toUpperCase();
  return ({HIGH:"высокое",MODERATE:"умеренное",LOW:"низкое",MODIFIER:"модификатор"} as Record<string,string>)[key] || value || "—";
}

export default async function ModelPage({params}:{params:Promise<{id:string}>}){
  const {id}=await params;
  const model=await apiGet<Model>(`/api/models/${encodeURIComponent(decodeURIComponent(id))}`);
  const genetics=model.genetics || {};
  const metadata=model.metadata || {};
  const priority=genetics.priority_variants || [];

  return <>
    <div className="breadcrumbs"><Link href="/models">Клеточные линии</Link><span>›</span><strong>{model.cell_line_name}</strong></div>

    <section className="model-hero">
      <div>
        <div className="eyebrow">КЛЕТОЧНАЯ МОДЕЛЬ · DEPMAP</div>
        <h1>{model.cell_line_name}</h1>
        <div className="model-id-pill">{model.model_id}</div>
      </div>
      <div className="model-hero-meta">
        <div><span>Тип модели</span><b>{model.depmap_model_type || "—"}</b></div>
        <div><span>OncoTree</span><b>{model.oncotree_code || model.oncotree_subtype || "—"}</b></div>
        <div><span>Мутационное профилирование</span><b>{model.sequencing_available ? "доступно" : "нет"}</b></div>
      </div>
    </section>

    <section className={styles.identityStrip}>
      <div><span>Опухоль</span><b>{model.oncotree_subtype || model.oncotree_primary_disease || "—"}</b></div>
      <div><span>Орган / lineage</span><b>{model.oncotree_lineage || "—"}</b></div>
      <div><span>Происхождение модели</span><b>{metadata.primary_or_metastasis || metadata.source_type || "—"}</b></div>
      <div><span>Источник образца</span><b>{metadata.sample_collection_site || metadata.tissue_origin || "—"}</b></div>
    </section>

    <section className={styles.scopeWarning}>
      <span>Модель</span>
      <div><b>Ниже показана генетика именно {model.cell_line_name}</b><p>Эти варианты характеризуют конкретную экспериментальную линию. Их нельзя трактовать как частоту или типичный генетический профиль всех пациентов с {model.oncotree_subtype || model.oncotree_primary_disease}.</p></div>
    </section>

    <section className="section">
      <div className="section-header">
        <div>
          <div className="eyebrow">ИНДИВИДУАЛЬНЫЙ МОЛЕКУЛЯРНЫЙ ПРОФИЛЬ</div>
          <h2>Генетика клеточной линии</h2>
          <div className="section-copy">Сначала показываются варианты с наибольшей потенциальной функциональной значимостью: драйверные, hotspot, вероятная потеря функции и варианты с высоким влиянием по VEP.</div>
        </div>
        <span className={`${styles.profileState} ${genetics.availability === "full" ? styles.full : styles.limited}`}>{genetics.availability_ru || "—"}</span>
      </div>

      <div className={styles.geneticMetrics}>
        <div><strong>{genetics.mutations_n ?? 0}</strong><span>мутационных записей</span></div>
        <div><strong>{genetics.mutated_genes_n ?? 0}</strong><span>затронутых генов</span></div>
        <div><strong>{genetics.driver_mutations_n ?? 0}</strong><span>драйверных вариантов</span></div>
        <div><strong>{genetics.hotspot_mutations_n ?? 0}</strong><span>hotspot-вариантов</span></div>
        <div><strong>{genetics.likely_lof_n ?? 0}</strong><span>вероятная потеря функции</span></div>
        <div><strong>{genetics.high_impact_n ?? 0}</strong><span>высокое влияние VEP</span></div>
      </div>

      {genetics.availability !== "full" && <div className={styles.indexPrompt}>
        <div><b>Сейчас карточка видит только варианты, использованные для формирования исследовательских групп.</b><p>Полный DepMap-файл уже есть локально. Один раз создайте компактный индекс — после этого здесь появится индивидуальный генетический профиль линии.</p></div>
        <code>{genetics.build_command}</code>
      </div>}

      {priority.length > 0 ? <div className={styles.variantTableWrap}>
        <table className={styles.variantTable}>
          <thead><tr><th>Ген</th><th>Изменение</th><th>Тип / следствие</th><th>Влияние</th><th>Драйвер</th><th>Hotspot</th><th>Потеря функции</th><th>Доля аллеля</th></tr></thead>
          <tbody>{priority.map((v:any,index:number)=><tr key={`${v.gene}-${v.protein_change}-${v.dna_change}-${index}`}>
            <td><Link href={`/genes/${encodeURIComponent(v.gene || "")}`} className={styles.geneLink}>{v.gene || "—"}</Link></td>
            <td><b>{v.protein_change || v.dna_change || "—"}</b>{v.dna_change && v.protein_change && <small>{v.dna_change}</small>}</td>
            <td>{v.variant_type || v.molecular_consequence || v.classification || "—"}</td>
            <td><span className={`${styles.impact} ${String(v.vep_impact||"").toLowerCase()}`}>{impactRu(v.vep_impact)}</span></td>
            <td>{yesNo(v.driver)}</td>
            <td>{yesNo(v.hotspot)}</td>
            <td>{yesNo(v.likely_lof)}</td>
            <td>{v.allele_fraction !== null && v.allele_fraction !== undefined ? Number(v.allele_fraction).toFixed(3) : "—"}</td>
          </tr>)}</tbody>
        </table>
      </div> : <div className="empty-state"><strong>Приоритетные варианты не найдены</strong><span>Это не означает отсутствие всех генетических изменений: проверьте статус полноты профиля выше.</span></div>}

      <p className={styles.sourceNote}><b>Источник:</b> {genetics.source || "—"}. {genetics.note}</p>
    </section>

    <section className="section split">
      <div className="card">
        <div className="eyebrow">ПРОИСХОЖДЕНИЕ И ОБРАЗЕЦ</div>
        <h3>{model.oncotree_subtype || model.oncotree_primary_disease || "—"}</h3>
        <div className="definition-list">
          <div><b>Первичное заболевание</b><span>{model.oncotree_primary_disease || "—"}</span></div>
          <div><b>Подтип</b><span>{model.oncotree_subtype || "—"}</span></div>
          <div><b>Первичная / метастатическая</b><span>{metadata.primary_or_metastasis || "—"}</span></div>
          <div><b>Место получения образца</b><span>{metadata.sample_collection_site || "—"}</span></div>
          <div><b>Возраст донора</b><span>{metadata.age ?? "—"}</span></div>
          <div><b>Пол</b><span>{metadata.sex || "—"}</span></div>
          <div><b>CCLE</b><span>{metadata.ccle_name || "—"}</span></div>
          <div><b>Sanger ID</b><span>{metadata.sanger_model_id || "—"}</span></div>
        </div>
      </div>
      <div className="callout">
        <div className="eyebrow">КАК ЧИТАТЬ ЭТУ КАРТОЧКУ</div>
        <h3>Контекстная мутация — только одна часть фона</h3>
        <p className="section-copy">Даже если две линии имеют одинаковый KRAS G12C, их сопутствующие изменения могут различаться и менять CRISPR-зависимости. Поэтому выводы группы нужно проверять на уровне отдельных моделей.</p>
      </div>
    </section>

    <section className="section">
      <div className="section-header"><div><div className="eyebrow">УЧАСТИЕ В КОНТЕКСТАХ</div><h2>Почему эта линия попала в анализ?</h2><div className="section-copy">Для каждого контекста показаны роль линии и варианты, которые использовались именно для назначения в целевую или контрольную группу.</div></div></div>
      <div className="membership-grid">
        {model.memberships?.map((m:any)=><article className="membership-card" key={`${m.cancer_id}-${m.model_id}`}>
          <div className="membership-card-top"><span className={`group-badge ${m.assigned_group}`}>{m.assigned_group_ru}</span><small>{m.cancer_id}</small></div>
          <h3>{m.cancer_ru}</h3>
          <div className="membership-molecular">{m.molecular_context || "Молекулярный контекст не указан"}</div>
          <div className="membership-variants"><span>Варианты, повлиявшие на назначение</span><b>{variantLabel(m.variants)}</b></div>
          {m.assignment_reason && <p className="technical-note">Причина назначения: {m.assignment_reason}</p>}
          <Link className="primary-link" href={`/atlas/${m.cancer_id}`}>Вернуться к опухолевому контексту →</Link>
        </article>)}
      </div>
    </section>

    <section className="section technical-note"><b>Следующий функциональный слой:</b> к этой карточке планируется добавить индивидуальный профиль CRISPR Gene Effect, экспрессию, изменения числа копий и оценку того, насколько модель репрезентативна для соответствующей опухоли у пациентов.</section>
  </>;
}
