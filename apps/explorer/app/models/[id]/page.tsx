import Link from "next/link";
import ModelDependencyPanel from "../../../components/ModelDependencyPanel";
import { ModelMultiOmicsPanel } from "../../../components/MultiOmicsPanels";
import PharmacologyPanel from "../../../components/PharmacologyPanel";
import { apiGet } from "../../../lib/api";
import styles from "./page.module.css";
import readiness from "./readiness.module.css";

type Model = Record<string, any>;
type MultiOmics = Record<string, any>;
type Dependencies = Record<string, any>;
type Pharmacology = Record<string, any>;
type Params = Record<string,string|string[]|undefined>;

function one(value:string|string[]|undefined){ return Array.isArray(value) ? value[0] : value; }

function variantLabel(variants:any[]|undefined){
  if(!variants?.length) return "Определяющие варианты не зарегистрированы";
  return variants.map(v=>[v.gene,v.protein_change].filter(Boolean).join(" ")).filter(Boolean).join(" · ");
}

function yesNo(value:any){ return value ? "да" : "—"; }

function impactRu(value:any){
  const key=String(value||"").toUpperCase();
  return ({HIGH:"высокое",MODERATE:"умеренное",LOW:"низкое",MODIFIER:"модификатор"} as Record<string,string>)[key] || value || "—";
}

export default async function ModelPage({params,searchParams}:{params:Promise<{id:string}>;searchParams:Promise<Params>}){
  const [{id},sp]=await Promise.all([params,searchParams]);
  const modelId=decodeURIComponent(id);
  const depQ=one(sp.dep_q)||"";
  const depType=one(sp.dep_type)||"";
  const depDomain=one(sp.dep_domain)||"";
  const depOffset=Math.max(0,Number(one(sp.dep_offset)||"0")||0);
  const depQuery=new URLSearchParams({limit:"200",offset:String(depOffset)});
  if(depQ)depQuery.set("search",depQ);
  if(depType)depQuery.set("dependency_type",depType);
  if(depDomain)depQuery.set("domain",depDomain);

  const [model, omics, dependencies, pharmacology]=await Promise.all([
    apiGet<Model>(`/api/models/${encodeURIComponent(modelId)}`),
    apiGet<MultiOmics>(`/api/models/${encodeURIComponent(modelId)}/multiomics?limit=30`),
    apiGet<Dependencies>(`/api/models/${encodeURIComponent(modelId)}/dependencies?${depQuery.toString()}`),
    apiGet<Pharmacology>(`/api/models/${encodeURIComponent(modelId)}/pharmacology?limit=100`),
  ]);
  const genetics=model.genetics || {};
  const metadata=model.metadata || {};
  const priority=genetics.priority_variants || [];
  const membershipsN=model.memberships?.length || 0;
  const fullGenetics=genetics.availability === "full";

  return <>
    <div className="breadcrumbs"><Link href="/models">Исследователь моделей</Link><span>›</span><strong>{model.cell_line_name}</strong></div>

    <section className="model-hero">
      <div>
        <div className="eyebrow">КЛЕТОЧНАЯ МОДЕЛЬ · DEPMAP</div>
        <h1>{model.cell_line_name}</h1>
        <div className="model-id-pill">{model.model_id}</div>
      </div>
      <div className="model-hero-meta">
        <div><span>Тип модели</span><b>{model.depmap_model_type || "—"}</b></div>
        <div><span>OncoTree</span><b>{model.oncotree_code || model.oncotree_subtype || "—"}</b></div>
        <div><span>CRISPR Gene Effect</span><b>{dependencies.available ? "полный профиль" : "нужно индексировать"}</b></div>
      </div>
    </section>

    <section className={styles.identityStrip}>
      <div><span>Опухоль</span><b>{model.oncotree_subtype || model.oncotree_primary_disease || "—"}</b></div>
      <div><span>Орган / lineage</span><b>{model.oncotree_lineage || "—"}</b></div>
      <div><span>Происхождение модели</span><b>{metadata.primary_or_metastasis || metadata.source_type || "—"}</b></div>
      <div><span>Источник образца</span><b>{metadata.sample_collection_site || metadata.tissue_origin || "—"}</b></div>
    </section>

    <section className="section">
      <div className="section-header">
        <div>
          <div className="eyebrow">ФУНКЦИОНАЛЬНЫЙ CRISPR-ПРОФИЛЬ</div>
          <h2>Что поддерживает жизнеспособность этой модели?</h2>
          <div className="section-copy">Все гены с доступным Chronos Gene Effect расположены от наиболее сильного влияния CRISPR-нокаута на рост и выживание клетки к наиболее слабому. Рядом показано, насколько зависимость распространена среди других моделей DepMap и к какому функциональному кластеру относится ген.</div>
        </div>
      </div>
      <ModelDependencyPanel
        payload={dependencies}
        modelId={modelId}
        filters={{q:depQ,type:depType,domain:depDomain,offset:depOffset}}
      />
    </section>

    <section className={styles.scopeWarning}>
      <span>Важно</span>
      <div><b>Сильная CRISPR-зависимость не равна готовой лекарственной мишени</b><p>Gene Effect показывает последствия потери функции гена в экспериментальной модели. Для выбора терапевтической мишени отдельно оцениваются селективность для опухоли, нормальные ткани, лекарственная доступность белка и воспроизводимость эффекта.</p></div>
    </section>

    <section className="section">
      <div className="section-header">
        <div>
          <div className="eyebrow">ФАРМАКОЛОГИЯ МОДЕЛИ</div>
          <h2>Какие вещества уже проверялись на этой клеточной линии?</h2>
          <div className="section-copy">MCL хранит экспериментальный ответ клетки отдельно от известных мишеней вещества. Если мишень аннотирована, рядом показывается Gene Effect этого белка в той же модели — это позволяет быстро увидеть согласованность фармакологии и CRISPR, не выдавая её за доказанный механизм.</div>
        </div>
      </div>
      <PharmacologyPanel payload={pharmacology}/>
    </section>

    <section className="section">
      <div className="section-header">
        <div>
          <div className="eyebrow">ПОКРЫТИЕ ДАННЫХ</div>
          <h2>Что мы реально знаем об этой модели?</h2>
          <div className="section-copy">MCL разделяет доступные и ещё не подключённые слои, чтобы отсутствие одного типа данных не выглядело как отрицательный биологический результат.</div>
        </div>
      </div>
      <div className={readiness.grid}>
        <article className={readiness.card}>
          <div className={readiness.cardTop}><span className={readiness.index}>01</span><span className={`${readiness.badge} ${readiness.available}`}>доступно</span></div>
          <h3>Идентичность модели</h3>
          <p>DepMap ID, тип модели, OncoTree-аннотация и опухолевое происхождение.</p>
          <strong>{model.oncotree_code || model.depmap_model_type || "DepMap model"}</strong>
        </article>
        <article className={readiness.card}>
          <div className={readiness.cardTop}><span className={readiness.index}>02</span><span className={`${readiness.badge} ${dependencies.available ? readiness.available : readiness.partial}`}>{dependencies.available ? "доступно" : "нужно перестроить"}</span></div>
          <h3>CRISPR-зависимости</h3>
          <p>Полногеномный Chronos Gene Effect с рангом каждого гена и оценкой распространённости зависимости.</p>
          <strong>{dependencies.available ? `${Number(dependencies.genes_measured_n||0).toLocaleString("ru-RU")} генов` : "Индекс не готов"}</strong>
        </article>
        <article className={readiness.card}>
          <div className={readiness.cardTop}><span className={readiness.index}>03</span><span className={`${readiness.badge} ${fullGenetics ? readiness.available : readiness.partial}`}>{fullGenetics ? "полный профиль" : "ограничено"}</span></div>
          <h3>Соматические варианты</h3>
          <p>Индивидуальный мутационный фон из DepMap OmicsSomaticMutations, а не только определяющая мутация контекста.</p>
          <strong>{genetics.mutations_n ?? 0} записей · {genetics.mutated_genes_n ?? 0} генов</strong>
        </article>
        <article className={readiness.card}>
          <div className={readiness.cardTop}><span className={readiness.index}>04</span><span className={`${readiness.badge} ${pharmacology.available ? readiness.available : readiness.pending}`}>{pharmacology.available ? "доступно" : "не загружено"}</span></div>
          <h3>Фармакология</h3>
          <p>Экспериментальные ответы на вещества и их отдельные target-аннотации с CRISPR-согласованностью.</p>
          <strong>{pharmacology.available ? `${Number(pharmacology.compounds_n||0).toLocaleString("ru-RU")} веществ` : "Источник нужно подключить"}</strong>
        </article>
      </div>
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

    <section className="section">
      <div className="section-header">
        <div>
          <div className="eyebrow">MULTI-OMICS ПРОФИЛЬ МОДЕЛИ</div>
          <h2>Что ещё объясняет индивидуальную зависимость?</h2>
          <div className="section-copy">Совмещаем индивидуальный CRISPR Gene Effect, RNA expression и относительное copy number, чтобы видеть функциональную зависимость на фоне экспрессии и геномной дозировки того же гена.</div>
        </div>
      </div>
      <ModelMultiOmicsPanel omics={omics}/>
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

    <section className="section technical-note"><b>Граница интерпретации:</b> мутации, RNA expression, относительное copy number, CRISPR Gene Effect и лекарственная чувствительность описывают экспериментальную модель. Они не заменяют пациентские когорты; совпадение чувствительности к препарату с CRISPR-зависимостью аннотированной мишени поддерживает, но не доказывает механизм действия.</section>
  </>;
}
