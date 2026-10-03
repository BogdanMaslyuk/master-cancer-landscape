import Link from "next/link";
import { apiGet } from "../../../lib/api";

type Model = Record<string, any>;

function variantLabel(variants:any[]|undefined){
  if(!variants?.length) return "Определяющие варианты не зарегистрированы";
  return variants.map(v=>[v.gene,v.protein_change].filter(Boolean).join(" ")).filter(Boolean).join(" · ");
}

export default async function ModelPage({params}:{params:Promise<{id:string}>}){
  const {id}=await params;
  const model=await apiGet<Model>(`/api/models/${encodeURIComponent(decodeURIComponent(id))}`);
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
        <div><span>Секвенирование</span><b>{model.sequencing_available ? "доступно" : "нет"}</b></div>
      </div>
    </section>

    <section className="section split">
      <div className="card">
        <div className="eyebrow">ПРОИСХОЖДЕНИЕ</div>
        <h3>{model.oncotree_subtype || model.oncotree_primary_disease || "—"}</h3>
        <div className="definition-list">
          <div><b>Орган / lineage</b><span>{model.oncotree_lineage || "—"}</span></div>
          <div><b>Первичное заболевание</b><span>{model.oncotree_primary_disease || "—"}</span></div>
          <div><b>Подтип</b><span>{model.oncotree_subtype || "—"}</span></div>
        </div>
      </div>
      <div className="callout">
        <div className="eyebrow">КАК ИСПОЛЬЗУЕТСЯ В MCL</div>
        <h3>Роль модели зависит от научного контекста</h3>
        <p className="section-copy">Одна клеточная линия может быть целевой, контрольной или исключённой в зависимости от того, какую молекулярную гипотезу мы проверяем. Ниже показаны все её текущие назначения.</p>
      </div>
    </section>

    <section className="section">
      <div className="section-header"><div><div className="eyebrow">УЧАСТИЕ В КОНТЕКСТАХ</div><h2>Где используется эта линия?</h2><div className="section-copy">Для каждого контекста показаны роль линии и молекулярные варианты, которые повлияли на назначение.</div></div></div>
      <div className="membership-grid">
        {model.memberships?.map((m:any)=><article className="membership-card" key={`${m.cancer_id}-${m.model_id}`}>
          <div className="membership-card-top"><span className={`group-badge ${m.assigned_group}`}>{m.assigned_group_ru}</span><small>{m.cancer_id}</small></div>
          <h3>{m.cancer_ru}</h3>
          <div className="membership-molecular">{m.molecular_context || "Молекулярный контекст не указан"}</div>
          <div className="membership-variants"><span>Зарегистрированные варианты</span><b>{variantLabel(m.variants)}</b></div>
          {m.assignment_reason && <p className="technical-note">Причина назначения: {m.assignment_reason}</p>}
          <Link className="primary-link" href={`/atlas/${m.cancer_id}`}>Открыть опухолевый контекст →</Link>
        </article>)}
      </div>
    </section>

    <section className="section technical-note"><b>Важно:</b> эта карточка отображает только данные, уже присутствующие в текущем DepMap-аудите MCL. Полная молекулярная характеристика клеточной линии будет расширяться по мере подключения новых слоёв данных.</section>
  </>;
}
