"use client";

import Link from "next/link";
import { useMemo, useState } from "react";

type Model = Record<string, any>;

export default function CRISPRModelTable({
  models,
  initialOrgan = "",
  initialCancer = "",
}: {
  models: Model[];
  initialOrgan?: string;
  initialCancer?: string;
}) {
  const [query, setQuery] = useState("");
  const [organ, setOrgan] = useState(initialOrgan);
  const [cancer, setCancer] = useState(initialCancer);
  const [subtype, setSubtype] = useState("");

  const organs = useMemo(() => {
    const map = new Map<string, string>();
    models.forEach((m) => {
      if (m.mcl_organ_id) map.set(String(m.mcl_organ_id), String(m.mcl_organ_ru || m.oncotree_lineage || m.mcl_organ_id));
    });
    return [...map.entries()].sort((a, b) => a[1].localeCompare(b[1], "ru"));
  }, [models]);

  const cancers = useMemo(() => {
    const map = new Map<string, string>();
    models.filter((m) => !organ || m.mcl_organ_id === organ).forEach((m) => {
      if (m.mcl_cancer_id) map.set(String(m.mcl_cancer_id), String(m.mcl_cancer_name || m.oncotree_primary_disease || m.mcl_cancer_id));
    });
    return [...map.entries()].sort((a, b) => a[1].localeCompare(b[1], "ru"));
  }, [models, organ]);

  const subtypes = useMemo(() => {
    const values = new Set<string>();
    models
      .filter((m) => (!organ || m.mcl_organ_id === organ) && (!cancer || m.mcl_cancer_id === cancer))
      .forEach((m) => { if (m.mcl_subtype_name) values.add(String(m.mcl_subtype_name)); });
    return [...values].sort((a, b) => a.localeCompare(b, "ru"));
  }, [models, organ, cancer]);

  const filtered = useMemo(() => {
    const needle = query.trim().toLowerCase();
    return models.filter((m) => {
      if (organ && m.mcl_organ_id !== organ) return false;
      if (cancer && m.mcl_cancer_id !== cancer) return false;
      if (subtype && m.mcl_subtype_name !== subtype) return false;
      if (!needle) return true;
      const haystack = [
        m.cell_line_name, m.model_id, m.mcl_organ_ru, m.mcl_cancer_name,
        m.mcl_subtype_name, m.oncotree_code, m.depmap_model_type,
      ].filter(Boolean).join(" ").toLowerCase();
      return haystack.includes(needle);
    });
  }, [models, query, organ, cancer, subtype]);

  function changeOrgan(value: string) {
    setOrgan(value);
    setCancer("");
    setSubtype("");
  }

  function changeCancer(value: string) {
    setCancer(value);
    setSubtype("");
  }

  return (
    <div className="model-browser">
      <div className="model-browser-toolbar">
        <label className="model-search">
          <span>Поиск</span>
          <input value={query} onChange={(e) => setQuery(e.target.value)} placeholder="Линия, DepMap ID, опухоль…" />
        </label>
        <label className="model-search">
          <span>Орган / система</span>
          <select value={organ} onChange={(e) => changeOrgan(e.target.value)}>
            <option value="">Все</option>
            {organs.map(([id, label]) => <option key={id} value={id}>{label}</option>)}
          </select>
        </label>
        <label className="model-search">
          <span>Опухоль</span>
          <select value={cancer} onChange={(e) => changeCancer(e.target.value)}>
            <option value="">Все</option>
            {cancers.map(([id, label]) => <option key={id} value={id}>{label}</option>)}
          </select>
        </label>
        <label className="model-search">
          <span>Подтип</span>
          <select value={subtype} onChange={(e) => setSubtype(e.target.value)}>
            <option value="">Все</option>
            {subtypes.map((value) => <option key={value} value={value}>{value}</option>)}
          </select>
        </label>
      </div>

      <div className="model-browser-meta">Показано {filtered.length} из {models.length} CRISPR-моделей</div>

      <div className="table-wrap model-table-wrap">
        <table className="model-table">
          <thead>
            <tr><th>Клеточная линия</th><th>Орган / система</th><th>Опухоль</th><th>Подтип</th><th>В MCL</th><th></th></tr>
          </thead>
          <tbody>
            {filtered.map((m) => (
              <tr key={m.model_id}>
                <td><div className="model-name"><strong>{m.cell_line_name || m.model_id}</strong><small>{m.model_id}</small></div></td>
                <td>{m.mcl_organ_ru || m.oncotree_lineage || "—"}</td>
                <td><div className="variant-copy">{m.mcl_cancer_name || m.oncotree_primary_disease || "—"}</div></td>
                <td><div className="variant-copy">{m.mcl_subtype_name || m.oncotree_subtype || "—"}</div></td>
                <td>{Number(m.curated_contexts_n || 0) > 0 ? <span className="status-inline yes">{m.curated_contexts_n} контекст.</span> : <span className="status-inline no">ещё не анализировалась</span>}</td>
                <td><Link className="row-link" href={`/models/${encodeURIComponent(m.model_id)}`}>Открыть →</Link></td>
              </tr>
            ))}
          </tbody>
        </table>
      </div>
      {filtered.length === 0 && <div className="empty-state">По текущим фильтрам CRISPR-модели не найдены.</div>}
    </div>
  );
}
