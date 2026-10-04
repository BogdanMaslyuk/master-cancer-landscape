"use client";

import Link from "next/link";
import { useMemo, useState } from "react";
import type { ModelListItem, ModelVariant } from "../lib/generated/api-types";

const tabs = [
  { id: "all", label: "Все" },
  { id: "context", label: "Целевая группа" },
  { id: "comparator", label: "Группа сравнения" },
  { id: "excluded", label: "Исключённые" },
];

function variantLabel(variants: ModelVariant[] | undefined) {
  if (!variants?.length) return "Определяющие варианты не зарегистрированы";
  return variants
    .map((v) => [v.gene, v.protein_change].filter(Boolean).join(" "))
    .filter(Boolean)
    .join(" · ");
}

export default function CellModelTable({ models, showContext = false }: { models: ModelListItem[]; showContext?: boolean }) {
  const [tab, setTab] = useState("all");
  const [query, setQuery] = useState("");

  const filtered = useMemo(() => {
    const needle = query.trim().toLowerCase();
    return models.filter((m) => {
      const inTab = tab === "all" || m.assigned_group === tab;
      if (!inTab) return false;
      if (!needle) return true;
      const haystack = `${m.cell_line_name || ""} ${m.model_id} ${m.cancer_ru || ""} ${m.organ_ru || ""} ${m.molecular_context || ""} ${variantLabel(m.variants)}`.toLowerCase();
      return haystack.includes(needle);
    });
  }, [models, query, tab]);

  const counts = useMemo(() => ({
    all: models.length,
    context: models.filter((m) => m.assigned_group === "context").length,
    comparator: models.filter((m) => m.assigned_group === "comparator").length,
    excluded: models.filter((m) => m.assigned_group === "excluded").length,
  }), [models]);

  return (
    <div className="model-browser">
      <div className="model-browser-toolbar">
        <div className="segmented-control" role="tablist" aria-label="Фильтр клеточных моделей">
          {tabs.map((item) => (
            <button
              type="button"
              role="tab"
              aria-selected={tab === item.id}
              className={tab === item.id ? "active" : ""}
              onClick={() => setTab(item.id)}
              key={item.id}
            >
              {item.label}<span>{counts[item.id as keyof typeof counts]}</span>
            </button>
          ))}
        </div>
        <label className="model-search">
          <span>Поиск</span>
          <input value={query} onChange={(e) => setQuery(e.target.value)} placeholder="Линия, DepMap ID, опухоль, мутация…" />
        </label>
      </div>

      <div className="model-browser-meta">Показано {filtered.length} из {models.length} записей</div>

      <div className="table-wrap model-table-wrap">
        <table className="model-table">
          <thead><tr><th>Клеточная линия</th>{showContext && <th>Опухолевый контекст</th>}<th>Роль</th><th>Молекулярные варианты</th><th>Секвенирование</th><th></th></tr></thead>
          <tbody>
            {filtered.map((m) => (
              <tr key={`${m.cancer_id}-${m.model_id}`} className={m.assigned_group === "excluded" ? "row-muted" : ""}>
                <td><div className="model-name"><strong>{m.cell_line_name || m.model_id}</strong><small>{m.model_id}</small></div></td>
                {showContext && <td><Link className="context-cell-link" href={`/atlas/${m.cancer_id}`}><strong>{m.cancer_ru || m.cancer_id}</strong><small>{m.molecular_context || "—"}</small></Link></td>}
                <td><span className={`group-badge ${m.assigned_group}`}>{m.assigned_group_ru || m.assigned_group}</span></td>
                <td><div className="variant-copy">{variantLabel(m.variants)}</div></td>
                <td>{m.sequencing_available ? <span className="status-inline yes">доступно</span> : <span className="status-inline no">нет</span>}</td>
                <td><Link className="row-link" href={`/models/${encodeURIComponent(m.model_id)}`}>Открыть →</Link></td>
              </tr>
            ))}
          </tbody>
        </table>
      </div>

      {filtered.length === 0 && <div className="empty-state">По текущим фильтрам клеточные модели не найдены.</div>}
    </div>
  );
}
