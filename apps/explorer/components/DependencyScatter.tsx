"use client";

import { useMemo, useState } from "react";
import { useRouter } from "next/navigation";
import {
  CartesianGrid,
  LabelList,
  ReferenceLine,
  ResponsiveContainer,
  Scatter,
  ScatterChart,
  Tooltip,
  XAxis,
  YAxis,
} from "recharts";
import styles from "./DependencyLandscape.module.css";

type Point = Record<string, any>;
type Mode = "landscape" | "contextual" | "top";

function finite(value: unknown): value is number {
  return typeof value === "number" && Number.isFinite(value);
}

function tooltipValue(value: unknown) {
  if (value === null || value === undefined || value === "") return "—";
  const n = Number(value);
  if (!Number.isFinite(n)) return String(value);
  return n.toFixed(3);
}

function CustomTooltip({ active, payload }: any) {
  if (!active || !payload?.length) return null;
  const p = payload[0]?.payload || {};
  return (
    <div className="chart-tooltip">
      <div className="chart-tooltip-gene">{p.gene_symbol}</div>
      <div><span>Целевая группа</span><b>{tooltipValue(p.context_median_gene_effect)}</b></div>
      <div><span>Группа сравнения</span><b>{tooltipValue(p.comparator_median_gene_effect)}</b></div>
      <div><span>Δ Gene Effect</span><b>{tooltipValue(p.delta_gene_effect)}</b></div>
      <div><span>Размер эффекта, Cliff&apos;s δ</span><b>{tooltipValue(p.cliffs_delta)}</b></div>
      <div><span>q-value / FDR</span><b>{tooltipValue(p.q_value)}</b></div>
      {String(p.present_all_thresholds).toLowerCase() === "true" && <div className="tooltip-stable">Устойчивый повторяющийся ген</div>}
    </div>
  );
}

export default function DependencyScatter({
  points,
  total,
  compact = false,
}: {
  points: Point[];
  total?: number;
  compact?: boolean;
}) {
  const router = useRouter();
  const [mode, setMode] = useState<Mode>("landscape");
  const [query, setQuery] = useState("");

  const clean = useMemo(
    () => points.filter(
      (p) => finite(p.context_median_gene_effect) && finite(p.comparator_median_gene_effect),
    ),
    [points],
  );

  const extent = useMemo(() => {
    const values = clean.flatMap((p) => [p.context_median_gene_effect, p.comparator_median_gene_effect]);
    const min = values.length ? Math.min(...values) : -2;
    const max = values.length ? Math.max(...values) : 1;
    const pad = Math.max((max - min) * 0.06, 0.05);
    return [min - pad, max + pad] as [number, number];
  }, [clean]);

  const selected = useMemo(() => {
    if (compact || mode === "landscape") return clean;
    if (mode === "contextual") {
      return clean.filter((p) => Number(p.delta_gene_effect) < 0);
    }
    return [...clean]
      .filter((p) => Number.isFinite(Number(p.delta_gene_effect)))
      .sort((a, b) => Number(a.delta_gene_effect) - Number(b.delta_gene_effect))
      .slice(0, 200);
  }, [clean, compact, mode]);

  const queryNormalized = query.trim().toUpperCase();
  const highlighted = useMemo(
    () => queryNormalized
      ? selected.filter((p) => String(p.gene_symbol || "").toUpperCase().includes(queryNormalized))
      : [],
    [selected, queryNormalized],
  );
  const highlightedSymbols = useMemo(
    () => new Set(highlighted.map((p) => String(p.gene_symbol))),
    [highlighted],
  );

  const stable = useMemo(
    () => selected.filter(
      (p) => !highlightedSymbols.has(String(p.gene_symbol)) && String(p.present_all_thresholds).toLowerCase() === "true",
    ),
    [selected, highlightedSymbols],
  );
  const regular = useMemo(
    () => selected.filter(
      (p) => !highlightedSymbols.has(String(p.gene_symbol)) && String(p.present_all_thresholds).toLowerCase() !== "true",
    ),
    [selected, highlightedSymbols],
  );

  function openGene(event: any) {
    const payload = event?.payload || event;
    const symbol = String(payload?.gene_symbol || "").trim();
    if (symbol) router.push(`/genes/${encodeURIComponent(symbol)}`);
  }

  const height = compact ? 270 : 430;

  return (
    <div className={`chart-shell ${compact ? styles.compactShell : ""}`}>
      {!compact && <>
        <div className={styles.toolbar}>
          <div className={styles.modeGroup} role="group" aria-label="Режим карты зависимостей">
            <button type="button" className={`${styles.modeButton} ${mode === "landscape" ? styles.modeButtonActive : ""}`} onClick={() => setMode("landscape")}>Все показанные</button>
            <button type="button" className={`${styles.modeButton} ${mode === "contextual" ? styles.modeButtonActive : ""}`} onClick={() => setMode("contextual")}>Δ &lt; 0</button>
            <button type="button" className={`${styles.modeButton} ${mode === "top" ? styles.modeButtonActive : ""}`} onClick={() => setMode("top")}>Top 200</button>
          </div>
          <label className={styles.search}>
            <span>Найти ген</span>
            <input value={query} onChange={(e) => setQuery(e.target.value)} placeholder="KRAS, DNM1L, FIS1…" />
          </label>
        </div>
        <div className={styles.metaRow}>
          <span><b>{selected.length}</b> точек показано{total ? ` · ${total} генов прошли серверные фильтры` : ""}</span>
          <span className={styles.hint}>Нажмите на точку, чтобы открыть карточку гена</span>
        </div>
      </>}

      <div className={`chart-legend ${compact ? styles.compactLegend : ""}`}>
        <span><i className="dot regular" /> гены после фильтров</span>
        <span><i className="dot stable" /> устойчивые повторяющиеся гены</span>
        {highlighted.length > 0 && <span className={styles.highlightLegend}><i className={styles.highlightDot} /> найденные гены</span>}
        {!compact && <span className="muted">Ниже диагонали — зависимость сильнее в целевой группе</span>}
      </div>
      <div style={{ width: "100%", height }}>
        <ResponsiveContainer>
          <ScatterChart margin={{ top: 16, right: 18, bottom: compact ? 16 : 30, left: 12 }}>
            <CartesianGrid strokeDasharray="3 4" stroke="#e8edf3" />
            <XAxis
              type="number"
              dataKey="comparator_median_gene_effect"
              domain={extent}
              tick={{ fontSize: compact ? 9 : 11, fill: "#748398" }}
              axisLine={{ stroke: "#d8e0e9" }}
              tickLine={{ stroke: "#d8e0e9" }}
              label={compact ? undefined : { value: "Медианный Gene Effect · группа сравнения", position: "insideBottom", offset: -17, fill: "#63758b", fontSize: 11 }}
            />
            <YAxis
              type="number"
              dataKey="context_median_gene_effect"
              domain={extent}
              tick={{ fontSize: compact ? 9 : 11, fill: "#748398" }}
              axisLine={{ stroke: "#d8e0e9" }}
              tickLine={{ stroke: "#d8e0e9" }}
              label={compact ? undefined : { value: "Медианный Gene Effect · целевая группа", angle: -90, position: "insideLeft", fill: "#63758b", fontSize: 11 }}
            />
            <ReferenceLine segment={[{ x: extent[0], y: extent[0] }, { x: extent[1], y: extent[1] }]} stroke="#9eacbc" strokeDasharray="5 5" />
            <Tooltip content={<CustomTooltip />} cursor={{ stroke: "#bdc9d7", strokeDasharray: "3 3" }} />
            <Scatter data={regular} fill="#7f9fca" fillOpacity={0.42} isAnimationActive={false} onClick={openGene} style={{ cursor: "pointer" }} />
            <Scatter data={stable} fill="#134a7d" fillOpacity={0.92} isAnimationActive={false} onClick={openGene} style={{ cursor: "pointer" }} />
            {highlighted.length > 0 && <Scatter data={highlighted} fill="#d97706" fillOpacity={1} isAnimationActive={false} onClick={openGene} style={{ cursor: "pointer" }}>
              <LabelList dataKey="gene_symbol" position="top" fill="#9a4f05" fontSize={10} fontWeight={800} />
            </Scatter>}
          </ScatterChart>
        </ResponsiveContainer>
      </div>
      {compact && <div className={styles.previewFooter}><span>Каждая точка — ген. Ниже диагонали находятся более выраженные зависимости целевой группы.</span><strong>{clean.length} точек в preview</strong></div>}
    </div>
  );
}
