"use client";

import { useMemo } from "react";
import {
  CartesianGrid,
  ReferenceLine,
  ResponsiveContainer,
  Scatter,
  ScatterChart,
  Tooltip,
  XAxis,
  YAxis,
} from "recharts";

type Point = Record<string, any>;

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
      {String(p.present_all_thresholds).toLowerCase() === "true" && <div className="tooltip-stable">Устойчивый повторяющийся ген</div>}
    </div>
  );
}

export default function DependencyScatter({ points }: { points: Point[] }) {
  const { regular, stable, extent } = useMemo(() => {
    const clean = points.filter(
      (p) => finite(p.context_median_gene_effect) && finite(p.comparator_median_gene_effect),
    );
    const stablePoints = clean.filter(
      (p) => String(p.present_all_thresholds).toLowerCase() === "true",
    );
    const regularPoints = clean.filter(
      (p) => String(p.present_all_thresholds).toLowerCase() !== "true",
    );
    const values = clean.flatMap((p) => [p.context_median_gene_effect, p.comparator_median_gene_effect]);
    const min = values.length ? Math.min(...values) : -2;
    const max = values.length ? Math.max(...values) : 1;
    const pad = Math.max((max - min) * 0.06, 0.05);
    return { regular: regularPoints, stable: stablePoints, extent: [min - pad, max + pad] as [number, number] };
  }, [points]);

  return (
    <div className="chart-shell">
      <div className="chart-legend">
        <span><i className="dot regular" /> гены после фильтров</span>
        <span><i className="dot stable" /> устойчивые повторяющиеся гены</span>
        <span className="muted">Ниже диагонали — зависимость сильнее в целевой группе</span>
      </div>
      <div style={{ width: "100%", height: 430 }}>
        <ResponsiveContainer>
          <ScatterChart margin={{ top: 16, right: 18, bottom: 30, left: 12 }}>
            <CartesianGrid strokeDasharray="3 4" stroke="#e8edf3" />
            <XAxis
              type="number"
              dataKey="comparator_median_gene_effect"
              domain={extent}
              tick={{ fontSize: 11, fill: "#748398" }}
              axisLine={{ stroke: "#d8e0e9" }}
              tickLine={{ stroke: "#d8e0e9" }}
              label={{ value: "Медианный Gene Effect · группа сравнения", position: "insideBottom", offset: -17, fill: "#63758b", fontSize: 11 }}
            />
            <YAxis
              type="number"
              dataKey="context_median_gene_effect"
              domain={extent}
              tick={{ fontSize: 11, fill: "#748398" }}
              axisLine={{ stroke: "#d8e0e9" }}
              tickLine={{ stroke: "#d8e0e9" }}
              label={{ value: "Медианный Gene Effect · целевая группа", angle: -90, position: "insideLeft", fill: "#63758b", fontSize: 11 }}
            />
            <ReferenceLine segment={[{ x: extent[0], y: extent[0] }, { x: extent[1], y: extent[1] }]} stroke="#9eacbc" strokeDasharray="5 5" />
            <Tooltip content={<CustomTooltip />} cursor={{ stroke: "#bdc9d7", strokeDasharray: "3 3" }} />
            <Scatter data={regular} fill="#7f9fca" fillOpacity={0.42} isAnimationActive={false} />
            <Scatter data={stable} fill="#134a7d" fillOpacity={0.9} isAnimationActive={false} />
          </ScatterChart>
        </ResponsiveContainer>
      </div>
    </div>
  );
}
