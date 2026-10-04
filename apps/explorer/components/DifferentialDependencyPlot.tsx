"use client";

import { useMemo } from "react";
import { useRouter } from "next/navigation";
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

function numberOrNull(value: unknown): number | null {
  const n = Number(value);
  return Number.isFinite(n) ? n : null;
}

function fmt(value: unknown, digits = 3) {
  const n = Number(value);
  if (!Number.isFinite(n)) return "—";
  return n.toFixed(digits);
}

function CustomTooltip({ active, payload }: any) {
  if (!active || !payload?.length) return null;
  const p = payload[0]?.payload || {};
  return <div className="chart-tooltip">
    <div className="chart-tooltip-gene">{p.gene_symbol}</div>
    <div><span>Δ Gene Effect</span><b>{fmt(p.delta_gene_effect)}</b></div>
    <div><span>q-value / FDR</span><b>{fmt(p.q_value, 5)}</b></div>
    <div><span>−log10(q)</span><b>{fmt(p._significance)}</b></div>
    <div><span>Cliff&apos;s δ</span><b>{fmt(p.cliffs_delta)}</b></div>
    {String(p.present_all_thresholds).toLowerCase() === "true" && <div className="tooltip-stable">Устойчивый повторяющийся ген</div>}
  </div>;
}

export default function DifferentialDependencyPlot({ points }: { points: Point[] }) {
  const router = useRouter();

  const { regular, significant, stable, xExtent, yMax } = useMemo(() => {
    const clean = points.map((p) => {
      const delta = numberOrNull(p.delta_gene_effect);
      const q = numberOrNull(p.q_value);
      if (delta === null || q === null || q <= 0) return null;
      return { ...p, _significance: -Math.log10(Math.max(q, 1e-300)) };
    }).filter(Boolean) as Point[];

    const stableRows = clean.filter((p) => String(p.present_all_thresholds).toLowerCase() === "true");
    const stableSymbols = new Set(stableRows.map((p) => String(p.gene_symbol)));
    const significantRows = clean.filter((p) => !stableSymbols.has(String(p.gene_symbol)) && Number(p.q_value) <= 0.05);
    const regularRows = clean.filter((p) => !stableSymbols.has(String(p.gene_symbol)) && Number(p.q_value) > 0.05);

    const xs = clean.map((p) => Number(p.delta_gene_effect));
    const min = xs.length ? Math.min(...xs, 0) : -1;
    const max = xs.length ? Math.max(...xs, 0) : 1;
    const pad = Math.max((max - min) * 0.07, 0.05);
    const maxY = clean.length ? Math.max(...clean.map((p) => Number(p._significance))) : 3;

    return {
      regular: regularRows,
      significant: significantRows,
      stable: stableRows,
      xExtent: [min - pad, max + pad] as [number, number],
      yMax: Math.max(2, maxY * 1.08),
    };
  }, [points]);

  function openGene(event: any) {
    const payload = event?.payload || event;
    const symbol = String(payload?.gene_symbol || "").trim();
    if (symbol) router.push(`/genes/${encodeURIComponent(symbol)}`);
  }

  const fdrLine = -Math.log10(0.05);

  return <div className="chart-shell">
    <div className="chart-legend">
      <span><i className="dot regular" /> FDR &gt; 0,05</span>
      <span><i style={{background:"#4d7fb7"}} className="dot" /> FDR ≤ 0,05</span>
      <span><i className="dot stable" /> устойчивые повторяющиеся гены</span>
      <span className="muted">Левее нуля — зависимость сильнее в целевой группе</span>
    </div>
    <div style={{width:"100%",height:390}}>
      <ResponsiveContainer>
        <ScatterChart margin={{top:16,right:18,bottom:30,left:12}}>
          <CartesianGrid strokeDasharray="3 4" stroke="#e8edf3" />
          <XAxis
            type="number"
            dataKey="delta_gene_effect"
            domain={xExtent}
            tick={{fontSize:11,fill:"#748398"}}
            axisLine={{stroke:"#d8e0e9"}}
            tickLine={{stroke:"#d8e0e9"}}
            label={{value:"Δ Gene Effect · целевая − группа сравнения",position:"insideBottom",offset:-17,fill:"#63758b",fontSize:11}}
          />
          <YAxis
            type="number"
            dataKey="_significance"
            domain={[0,yMax]}
            tick={{fontSize:11,fill:"#748398"}}
            axisLine={{stroke:"#d8e0e9"}}
            tickLine={{stroke:"#d8e0e9"}}
            label={{value:"−log10(q-value)",angle:-90,position:"insideLeft",fill:"#63758b",fontSize:11}}
          />
          <ReferenceLine x={0} stroke="#9eacbc" strokeDasharray="5 5" />
          <ReferenceLine y={fdrLine} stroke="#9eacbc" strokeDasharray="5 5" />
          <Tooltip content={<CustomTooltip />} cursor={{stroke:"#bdc9d7",strokeDasharray:"3 3"}} />
          <Scatter data={regular} fill="#9bb2cc" fillOpacity={0.36} isAnimationActive={false} onClick={openGene} style={{cursor:"pointer"}} />
          <Scatter data={significant} fill="#4d7fb7" fillOpacity={0.72} isAnimationActive={false} onClick={openGene} style={{cursor:"pointer"}} />
          <Scatter data={stable} fill="#134a7d" fillOpacity={0.95} isAnimationActive={false} onClick={openGene} style={{cursor:"pointer"}} />
        </ScatterChart>
      </ResponsiveContainer>
    </div>
  </div>;
}
