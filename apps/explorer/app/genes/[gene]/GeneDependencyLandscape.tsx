import Link from "next/link";
import styles from "./GeneDependencyLandscape.module.css";

type CancerRow = {
  cancer_id?: string;
  cancer_name?: string;
  organ_ru?: string;
  system_ru?: string;
  models_n?: number;
  dependent_models_n?: number;
  dependent_fraction?: number;
  median_gene_effect?: number | null;
  q1_gene_effect?: number | null;
  q3_gene_effect?: number | null;
  specificity_vs_global?: number | null;
  eligible_for_label?: boolean;
};

type Payload = {
  available?: boolean;
  summary?: Record<string, any>;
  cancers?: CancerRow[];
  interpretation?: Record<string, string>;
};

function pct(value: unknown, digits = 0) {
  const n = Number(value);
  return Number.isFinite(n) ? `${(n * 100).toFixed(digits)}%` : "—";
}
function num(value: unknown, digits = 2) {
  const n = Number(value);
  return Number.isFinite(n) ? n.toFixed(digits) : "—";
}
function specificity(row: CancerRow) {
  const value = Number(row.specificity_vs_global);
  if (!Number.isFinite(value)) return "не определена";
  if (value >= .5) return "очень высокая";
  if (value >= .3) return "высокая";
  if (value >= .15) return "средняя";
  return "низкая";
}

export default function GeneDependencyLandscape({symbol,payload}:{symbol:string;payload:Payload}){
  if(!payload?.available) return <div className={styles.empty}>Полный профиль CRISPR-зависимости для {symbol} пока недоступен.</div>;
  const summary=payload.summary||{};
  const rows=(payload.cancers||[]).slice().sort((a,b)=>Number(b.dependent_fraction||0)-Number(a.dependent_fraction||0));
  const visible=rows.slice(0,24);
  return <>
    <div className={styles.summaryGrid}>
      <div className={styles.summaryCard}><span>Зависимые модели</span><strong>{summary.dependent_models_n ?? 0} / {summary.models_n ?? 0}</strong><small>{pct(summary.dependent_fraction,1)} при Gene Effect ≤ −0,5</small></div>
      <div className={styles.summaryCard}><span>Тип зависимости</span><strong className={styles.compact}>{summary.dependency_type_ru || "не классифицирована"}</strong><small>описательный CRISPR-профиль</small></div>
      <div className={styles.summaryCard}><span>Наиболее выражена</span><strong className={styles.compact}>{summary.best_cancer_name || "—"}</strong><small>{summary.best_cancer_organ_ru || ""} · {pct(summary.best_cancer_dependency_fraction,0)} моделей</small></div>
      <div className={styles.summaryCard}><span>Специфичность</span><strong>{summary.specificity_label_ru || "—"}</strong><small>оценка по распределению 1208 моделей</small></div>
    </div>

    <div className={styles.explainer}>
      <b>Как читать:</b> высокая доля означает, что во многих моделях этой опухоли выключение {symbol} сопровождается выраженной потерей жизнеспособности. Специфичность показывает, насколько эта доля выше общего фона среди всех CRISPR-моделей.
    </div>

    <div className={styles.tableWrap}>
      <table className={styles.table}>
        <thead><tr><th>Опухоль</th><th>Моделей</th><th>Зависимых</th><th>Доля</th><th>Медиана GE</th><th>Специфичность</th></tr></thead>
        <tbody>{visible.map((row,index)=>{
          const fraction=Math.max(0,Math.min(1,Number(row.dependent_fraction||0)));
          return <tr key={`${row.cancer_id}-${index}`}>
            <td><div className={styles.cancerName}>{row.cancer_name || row.cancer_id || "—"}</div><span>{row.organ_ru || row.system_ru || "—"}</span></td>
            <td>{row.models_n ?? 0}</td>
            <td><b>{row.dependent_models_n ?? 0}</b></td>
            <td><div className={styles.fraction}><b>{pct(fraction,0)}</b><span><i style={{width:`${fraction*100}%`}}/></span></div></td>
            <td className={Number(row.median_gene_effect)<0?styles.negative:""}>{num(row.median_gene_effect)}</td>
            <td><span className={`${styles.spec} ${!row.eligible_for_label?styles.lowN:""}`}>{row.eligible_for_label?specificity(row):"малое n"}</span>{row.eligible_for_label&&<small>{Number(row.specificity_vs_global)>=0?"+":""}{pct(row.specificity_vs_global,0)} к общему фону</small>}</td>
          </tr>;
        })}</tbody>
      </table>
    </div>
    {rows.length>visible.length&&<div className={styles.note}>Показаны 24 опухолевые группы с наибольшей долей зависимых моделей из {rows.length}. Полный список остаётся доступен через API и будет подключён к интерактивной фильтрации.</div>}
    <div className={styles.foot}>{payload.interpretation?.dependency_definition} {payload.interpretation?.minimum_group} {payload.interpretation?.guardrail}</div>
  </>;
}
