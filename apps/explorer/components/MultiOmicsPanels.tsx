import Link from "next/link";
import styles from "./MultiOmics.module.css";

type Payload = Record<string, any>;

function fmt(value: any, digits = 2){
  if(value === null || value === undefined || value === "") return "—";
  const n = Number(value);
  if(!Number.isFinite(n)) return String(value);
  return n.toFixed(digits).replace(/\.00$/, "").replace(/(\.\d)0$/, "$1");
}

function layerRu(layer:string){
  return ({
    expression:"RNA expression",
    copy_number:"Copy number",
    gene_effect:"CRISPR Gene Effect",
  } as Record<string,string>)[layer] || layer;
}

function LayerCards({availability, modelLayers}:{availability:any; modelLayers?:any}){
  const layers=availability?.layers || {};
  return <div className={styles.layerGrid}>
    {["expression","copy_number","gene_effect"].map((layer)=>{
      const meta=layers[layer] || {};
      const modelState=modelLayers?.[layer];
      const available=Boolean(meta.available) && (modelState ? Boolean(modelState.model_present) : true);
      return <div className={styles.layerCard} key={layer}>
        <span className={`${styles.layerState} ${available ? styles.ready : styles.missing}`}>{available ? "доступно" : "не подключено"}</span>
        <strong>{layerRu(layer)}</strong>
        <span>{meta.value_semantics || (layer === "expression" ? "log2(TPM + 1)" : layer === "copy_number" ? "relative linear" : "Chronos Gene Effect")}</span>
      </div>;
    })}
  </div>;
}

export function ModelMultiOmicsPanel({omics}:{omics:Payload}){
  const availability=omics?.availability || {};
  const panel=omics?.candidate_panel || [];
  const top=omics?.top_dependencies || [];

  return <div className={styles.panel}>
    <LayerCards availability={availability} modelLayers={omics?.model_layers}/>
    <div className={styles.content}>
      {!availability.available ? <div className={styles.empty}>
        <strong>Multi-omics индекс ещё не построен</strong>
        <p>После индексации здесь появятся RNA expression, относительное число копий и индивидуальный CRISPR Gene Effect для этой модели.</p>
        <code>{availability.build_command || ".\\.venv\\Scripts\\python.exe .\\scripts\\build_depmap_multiomics.py"}</code>
      </div> : <>
        <p className={styles.note}>Один и тот же ген здесь читается сразу в трёх слоях: насколько он нужен клетке после CRISPR-выключения, насколько выражена его РНК и как выглядит относительное число копий. Эти слои помогают объяснять зависимость, но сами по себе не доказывают лекарственную уязвимость.</p>

        {panel.length > 0 && <>
          <div className={styles.subheader}><h3>Устойчивые кандидаты MCL в этой модели</h3><span>единая шкала не используется — каждый столбец интерпретируется отдельно</span></div>
          <div className={styles.tableWrap}>
            <table className={styles.table}>
              <thead><tr><th>Ген</th><th>Gene Effect</th><th>RNA expression</th><th>Relative copy number</th></tr></thead>
              <tbody>{panel.map((row:any)=><tr key={row.gene}>
                <td><Link href={`/genes/${encodeURIComponent(row.gene)}`} className={styles.geneLink}>{row.gene}</Link></td>
                <td className={`${styles.numeric} ${Number(row.gene_effect) <= -0.5 ? styles.strongDependency : ""}`}>{fmt(row.gene_effect,3)}</td>
                <td className={styles.numeric}>{fmt(row.expression,2)}</td>
                <td className={styles.numeric}>{fmt(row.copy_number,2)}</td>
              </tr>)}</tbody>
            </table>
          </div>
        </>}

        {top.length > 0 && <>
          <div className={styles.subheader}><h3>Сильнейшие индивидуальные CRISPR-зависимости</h3><span>наиболее отрицательный Gene Effect в этой клеточной линии</span></div>
          <div className={styles.tableWrap}>
            <table className={styles.table}>
              <thead><tr><th>Ген</th><th>Gene Effect</th><th>RNA expression</th><th>Relative copy number</th></tr></thead>
              <tbody>{top.slice(0,20).map((row:any)=><tr key={row.gene}>
                <td><Link href={`/genes/${encodeURIComponent(row.gene)}`} className={styles.geneLink}>{row.gene}</Link></td>
                <td className={`${styles.numeric} ${Number(row.gene_effect) <= -0.5 ? styles.strongDependency : ""}`}>{fmt(row.gene_effect,3)}</td>
                <td className={styles.numeric}>{fmt(row.expression,2)}</td>
                <td className={styles.numeric}>{fmt(row.copy_number,2)}</td>
              </tr>)}</tbody>
            </table>
          </div>
        </>}

        <div className={styles.interpretation}>
          <div><b>Gene Effect</b><span>{omics?.interpretation?.gene_effect || "Более отрицательное значение означает более сильную CRISPR-зависимость."}</span></div>
          <div><b>RNA expression</b><span>{omics?.interpretation?.expression || "Экспрессия РНК не равна активности белка."}</span></div>
          <div><b>Copy number</b><span>{omics?.interpretation?.copy_number || "Относительное число копий не трактуется как абсолютный клинический CNV-вызов."}</span></div>
        </div>
      </>}
    </div>
  </div>;
}

export function ContextMultiOmicsPanel({omics}:{omics:Payload}){
  const availability=omics?.availability || {};
  const coverage=omics?.coverage || {};
  const panel=omics?.candidate_panel || [];

  return <div className={styles.panel}>
    <LayerCards availability={availability}/>
    <div className={styles.content}>
      {!availability.available ? <div className={styles.empty}>
        <strong>Multi-omics слой подготовлен архитектурно, но локальные индексы ещё не построены</strong>
        <p>Нужны матрицы RNA expression и copy number того же релиза DepMap; CRISPR Gene Effect уже используется в функциональном анализе и также войдёт в единый индекс.</p>
        <code>{availability.build_command || ".\\.venv\\Scripts\\python.exe .\\scripts\\build_depmap_multiomics.py"}</code>
      </div> : <>
        <div className={styles.coverageGrid}>
          {["expression","copy_number","gene_effect"].map((layer)=>{
            const row=coverage[layer] || {};
            return <div className={styles.coverageCard} key={layer}>
              <strong>{layerRu(layer)}</strong>
              <b>{row.models_n ?? 0}/{row.models_total_n ?? 0}</b>
              <span>моделей контекста имеют этот слой · целевая {row.context_models_n ?? 0}/{row.context_total_n ?? 0} · контроль {row.comparator_models_n ?? 0}/{row.comparator_total_n ?? 0}</span>
            </div>;
          })}
        </div>

        {panel.length > 0 && <>
          <div className={styles.subheader}><h3>Multi-omics проверка устойчивых кандидатов</h3><span>медианы целевой группы и группы сравнения</span></div>
          <div className={styles.tableWrap}>
            <table className={styles.table}>
              <thead><tr><th>Ген</th><th>Δ Gene Effect</th><th>Gene Effect target / control</th><th>RNA target / control</th><th>CN target / control</th></tr></thead>
              <tbody>{panel.map((row:any)=><tr key={row.gene}>
                <td><Link href={`/genes/${encodeURIComponent(row.gene)}`} className={styles.geneLink}>{row.gene}</Link></td>
                <td className={`${styles.numeric} ${Number(row.gene_effect_delta) < 0 ? styles.strongDependency : ""}`}>{fmt(row.gene_effect_delta,3)}</td>
                <td className={styles.numeric}>{fmt(row.gene_effect_context_median,3)} / {fmt(row.gene_effect_comparator_median,3)}</td>
                <td className={styles.numeric}>{fmt(row.expression_context_median,2)} / {fmt(row.expression_comparator_median,2)}</td>
                <td className={styles.numeric}>{fmt(row.copy_number_context_median,2)} / {fmt(row.copy_number_comparator_median,2)}</td>
              </tr>)}</tbody>
            </table>
          </div>
        </>}
        <p className={styles.note} style={{marginTop:14,marginBottom:0}}>{omics.note}</p>
      </>}
    </div>
  </div>;
}
