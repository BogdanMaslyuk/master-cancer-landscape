"use client";

import styles from "./GeneInsightCharts.module.css";

type Point = {
  model_id?: string;
  cell_line_name?: string;
  gene_effect?: number | null;
  expression?: number | null;
  copy_number?: number | null;
  oncotree_subtype?: string | null;
  cancer_ids?: string[];
};

type Relationship = {
  layer: "expression" | "copy_number" | string;
  label: string;
  available: boolean;
  rho?: number | null;
  p_value?: number | null;
  n?: number;
};

function fmt(value: unknown, digits = 2) {
  const n = Number(value);
  return Number.isFinite(n) ? n.toFixed(digits) : "—";
}

function Scatter({layer, title, points, relationship}:{layer:"expression"|"copy_number";title:string;points:Point[];relationship?:Relationship}){
  const usable=points.filter((p)=>Number.isFinite(Number(p[layer]))&&Number.isFinite(Number(p.gene_effect)));
  if(!usable.length) return <div className={styles.empty}>Слой {title} для этого гена не доступен.</div>;
  const xs=usable.map((p)=>Number(p[layer]));
  const ys=usable.map((p)=>Number(p.gene_effect));
  const xmin=Math.min(...xs), xmax=Math.max(...xs), ymin=Math.min(...ys), ymax=Math.max(...ys);
  const dx=(xmax-xmin)||1, dy=(ymax-ymin)||1;
  const left=48,right=16,top=18,bottom=34,width=520,height=280;
  const sx=(x:number)=>left+((x-xmin)/dx)*(width-left-right);
  const sy=(y:number)=>top+(1-(y-ymin)/dy)*(height-top-bottom);
  const zeroY=(ymin<=0&&ymax>=0)?sy(0):null;
  return <article className={styles.chartCard}>
    <div className={styles.chartHead}>
      <div><b>{title} ↔ Gene Effect</b><span>Spearman ρ {fmt(relationship?.rho,3)} · p {fmt(relationship?.p_value,4)} · n={relationship?.n??usable.length}</span></div>
      <span className={styles.badge}>исследовательская корреляция</span>
    </div>
    <svg className={styles.svg} viewBox={`0 0 ${width} ${height}`} role="img" aria-label={`${title} versus Gene Effect`}>
      <line x1={left} x2={width-right} y1={height-bottom} y2={height-bottom} className={styles.axis}/>
      <line x1={left} x2={left} y1={top} y2={height-bottom} className={styles.axis}/>
      {zeroY!==null&&<line x1={left} x2={width-right} y1={zeroY} y2={zeroY} className={styles.zero}/>}      
      <text x={left} y={height-8} className={styles.tick}>{fmt(xmin)}</text>
      <text x={width-right} y={height-8} textAnchor="end" className={styles.tick}>{fmt(xmax)}</text>
      <text x={7} y={top+4} className={styles.tick}>{fmt(ymax)}</text>
      <text x={7} y={height-bottom} className={styles.tick}>{fmt(ymin)}</text>
      <text x={(left+width-right)/2} y={height-8} textAnchor="middle" className={styles.label}>{layer==="expression"?"RNA · log2(TPM+1)":"Relative copy number"}</text>
      <text x={12} y={(top+height-bottom)/2} transform={`rotate(-90 12 ${(top+height-bottom)/2})`} textAnchor="middle" className={styles.label}>Gene Effect</text>
      {usable.map((p,index)=>{
        const x=Number(p[layer]), y=Number(p.gene_effect);
        const titleText=`${p.cell_line_name||p.model_id||"model"} · ${p.oncotree_subtype||""} · ${layer}=${fmt(x,3)} · GE=${fmt(y,3)}`;
        return <circle key={`${p.model_id}-${index}`} cx={sx(x)} cy={sy(y)} r="4.2" className={styles.point}><title>{titleText}</title></circle>;
      })}
    </svg>
    <p>Каждая точка — клеточная модель. Более отрицательный Gene Effect означает более сильную зависимость после CRISPR-выключения.</p>
  </article>;
}

export default function GeneInsightCharts({points,relationships}:{points:Point[];relationships:Relationship[]}){
  const expression=relationships.find((x)=>x.layer==="expression");
  const copyNumber=relationships.find((x)=>x.layer==="copy_number");
  return <div className={styles.grid}>
    <Scatter layer="expression" title="RNA expression" points={points} relationship={expression}/>
    <Scatter layer="copy_number" title="Copy number" points={points} relationship={copyNumber}/>
  </div>;
}
