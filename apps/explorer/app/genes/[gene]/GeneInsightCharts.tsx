"use client";

import Link from "next/link";
import styles from "./GeneInsightCharts.module.css";

type Point = {
  model_id?: string;
  cell_line_name?: string;
  gene_effect?: number | null;
  expression?: number | null;
  copy_number?: number | null;
  oncotree_subtype?: string | null;
  mcl_cancer_name?: string | null;
  mcl_cancer_id?: string | null;
};

type CancerRelationship={cancer_id?:string;cancer_name?:string;organ_ru?:string;rho?:number|null;p_value?:number|null;n?:number;strength?:string};
type Relationship = {
  layer?: "expression" | "copy_number" | string;
  label?: string;
  available?: boolean;
  rho?: number | null;
  p_value?: number | null;
  n?: number;
  interpretation?:{strength?:string;direction?:string;statistically_clear?:boolean;conclusion_ru?:string};
  within_cancers?:CancerRelationship[];
};

function fmt(value: unknown, digits = 2) {
  const n = Number(value);
  return Number.isFinite(n) ? n.toFixed(digits) : "—";
}
function pFmt(value:unknown){const n=Number(value);if(!Number.isFinite(n))return "—";if(n<.0001)return "< 0,0001";return n.toLocaleString("ru-RU",{maximumFractionDigits:4});}
function quantile(values:number[],q:number){if(!values.length)return 0;const a=[...values].sort((x,y)=>x-y);const pos=(a.length-1)*q;const lo=Math.floor(pos),hi=Math.ceil(pos);return lo===hi?a[lo]:a[lo]+(a[hi]-a[lo])*(pos-lo);}
function linearFit(xs:number[],ys:number[]){const n=xs.length;if(n<2)return null;const mx=xs.reduce((a,b)=>a+b,0)/n,my=ys.reduce((a,b)=>a+b,0)/n;let num=0,den=0;for(let i=0;i<n;i++){num+=(xs[i]-mx)*(ys[i]-my);den+=(xs[i]-mx)*(xs[i]-mx);}if(!den)return null;const slope=num/den;return {slope,intercept:my-slope*mx};}

function Scatter({layer,title,points,relationship}:{layer:"expression"|"copy_number";title:string;points:Point[];relationship?:Relationship}){
  const all=points.filter((p)=>Number.isFinite(Number(p[layer]))&&Number.isFinite(Number(p.gene_effect)));
  if(!all.length) return <div className={styles.empty}>Слой {title} для этого гена не доступен.</div>;
  const allX=all.map((p)=>Number(p[layer]));
  const clipLow=layer==="copy_number"?quantile(allX,.01):Math.min(...allX);
  const clipHigh=layer==="copy_number"?quantile(allX,.99):Math.max(...allX);
  const usable=all.filter((p)=>Number(p[layer])>=clipLow&&Number(p[layer])<=clipHigh);
  const xs=usable.map((p)=>Number(p[layer]));
  const ys=usable.map((p)=>Number(p.gene_effect));
  const xmin=Math.min(...xs), xmax=Math.max(...xs), ymin=Math.min(...ys,-1), ymax=Math.max(...ys,0);
  const dx=(xmax-xmin)||1, dy=(ymax-ymin)||1;
  const left=52,right=16,top=18,bottom=36,width=520,height=286;
  const sx=(x:number)=>left+((x-xmin)/dx)*(width-left-right);
  const sy=(y:number)=>top+(1-(y-ymin)/dy)*(height-top-bottom);
  const fit=linearFit(xs,ys);
  const clipped=all.length-usable.length;
  const lines=[{y:0,label:"нет выраженной зависимости"},{y:-.5,label:"порог MCL −0,5"},{y:-1,label:"сильная зависимость"}].filter(x=>x.y>=ymin&&x.y<=ymax);
  return <article className={styles.chartCard}>
    <div className={styles.chartHead}>
      <div><b>{title} → CRISPR-зависимость</b><span>Spearman ρ {fmt(relationship?.rho,3)} · p {pFmt(relationship?.p_value)} · n={relationship?.n??all.length}</span></div>
      <span className={styles.badge}>{relationship?.interpretation?.strength || "исследовательская связь"}</span>
    </div>
    <div className={styles.interpretation}>
      <strong>{relationship?.interpretation?.direction || "Направление связи не определено"}</strong>
      <p>{relationship?.interpretation?.conclusion_ru || "Корреляция оценивает совместное изменение показателей, но не доказывает причинность."}</p>
    </div>
    <svg className={styles.svg} viewBox={`0 0 ${width} ${height}`} role="img" aria-label={`${title} versus Gene Effect`}>
      <line x1={left} x2={width-right} y1={height-bottom} y2={height-bottom} className={styles.axis}/>
      <line x1={left} x2={left} y1={top} y2={height-bottom} className={styles.axis}/>
      {lines.map((line)=><g key={line.y}><line x1={left} x2={width-right} y1={sy(line.y)} y2={sy(line.y)} className={line.y===-.5?styles.threshold:styles.guide}/><text x={width-right-2} y={sy(line.y)-4} textAnchor="end" className={styles.guideLabel}>{line.label}</text></g>)}
      <text x={left} y={height-8} className={styles.tick}>{fmt(xmin)}</text>
      <text x={width-right} y={height-8} textAnchor="end" className={styles.tick}>{fmt(xmax)}</text>
      <text x={7} y={top+4} className={styles.tick}>{fmt(ymax)}</text>
      <text x={7} y={height-bottom} className={styles.tick}>{fmt(ymin)}</text>
      <text x={(left+width-right)/2} y={height-8} textAnchor="middle" className={styles.label}>{layer==="expression"?"RNA · log2(TPM+1)":"Относительное число копий"}</text>
      <text x={12} y={(top+height-bottom)/2} transform={`rotate(-90 12 ${(top+height-bottom)/2})`} textAnchor="middle" className={styles.label}>Gene Effect</text>
      {fit&&<line x1={sx(xmin)} y1={sy(fit.intercept+fit.slope*xmin)} x2={sx(xmax)} y2={sy(fit.intercept+fit.slope*xmax)} className={styles.trend}/>}      
      {usable.map((p,index)=>{
        const x=Number(p[layer]), y=Number(p.gene_effect);
        const titleText=`${p.cell_line_name||p.model_id||"model"} · ${p.mcl_cancer_name||p.oncotree_subtype||""} · ${layer}=${fmt(x,3)} · GE=${fmt(y,3)}`;
        return <circle key={`${p.model_id}-${index}`} cx={sx(x)} cy={sy(y)} r="3.1" className={styles.point}><title>{titleText}</title></circle>;
      })}
    </svg>
    <div className={styles.chartFoot}>{clipped>0?`Для читаемости скрыто ${clipped} крайних значений по оси X; расчёт ρ использует все ${all.length} моделей.`:`Показаны все ${all.length} моделей с доступными парными измерениями.`}</div>
    {!!relationship?.within_cancers?.length&&<div className={styles.cancerLinks}><span>Где связь выражена сильнее:</span>{relationship.within_cancers.slice(0,4).map((row)=><div key={row.cancer_id}><b>{row.cancer_name||row.cancer_id}</b><small>{row.organ_ru||""} · ρ {fmt(row.rho,2)} · n={row.n}</small></div>)}</div>}
  </article>;
}

export default function GeneInsightCharts({points,relationships}:{points:Point[];relationships:Relationship[]}){
  const expression=relationships.find((x)=>x.layer==="expression");
  const copyNumber=relationships.find((x)=>x.layer==="copy_number");
  const available=[expression,copyNumber].filter(Boolean) as Relationship[];
  const weak=available.every((x)=>Math.abs(Number(x.rho||0))<.4);
  const directions=available.map(x=>Number(x.rho||0)<0?"stronger":"weaker");
  const sameDirection=directions.length>1&&directions.every(x=>x===directions[0]);
  return <>
    <div className={styles.overall}>
      <div><span>Автоматическая интерпретация</span><strong>{weak?"RNA и число копий объясняют лишь часть различий между моделями":"Есть выраженная связь хотя бы с одним молекулярным слоем"}</strong></div>
      <p>{sameDirection&&directions[0]==="stronger"?"В обоих слоях более высокие значения связаны с более сильной CRISPR-зависимостью. ":sameDirection?"В обоих слоях более высокие значения связаны с более слабой зависимостью. ":"RNA и число копий показывают разное направление связи. "}Это описательная ассоциация; тип опухоли, мутации и другие факторы могут объяснять часть наблюдаемого паттерна.</p>
    </div>
    <div className={styles.grid}>
      <Scatter layer="expression" title="Экспрессия RNA" points={points} relationship={expression}/>
      <Scatter layer="copy_number" title="Число копий" points={points} relationship={copyNumber}/>
    </div>
  </>;
}
