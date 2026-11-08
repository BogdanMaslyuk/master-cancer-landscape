import Link from "next/link";
import { apiGet } from "../../../lib/api";
import styles from "../pyz.module.css";

export const dynamic = "force-dynamic";

type Params=Record<string,string|string[]|undefined>;
function one(v:string|string[]|undefined){return Array.isArray(v)?v[0]:v;}
function cellClass(value:number){if(value>=.50)return styles.sim4;if(value>=.35)return styles.sim3;if(value>=.20)return styles.sim2;if(value>0)return styles.sim1;return styles.sim0;}

export default async function PyzMatrixPage({searchParams}:{searchParams:Promise<Params>}){
  const sp=await searchParams;
  const q=one(sp.q)||"";
  const shortlist=one(sp.shortlist_only)==="true"||one(sp.shortlist_only)==="1";
  const params=new URLSearchParams();
  if(q)params.set("q",q);
  if(shortlist)params.set("shortlist_only","true");
  const data=await apiGet<Record<string,any>>(`/api/pyz/matrix?${params.toString()}`);
  const targets=(data.targets||[]) as string[];
  const compounds=(data.compounds||[]) as any[];
  const items=(data.items||[]) as any[];
  const byKey=new Map<string,any>();
  for(const row of items)byKey.set(`${row.own_compound_id}|${row.target_gene}`,row);
  const visibleIds=q?Array.from(new Set(items.map((x:any)=>String(x.own_compound_id)))):compounds.map((x:any)=>String(x.own_compound_id));
  const visible=compounds.filter((x:any)=>visibleIds.includes(String(x.own_compound_id)));

  return <>
    <Link href="/pyz" className={styles.backLink}>← Молекулы PYZ</Link>
    <section className={styles.pageHeader}>
      <div><div className="eyebrow">65 МОЛЕКУЛ × 13 CORE-МИШЕНЕЙ</div><h1>Матрица PYZ × мишень</h1><p>Каждая ячейка показывает лучший Morgan/Tanimoto между нашей молекулой и экспериментально измеренным лигандом конкретной мишени. Матрица нужна для навигации и отбора гипотез, а не для утверждения механизма.</p></div>
      <div className={styles.heroStats}><div className={styles.heroStat}><strong>845</strong><span>полных сравнений</span></div><div className={styles.heroStat}><strong>10</strong><span>пар ≥0,35</span></div><div className={styles.heroStat}><strong>0</strong><span>пар ≥0,50</span></div><div className={styles.heroStat}><strong>CDK4 · BCL2L1</strong><span>мишени shortlist</span></div></div>
    </section>

    <div className={styles.warning}><b>Как читать:</b> значение <b>&lt;0,35</b> не исключает связывание другим хемотипом, но не позволяет переносить механизм ближайшего известного лиганда. Значения 0,35–0,50 — только исследовательский shortlist.</div>

    <form action="/pyz/matrix" method="get" className={styles.toolbar} style={{marginTop:16}}>
      <label className={styles.field}><span>Поиск PYZ или мишени</span><input name="q" defaultValue={q} placeholder="PYZ-021 или CDK4"/></label>
      <label className={styles.field}><span>Режим</span><select name="shortlist_only" defaultValue={shortlist?"true":"false"}><option value="false">полная матрица</option><option value="true">только ячейки ≥0,35</option></select></label>
      <button className={styles.button} type="submit">Применить</button><Link href="/pyz/matrix" className={styles.reset}>Сбросить</Link>
    </form>

    <div className={styles.legend}><span className={styles.sim0}>почти нет сходства</span><span className={styles.sim1}>слабое</span><span className={styles.sim2}>0,20–0,35</span><span className={styles.sim3}>0,35–0,50 · shortlist</span><span className={styles.sim4}>≥0,50</span></div>

    <div className={styles.matrixWrap}><table className={styles.matrix}><thead><tr><th>PYZ</th>{targets.map(g=><th key={g}>{g}</th>)}</tr></thead><tbody>{visible.map((compound:any)=><tr key={compound.own_compound_id}><td><Link href={`/pyz/${encodeURIComponent(compound.own_compound_id)}`}>{compound.own_compound_id}</Link><div className={styles.muted}>{compound.chemotype_ru}</div></td>{targets.map(g=>{const row=byKey.get(`${compound.own_compound_id}|${g}`);if(!row)return <td key={g}><span className={`${styles.matrixCell} ${styles.sim0}`}>—</span></td>;const value=Number(row.best_tanimoto_morgan_r2_2048);return <td key={g}><Link title={`${compound.own_compound_id} → ${g}: Tanimoto ${value.toFixed(3)}`} href={`/pyz/${encodeURIComponent(compound.own_compound_id)}?target_gene=${encodeURIComponent(g)}`} className={`${styles.matrixCell} ${cellClass(value)}`}>{value.toFixed(3)}</Link></td>})}</tr>)}</tbody></table></div>
    <div className={styles.notice} style={{marginTop:12}}>{data.interpretation_ru}</div>
  </>;
}
