import Link from "next/link";
import { apiGet } from "../../lib/api";
import styles from "../pharmacology.module.css";

export const dynamic = "force-dynamic";

type Params = Record<string,string|string[]|undefined>;
type Catalog = {available:boolean;total:number;limit:number;offset:number;items:any[];note_ru?:string;build_command?:string};
type Summary = Record<string,any>;

function one(value:string|string[]|undefined){return Array.isArray(value)?value[0]:value;}
function yes(value:string|undefined){return value==="1"||value==="true"||value==="on";}
function n(value:any){const x=Number(value);return Number.isFinite(x)?x.toLocaleString("ru-RU"):"—";}
function short(value:any,max=86){const text=String(value||"");return text.length>max?`${text.slice(0,max)}…`:text;}
function pageHref(sp:Params,offset:number){const q=new URLSearchParams();for(const key of ["q","target_gene"]){const v=one(sp[key]);if(v)q.set(key,v);}if(yes(one(sp.has_smiles)))q.set("has_smiles","true");q.set("offset",String(Math.max(0,offset)));q.set("limit",one(sp.limit)||"100");return `/compounds?${q.toString()}`;}

export default async function CompoundsPage({searchParams}:{searchParams:Promise<Params>}){
  const sp=await searchParams;
  const limit=Math.max(1,Math.min(500,Number(one(sp.limit)||"100")||100));
  const offset=Math.max(0,Number(one(sp.offset)||"0")||0);
  const q=new URLSearchParams({limit:String(limit),offset:String(offset)});
  const search=one(sp.q)||"";const target=one(sp.target_gene)||"";const hasSmiles=yes(one(sp.has_smiles));
  if(search)q.set("q",search);if(target)q.set("target_gene",target);if(hasSmiles)q.set("has_smiles","true");
  const [catalog,summary]=await Promise.all([
    apiGet<Catalog>(`/api/compounds?${q.toString()}`),
    apiGet<Summary>("/api/pharmacology/summary"),
  ]);
  const items=catalog.items||[];
  const end=Math.min(catalog.total||0,offset+items.length);
  return <>
    <section className={styles.pageHeader}>
      <div>
        <div className="eyebrow">ФАРМАКОЛОГИЯ · ВЕЩЕСТВО → МОДЕЛЬ → МИШЕНЬ</div>
        <h1>Вещества</h1>
        <p>Единый каталог соединений, протестированных на клеточных моделях. Здесь химическая идентичность, экспериментальный ответ клетки и аннотация мишени хранятся раздельно и связываются только через явные доказательства.</p>
      </div>
      <div className={styles.heroStats}>
        <div className={styles.heroStat}><strong>{n(summary.compounds_n)}</strong><span>compound-профилей</span></div>
        <div className={styles.heroStat}><strong>{n(summary.models_with_response_n)}</strong><span>моделей с фармакологией</span></div>
        <div className={styles.heroStat}><strong>{n(summary.responses_n)}</strong><span>наблюдений</span></div>
      </div>
    </section>

    <form action="/compounds" method="get" className={styles.toolbar}>
      <label className={styles.field}><span>Поиск вещества</span><input name="q" defaultValue={search} placeholder="название, ID, SMILES, ChEMBL…"/></label>
      <label className={styles.field}><span>Белковая мишень / ген</span><input name="target_gene" defaultValue={target} placeholder="EGFR, BRAF, PARP1…"/></label>
      <label className={styles.checkbox}><input type="checkbox" name="has_smiles" value="true" defaultChecked={hasSmiles}/>только с SMILES</label>
      <input type="hidden" name="limit" value={String(limit)}/>
      <button className={styles.button} type="submit">Применить</button>
      <Link className={styles.reset} href="/compounds">Сбросить</Link>
    </form>

    {catalog.note_ru&&<div className={styles.note}>{catalog.note_ru}</div>}

    {!catalog.available?<div className={styles.empty}><b>Каталог веществ ещё не материализован.</b><br/>{catalog.build_command||"Перестройте фармакологический слой."}</div>:
    <>
      <div className={styles.tableWrap}><table className={styles.table}>
        <thead><tr><th>Вещество</th><th>Химическая запись</th><th>Мишени</th><th>Модели</th><th>Наблюдения</th><th>Источники</th></tr></thead>
        <tbody>{items.map((row:any)=>{
          const targets=(row.target_genes||[]) as string[];
          return <tr key={row.compound_id}>
            <td><Link className={styles.primary} href={`/compounds/${encodeURIComponent(row.compound_id)}`}>{row.preferred_name||row.compound_id}</Link><span className={styles.sub}>{row.compound_id}</span>{row.chembl_id&&<span className={styles.sub}>{row.chembl_id}</span>}</td>
            <td>{row.canonical_smiles?<div className={styles.smiles}>{short(row.canonical_smiles)}</div>:<span className={styles.sub}>SMILES не указан</span>}</td>
            <td>{targets.length?<div className={styles.chips}>{targets.slice(0,4).map((gene:string)=><Link className={styles.chip} key={gene} href={`/targets/${encodeURIComponent(gene)}`}>{gene}</Link>)}{targets.length>4&&<span className={styles.chip}>+{targets.length-4}</span>}</div>:<span className={styles.sub}>нет аннотации</span>}</td>
            <td><b>{n(row.models_n)}</b><span className={styles.sub}>клеточных моделей</span></td>
            <td><b>{n(row.observations_n)}</b><span className={styles.sub}>{(row.endpoints||[]).slice(0,3).join(" · ")||"—"}</span></td>
            <td><div className={styles.chips}>{(row.sources||[]).slice(0,3).map((source:string)=><span className={styles.chip} key={source}>{source}</span>)}</div></td>
          </tr>})}</tbody>
      </table></div>
      <div className={styles.pagination}>
        <Link className={offset<=0?styles.disabled:""} href={pageHref(sp,Math.max(0,offset-limit))}>← Назад</Link>
        <span className={styles.sub}>{catalog.total?`${offset+1}–${end} из ${n(catalog.total)}`:"0 записей"}</span>
        <Link className={end>=catalog.total?styles.disabled:""} href={pageHref(sp,offset+limit)}>Дальше →</Link>
      </div>
    </>}
  </>;
}
