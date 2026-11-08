import Link from "next/link";
import { apiGet } from "../../lib/api";
import styles from "../pharmacology.module.css";

export const dynamic = "force-dynamic";

type Params = Record<string,string|string[]|undefined>;
type Catalog = {available:boolean;total:number;limit:number;offset:number;items:any[];note_ru?:string;build_command?:string;protein_registry_available?:boolean};
type Summary = Record<string,any>;
function one(value:string|string[]|undefined){return Array.isArray(value)?value[0]:value;}
function n(value:any){const x=Number(value);return Number.isFinite(x)?x.toLocaleString("ru-RU"):"—";}
function pct(value:any){if(value===null||value===undefined||value==="")return "—";const x=Number(value);return Number.isFinite(x)?`${Math.round(x*100)}%`:"—";}
function pageHref(sp:Params,offset:number){const q=new URLSearchParams();const v=one(sp.q);if(v)q.set("q",v);q.set("offset",String(Math.max(0,offset)));q.set("limit",one(sp.limit)||"100");return `/targets?${q.toString()}`;}

export default async function TargetsPage({searchParams}:{searchParams:Promise<Params>}){
  const sp=await searchParams;
  const limit=Math.max(1,Math.min(500,Number(one(sp.limit)||"100")||100));
  const offset=Math.max(0,Number(one(sp.offset)||"0")||0);
  const search=one(sp.q)||"";
  const query=new URLSearchParams({limit:String(limit),offset:String(offset)});if(search)query.set("q",search);
  const [catalog,summary]=await Promise.all([
    apiGet<Catalog>(`/api/targets?${query.toString()}`),
    apiGet<Summary>("/api/pharmacology/summary"),
  ]);
  const items=catalog.items||[];const end=Math.min(catalog.total||0,offset+items.length);
  return <>
    <section className={styles.pageHeader}>
      <div>
        <div className="eyebrow">БЕЛОК ↔ КОДИРУЮЩИЙ ГЕН ↔ ВЕЩЕСТВО ↔ КЛЕТОЧНАЯ МОДЕЛЬ</div>
        <h1>Белковые мишени</h1>
        <p>Белок и ген показаны как разные сущности. Фармакологический источник обычно задаёт мишень через символ гена, а MCL отдельно сопоставляет его с reviewed-записью UniProtKB/Swiss-Prot. Такое сопоставление не означает, что исходный эксперимент различал конкретную изоформу.</p>
      </div>
      <div className={styles.heroStats}>
        <div className={styles.heroStat}><strong>{n(summary.targets_n)}</strong><span>геновых аннотаций мишеней</span></div>
        <div className={styles.heroStat}><strong>{n(summary.protein_target_registry_n)}</strong><span>записей в реестре белков</span></div>
        <div className={styles.heroStat}><strong>{n(summary.compound_target_pairs_n)}</strong><span>пар вещество × мишень</span></div>
      </div>
    </section>

    <form action="/targets" method="get" className={styles.toolbar}>
      <label className={styles.field}><span>Поиск белка, гена или UniProt</span><input name="q" defaultValue={search} placeholder="5-hydroxytryptamine receptor 2A, HTR2A, P28223…"/></label>
      <input type="hidden" name="limit" value={String(limit)}/>
      <button className={styles.button} type="submit">Найти</button>
      <Link className={styles.reset} href="/targets">Сбросить</Link>
    </form>
    {catalog.note_ru&&<div className={styles.note}>{catalog.note_ru}</div>}

    {!catalog.available?<div className={styles.empty}><b>Каталог мишеней ещё не материализован.</b><br/>{catalog.build_command||"Перестройте фармакологический слой."}</div>:
    <>
      <div className={styles.tableWrap}><table className={styles.table}>
        <thead><tr><th>Белковая мишень</th><th>Кодирующий ген</th><th>Веществ</th><th>Моделей</th><th>CRISPR-зависимость гена</th><th>Аннотации действия</th></tr></thead>
        <tbody>{items.map((row:any)=>{const protein=row.protein_preferred_name||null;return <tr key={row.target_gene}>
          <td><Link className={styles.primary} href={`/targets/${encodeURIComponent(row.target_gene)}`}>{protein||`Белковая мишень, связанная с ${row.target_gene}`}</Link><span className={styles.sub}>{row.uniprot_primary_accession?`UniProt ${row.uniprot_primary_accession}`:row.protein_mapping_status==="multiple_swissprot"?"несколько reviewed-белковых записей":"белковое сопоставление не разрешено"}</span></td>
          <td><Link className={styles.primary} href={`/genes/${encodeURIComponent(row.target_gene)}`}>{row.target_gene} →</Link><span className={styles.sub}>символ гена из фармакологической аннотации</span></td>
          <td><b>{n(row.compounds_n)}</b><span className={styles.sub}>{n(row.target_evidence_n)} строк доказательств</span></td>
          <td><b>{n(row.models_n)}</b><span className={styles.sub}>с фармакологическим покрытием</span></td>
          <td><b>{pct(row.dependency_fraction)}</b><span className={styles.sub}>{n(row.dependent_models_n)} / {n(row.crispr_models_n)} моделей с GE ≤ −0,5</span></td>
          <td>{(row.actions||[]).length?<div className={styles.chips}>{row.actions.slice(0,4).map((x:string)=><span className={styles.chip} key={x}>{x}</span>)}</div>:<span className={styles.sub}>действие не уточнено</span>}</td>
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
