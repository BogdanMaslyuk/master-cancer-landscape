"use client";

import { useEffect, useMemo, useState } from "react";
import { apiBase } from "../../lib/api";
import styles from "./phenotype.module.css";

type Model = {cell_line_id: string; display_name: string; model_role: string; inventory_status: string; assay_ready: boolean; note_ru: string};
type Compound = {compound_id: string};
type Flag = {severity: string; code: string; scope?: string; well?: string; compound_id?: string};
type Imported = {sha256: string; filename: string; wells_n: number; reused_existing_import: boolean; qc: {status: string; flags: Flag[]}};
type Experiment = {experiment_id: string; title: string; experiment_status: string; cell_line_ids: string[]; compound_ids: string[]; imported_files_n: number; raw_wells_n: number; selectivity_status: string};

async function getJson<T>(path: string, init?: RequestInit): Promise<T> {
  const response = await fetch(apiBase + path, {...init, cache: "no-store"});
  if (!response.ok) {
    const raw = await response.text();
    let detail = raw;
    try {const parsed = JSON.parse(raw); if (typeof parsed.detail === "string") detail = parsed.detail;} catch {}
    throw new Error(detail);
  }
  return response.json() as Promise<T>;
}

export default function PhenotypePage() {
  const [models, setModels] = useState<Model[]>([]);
  const [compounds, setCompounds] = useState<Compound[]>([]);
  const [variant, setVariant] = useState("HEK293T");
  const [mrc, setMrc] = useState(true);
  const [bj, setBj] = useState(false);
  const [selected, setSelected] = useState<string[]>(["PYZ-035","PYZ-055","PYZ-031","PYZ-036","PYZ-021"]);
  const [search, setSearch] = useState("");
  const [title, setTitle] = useState("Первичный скрининг PYZ");
  const [protocol, setProtocol] = useState("DRAFT");
  const [hours, setHours] = useState(72);
  const [experimentId, setExperimentId] = useState("");
  const [experiment, setExperiment] = useState<Experiment | null>(null);
  const [file, setFile] = useState<File | null>(null);
  const [lastImport, setLastImport] = useState<Imported | null>(null);
  const [busy, setBusy] = useState(false);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState("");
  const [notice, setNotice] = useState("");

  useEffect(() => {
    Promise.all([
      getJson<{items: Model[]}>("/api/phenotype/models"),
      getJson<{items: Compound[]}>("/api/phenotype/compounds"),
    ]).then(([m,c]) => {setModels(m.items); setCompounds(c.items);})
      .catch(e => setError(e instanceof Error ? e.message : "Ошибка API"))
      .finally(()=>setLoading(false));
  }, []);

  const filtered = useMemo(()=>compounds.filter(x=>x.compound_id.toLowerCase().includes(search.toLowerCase())),[compounds,search]);
  const current = models.find(x=>x.cell_line_id===variant);
  const selectedModels = [variant, ...(mrc?["MRC-5"]:[]), ...(bj?["BJ5ta"]:[])];
  function toggle(id: string) {setSelected(old=>old.includes(id)?old.filter(x=>x!==id):[...old,id]);}
  function reportError(e: unknown) {setError(e instanceof Error?e.message:"Не удалось выполнить операцию");}

  async function create() {
    setBusy(true);setError("");setNotice("");
    try {
      const created=await getJson<{experiment_id:string}>("/api/phenotype/experiments",{
        method:"POST", headers:{"Content-Type":"application/json"},
        body:JSON.stringify({title,protocol_id:protocol,cell_line_ids:selectedModels,compound_ids:selected,exposure_h:hours,assay_type:"MTT"}),
      });
      const detail=await getJson<Experiment>("/api/phenotype/experiments/"+created.experiment_id);
      setExperimentId(created.experiment_id);setExperiment(detail);setLastImport(null);
      setNotice("Проект создан. Наличие культур и активность PYZ не подтверждены.");
    } catch(e) {reportError(e);} finally {setBusy(false);}
  }
  async function open() {
    setBusy(true);setError("");setNotice("");
    try {setExperiment(await getJson<Experiment>("/api/phenotype/experiments/"+encodeURIComponent(experimentId.trim())));setLastImport(null);}
    catch(e){reportError(e);}finally{setBusy(false);}
  }
  async function upload() {
    if(!file||!experiment)return;
    setBusy(true);setError("");setNotice("");
    try {
      if(!file.name.toLowerCase().endsWith(".csv"))throw new Error("Пока поддерживается только CSV UTF-8");
      if(file.size > 8_000_000) throw new Error("CSV не должен превышать 8 МБ");
      const bytes = new Uint8Array(await file.arrayBuffer());
      let binary = "";
      for(let i=0;i<bytes.length;i+=16384) binary += String.fromCharCode(...bytes.subarray(i,i+16384));
      const imported=await getJson<Imported>("/api/phenotype/experiments/"+experiment.experiment_id+"/import",{
        method:"POST",headers:{"Content-Type":"application/json"},
        body:JSON.stringify({filename:file.name,file_base64:btoa(binary)}),
      });
      setLastImport(imported);
      setExperiment(await getJson<Experiment>("/api/phenotype/experiments/"+experiment.experiment_id));
      setFile(null);
      setNotice(imported.reused_existing_import?"Файл уже импортирован; данные не изменены.":"Исходный файл сохранён. QC и научное заключение пока не утверждены.");
    }catch(e){reportError(e);}finally{setBusy(false);}
  }

  return <div className={styles.page}>
    <header className={styles.hero}>
      <span className={styles.kicker}>MCL SCIENTIFIC CORE v2.1 · P0 / P1</span>
      <h1>Исследование клеточной активности</h1>
      <p>Структура PYZ → собственные измерения → контроль качества. Доступно без предполагаемой мишени или порога Tanimoto.</p>
      <div className={styles.warning}>HEK293T/FT — трансформированные исследовательские линии, а не рак почки.
        MRC-5 — фибробласты лёгкого, а не универсальная модель безопасности.</div>
    </header>
    {error&&<div className={styles.error} role="alert">{error}</div>}
    {notice&&<div className={styles.notice} role="status">{notice}</div>}
    {loading&&<section className={styles.card}>Загружаем реестры…</section>}
    <div className={styles.grid}>
      <section className={styles.card}>
        <span className={styles.kicker}>01 · МОДЕЛИ</span><h2>Исследуемые клетки</h2>
        <label>Точный вариант HEK293</label>
        <select className={styles.input} value={variant} onChange={e=>setVariant(e.target.value)}>
          {models.filter(m=>m.cell_line_id.startsWith("HEK")).map(m=><option key={m.cell_line_id} value={m.cell_line_id}>{m.display_name}</option>)}
        </select>
        {current&&<p className={styles.muted}>{current.note_ru} <strong>{current.inventory_status}</strong></p>}
        <label className={styles.check}><input type="checkbox" checked={mrc} onChange={e=>setMrc(e.target.checked)}/>
          MRC-5 — нормальные фибробласты (наличие не подтверждено)</label>
        <label className={styles.check}><input type="checkbox" checked={bj} onChange={e=>setBj(e.target.checked)}/>
          BJ5ta — дополнительное сравнение</label>
        <p className={styles.muted}>До начала опыта нужны идентификация культур, проверка микоплазмы, пассаж и подтверждение физического наличия.</p>
      </section>
      <section className={styles.card}>
        <span className={styles.kicker}>02 · СОЕДИНЕНИЯ</span><h2>Библиотека PYZ</h2>
        <p className={styles.muted}>Предвыбранные вещества — примеры для планирования, не доказанные лидеры.</p>
        <div className={styles.row}>
          <input className={styles.input} placeholder="Поиск PYZ-…" value={search} onChange={e=>setSearch(e.target.value)}/>
          <button className={styles.small} onClick={()=>setSelected(compounds.map(x=>x.compound_id))}>Все 65</button>
          <button className={styles.small} onClick={()=>setSelected([])}>Снять</button>
        </div>
        <div className={styles.compounds}>
          {filtered.map(c=><label className={styles.check} key={c.compound_id}>
            <input type="checkbox" checked={selected.includes(c.compound_id)} onChange={()=>toggle(c.compound_id)}/>
            {c.compound_id}
          </label>)}
        </div>
        <p className={styles.muted}>Выбрано: {selected.length}. Наличие физических партий и их чистота пока не подтверждены.</p>
      </section>
    </div>
    <section className={styles.card}>
      <span className={styles.kicker}>03 · ПЛАН ОПЫТА</span><h2>Создать эксперимент</h2>
      <div className={styles.formGrid}>
        <label>Название<input className={styles.input} value={title} onChange={e=>setTitle(e.target.value)}/></label>
        <label>Идентификатор протокола<input className={styles.input} value={protocol} onChange={e=>setProtocol(e.target.value)}/></label>
        <label>Время воздействия, ч<input className={styles.input} type="number" min={1} max={720} value={hours} onChange={e=>setHours(Number(e.target.value))}/></label>
      </div>
      <button className={styles.primary} onClick={create} disabled={busy||loading||!models.length||!selected.length}>
        {busy?"Обработка…":"Создать проект"}</button>
    </section>
    <section className={styles.card}>
      <span className={styles.kicker}>04 · ИСХОДНЫЕ ИЗМЕРЕНИЯ</span><h2>Импорт 96-луночного планшета</h2>
      <div className={styles.row}>
        <input className={styles.input} placeholder="PHENO-xxxxxxxxxxxx" value={experimentId} onChange={e=>setExperimentId(e.target.value)}/>
        <button className={styles.small} onClick={open} disabled={busy||!experimentId.trim()}>Открыть</button>
      </div>
      <p className={styles.muted}>Принимается CSV UTF-8. Данные хранятся с SHA-256, повторы файлов и лунок проверяются.</p>
      <p><a className={styles.template} href="/templates/phenotype_96well_template.csv" download>Скачать шаблон CSV</a></p>
      {experiment&&<div className={styles.summary}>
        <strong>{experiment.title}</strong>
        <code>{experiment.experiment_id}</code>
        <span>Культуры: {experiment.cell_line_ids.join(", ")}</span>
        <span>Файлов: {experiment.imported_files_n} · Лунок: {experiment.raw_wells_n} · Статус: {experiment.experiment_status}</span>
        <div className={styles.row}>
          <input type="file" accept=".csv,text/csv" onChange={e=>setFile(e.target.files?.[0]||null)}/>
          <button className={styles.primary} disabled={!file||busy} onClick={upload}>Импортировать CSV</button>
        </div>
        {lastImport&&<div className={styles.qc}>
          <h3>Контроль качества: {lastImport.qc.status}</h3>
          <p>Лунок: {lastImport.wells_n}. <code>SHA-256: {lastImport.sha256}</code></p>
          <ul>{lastImport.qc.flags.map((f,i)=><li key={i}>
            {f.severity}: {f.code}{f.scope?" · "+f.scope:""}{f.well?" · "+f.well:""}
          </li>)}</ul>
          <p>Независимое подтверждение и пороги лабораторного протокола пока отсутствуют.</p>
        </div>}
        <div className={styles.warning}>IC₅₀/GR₅₀: не рассчитаны. Избирательность: {experiment.selectivity_status}.
          Загрузка измерений не доказывает гибель клеток или действие через CDK4.</div>
      </div>}
    </section>
  </div>;
}
