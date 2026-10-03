import { apiGet } from "../../lib/api";
type QC={counts:Record<string,number>;records:Record<string,any>[]};

export default async function QCPage(){
  const qc=await apiGet<QC>("/api/qc");
  const warnings=qc.records.filter((r:any)=>String(r.severity).toUpperCase()==="WARNING");
  const errors=qc.records.filter((r:any)=>String(r.severity).toUpperCase()==="ERROR");
  const infos=qc.records.filter((r:any)=>String(r.severity).toUpperCase()==="INFO");
  return <>
    <div className="section-header" style={{marginBottom:18}}>
      <div>
        <div className="eyebrow">QUALITY CONTROL</div>
        <h1 style={{margin:"4px 0 8px"}}>QC Dashboard</h1>
        <div className="section-copy">ERROR, WARNING и INFO разделены. WARNING — ограничение интерпретации или качества данных, а не автоматически сломанный анализ.</div>
      </div>
    </div>

    <section className="qc-summary">
      <div className="card qc-card error"><div className="label">ERROR</div><div className="value">{qc.counts.ERROR||0}</div><div className="section-copy">Блокирующие проблемы анализа.</div></div>
      <div className="card qc-card warning"><div className="label">WARNING</div><div className="value">{qc.counts.WARNING||0}</div><div className="section-copy">Ограничения, которые нужно учитывать при интерпретации.</div></div>
      <div className="card qc-card info"><div className="label">INFO</div><div className="value">{qc.counts.INFO||0}</div><div className="section-copy">Контрольные и описательные записи.</div></div>
    </section>

    {errors.length>0 && <section className="section"><h2>Errors</h2><div className="table-wrap" style={{marginTop:12}}><table><thead><tr><th>Check</th><th>Entity</th><th>Observed</th><th>Expected</th><th>Message</th></tr></thead><tbody>{errors.map((r:any,i:number)=><tr key={i}><td>{r.check}</td><td>{r.entity_id}</td><td>{r.observed}</td><td>{r.expected}</td><td>{r.message}</td></tr>)}</tbody></table></div></section>}

    <section className="section">
      <div className="section-header"><div><h2>Warnings requiring interpretation</h2><div className="section-copy">Сначала показываются только предупреждения — это наиболее полезный слой для научного контроля.</div></div></div>
      <div className="table-wrap"><table><thead><tr><th>Check</th><th>Entity</th><th>Observed</th><th>Expected</th><th>Message</th><th>Source</th></tr></thead><tbody>{warnings.map((r:any,i:number)=><tr key={i}><td><span className="badge warning">WARNING</span><div style={{marginTop:6,fontWeight:800}}>{r.check}</div></td><td>{r.entity_id}</td><td>{r.observed}</td><td>{r.expected}</td><td>{r.message}</td><td className="muted">{r.source_file}</td></tr>)}</tbody></table></div>
    </section>

    <section className="section">
      <details className="card">
        <summary style={{cursor:"pointer",fontWeight:850}}>Показать {infos.length} INFO-записей</summary>
        <div className="table-wrap" style={{marginTop:14}}><table><thead><tr><th>Check</th><th>Entity</th><th>Observed</th><th>Expected</th><th>Message</th></tr></thead><tbody>{infos.map((r:any,i:number)=><tr key={i}><td>{r.check}</td><td>{r.entity_id}</td><td>{r.observed}</td><td>{r.expected}</td><td>{r.message}</td></tr>)}</tbody></table></div>
      </details>
    </section>
  </>;
}
