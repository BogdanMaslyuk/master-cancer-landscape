import { apiGet } from "../../lib/api";
type QC={counts:Record<string,number>;records:Record<string,any>[]};

function severityLabel(value:string){
  const s=String(value||"INFO").toUpperCase();
  if(s==="ERROR") return "Критическая ошибка";
  if(s==="WARNING") return "Ограничение";
  return "Информационная проверка";
}

export default async function QCPage(){
  const qc=await apiGet<QC>("/api/qc");
  const warnings=qc.records.filter((r:any)=>String(r.severity).toUpperCase()==="WARNING");
  const errors=qc.records.filter((r:any)=>String(r.severity).toUpperCase()==="ERROR");
  const infos=qc.records.filter((r:any)=>String(r.severity).toUpperCase()==="INFO");
  const overall=errors.length?"Есть критические ошибки":warnings.length?"Анализ выполнен, есть ограничения":"Проверки пройдены";
  return <>
    <div className="section-header" style={{marginBottom:18}}>
      <div>
        <div className="eyebrow">КОНТРОЛЬ КАЧЕСТВА</div>
        <h1 style={{margin:"4px 0 8px"}}>Насколько можно доверять текущему анализу?</h1>
        <div className="section-copy">Здесь отдельно показаны критические ошибки, ограничения интерпретации и служебные проверки. Наличие предупреждения не означает, что весь анализ недействителен.</div>
      </div>
    </div>

    <section className="callout">
      <div className="eyebrow">ИТОГ</div>
      <h3>{overall}</h3>
      <p className="section-copy">Для научной интерпретации сначала смотрите раздел «Ограничения». Полный технический журнал нужен для воспроизводимости и аудита.</p>
    </section>

    <section className="qc-summary section">
      <div className="card qc-card error"><div className="label">Критические ошибки</div><div className="value">{qc.counts.ERROR||0}</div><div className="section-copy">Проблемы, которые могут блокировать интерпретацию результата.</div></div>
      <div className="card qc-card warning"><div className="label">Ограничения</div><div className="value">{qc.counts.WARNING||0}</div><div className="section-copy">Факторы, которые нужно учитывать при формулировке вывода.</div></div>
      <div className="card qc-card info"><div className="label">Информационные проверки</div><div className="value">{qc.counts.INFO||0}</div><div className="section-copy">Подтверждения объёма данных, числа генов и других параметров запуска.</div></div>
    </section>

    {errors.length>0 && <section className="section"><h2>Критические ошибки</h2><div className="table-wrap" style={{marginTop:12}}><table><thead><tr><th>Проверка</th><th>Объект</th><th>Наблюдалось</th><th>Ожидалось</th><th>Сообщение</th></tr></thead><tbody>{errors.map((r:any,i:number)=><tr key={i}><td>{r.check}</td><td>{r.entity_id}</td><td>{r.observed}</td><td>{r.expected}</td><td>{r.message}</td></tr>)}</tbody></table></div></section>}

    <section className="section">
      <div className="section-header"><div><div className="eyebrow">СНАЧАЛА СМОТРЕТЬ СЮДА</div><h2>Ограничения интерпретации</h2><div className="section-copy">Это наиболее важная часть контроля качества для исследователя. Она показывает, где вывод нужно формулировать осторожнее.</div></div></div>
      {warnings.length ? <div className="table-wrap"><table><thead><tr><th>Тип</th><th>Проверка</th><th>Объект</th><th>Наблюдалось</th><th>Ожидалось</th><th>Что сообщает система</th></tr></thead><tbody>{warnings.map((r:any,i:number)=><tr key={i}><td><span className="badge warning">{severityLabel(r.severity)}</span></td><td>{r.check}</td><td>{r.entity_id}</td><td>{r.observed}</td><td>{r.expected}</td><td>{r.message}</td></tr>)}</tbody></table></div> : <div className="empty-state">Ограничения уровня WARNING отсутствуют.</div>}
    </section>

    <section className="section">
      <details className="details-block">
        <summary>Показать {infos.length} служебных информационных записей</summary>
        <div className="table-wrap" style={{marginTop:14}}><table><thead><tr><th>Проверка</th><th>Объект</th><th>Наблюдалось</th><th>Ожидалось</th><th>Сообщение</th><th>Файл-источник</th></tr></thead><tbody>{infos.map((r:any,i:number)=><tr key={i}><td>{r.check}</td><td>{r.entity_id}</td><td>{r.observed}</td><td>{r.expected}</td><td>{r.message}</td><td className="muted">{r.source_file}</td></tr>)}</tbody></table></div>
      </details>
    </section>
  </>;
}
