import Link from "next/link";
import ResearchTrail from "../../../components/ResearchTrail";
import { apiGet, formatNumber } from "../../../lib/api";

type Payload = { identity:Record<string,any>; stability:Record<string,any>|null; comparisons:Record<string,any>[]; pathways:Record<string,any>[] };

function yes(value:unknown){ return String(value).toLowerCase()==="true"; }

export default async function GenePage({params}:{params:Promise<{gene:string}>}){
  const {gene}=await params;
  const data=await apiGet<Payload>(`/api/genes/${encodeURIComponent(gene)}`);
  const symbol=data.identity.gene_symbol || gene;
  const stable=yes(data.stability?.present_all_thresholds);
  const negative=data.comparisons.filter((r:any)=>Number(r.delta_gene_effect)<0).length;
  const significantPathways=data.pathways.filter((p:any)=>yes(p.significant));
  const broad=data.comparisons.some((r:any)=>yes(r.broad_dependency_warning));
  const lowSample=data.comparisons.some((r:any)=>yes(r.low_sample_size));

  return <>
    <section className="hero">
      <div className="eyebrow hero-eyebrow">ГЕН-КАНДИДАТ</div>
      <h1>{symbol}</h1>
      <p>Карточка объединяет CRISPR-зависимость, устойчивость между порогами отбора и связанные функциональные сигналы. Она помогает понять, почему ген заслуживает дальнейшей проверки, но не заменяет фармакологическую оценку.</p>
    </section>

    <ResearchTrail current={3} />

    <section className="grid cards">
      <div className="card"><div className="label">Статус устойчивости</div><div className="value" style={{fontSize:22}}>{stable ? "Устойчивый повторяющийся" : "Контекстный кандидат"}</div><div className="section-copy">{stable?"Сохраняется при Top-50, Top-100 и Top-200":"Не проходит критерий устойчивости 3/3"}</div></div>
      <div className="card"><div className="label">Сравнений с более сильной зависимостью в целевой группе</div><div className="value">{negative}</div><div className="section-copy">из {data.comparisons.length} доступных сравнений</div></div>
      <div className="card"><div className="label">Значимых функциональных связей</div><div className="value">{significantPathways.length}</div><div className="section-copy">наблюдений обогащения при разных Top-N</div></div>
      <div className="card"><div className="label">Технические идентификаторы</div><div className="technical-note"><b>HGNC:</b> {data.identity.hgnc_id || data.identity.HGNC_ID || "—"}<br/><b>Ensembl:</b> {data.identity.ensembl_gene_id || data.identity.ensembl_id || "—"}</div></div>
    </section>

    <section className="section split">
      <div className="callout">
        <div className="eyebrow">ПОЧЕМУ ГЕН ПОПАЛ В ПОЛЕ ЗРЕНИЯ</div>
        <h3>{stable ? "Сигнал устойчив к изменению порога отбора" : "Сигнал требует более осторожной интерпретации"}</h3>
        <div className="definition-list">
          <div><b>Устойчивость</b><span>{stable ? "ген остаётся повторяющимся при Top-50, Top-100 и Top-200" : "ген не сохраняет повторяющийся статус при всех трёх порогах"}</span></div>
          <div><b>Сравнения</b><span>в {negative} из {data.comparisons.length} сравнений Δ Gene Effect отрицательна, то есть зависимость сильнее в целевой группе</span></div>
          <div><b>Функциональный контекст</b><span>{significantPathways.length ? `ген входит в ${significantPathways.length} значимых наблюдений функционального обогащения` : "значимых связей с функциональными терминами пока не найдено"}</span></div>
        </div>
      </div>
      <div className="card">
        <div className="eyebrow">ЧТО ОГРАНИЧИВАЕТ ВЫВОД</div>
        <h3>Что пока не доказано</h3>
        <ul className="plain-list">
          <li>CRISPR-нокаут не равен фармакологическому ингибированию белка.</li>
          <li>Текущие данные не доказывают безопасность воздействия на нормальные ткани.</li>
          <li>Лекарственная достижимость и наличие подходящего сайта связывания ещё не оценены.</li>
          {broad && <li>Для гена отмечалась широкая клеточная зависимость — селективность требует особого внимания.</li>}
          {lowSample && <li>В одном или нескольких сравнениях есть ограничение по размеру выборки.</li>}
        </ul>
      </div>
    </section>

    <section className="section">
      <div className="section-header"><div><h2>Результаты по опухолевым сравнениям</h2><div className="section-copy"><abbr title="Gene Effect — изменение жизнеспособности клеток после CRISPR-выключения гена">Gene Effect</abbr> оставлен на английском, потому что это стандартное название метрики DepMap. Остальные поля снабжены русским смысловым описанием.</div></div></div>
      <div className="table-wrap"><table><thead><tr><th>Сравнение</th><th><abbr title="Разница медианного Gene Effect: целевая группа минус группа сравнения">Δ Gene Effect</abbr></th><th>Целевая группа</th><th>Группа сравнения</th><th><abbr title="Размер эффекта между двумя группами">Cliff's δ</abbr></th><th><abbr title="Значимость после поправки на множественные проверки">q-value / FDR</abbr></th><th>Широкая зависимость</th><th>Малая выборка</th></tr></thead><tbody>{data.comparisons.map((r:any)=><tr key={r.comparison_id}><td>{r.comparison_label}</td><td>{formatNumber(r.delta_gene_effect)}</td><td>{formatNumber(r.context_median_gene_effect)}</td><td>{formatNumber(r.comparator_median_gene_effect)}</td><td>{formatNumber(r.cliffs_delta)}</td><td>{formatNumber(r.q_value,4)}</td><td>{yes(r.broad_dependency_warning)?"да":"нет"}</td><td>{yes(r.low_sample_size)?"да":"нет"}</td></tr>)}</tbody></table></div>
    </section>

    <section className="section">
      <div className="section-header"><div><h2>Функциональные связи</h2><div className="section-copy">Показаны только статистически значимые наблюдения обогащения, в пересечение которых входит {symbol}. Повторение одного термина при разных Top-N отражает его устойчивость.</div></div></div>
      {significantPathways.length ? <div className="table-wrap"><table><thead><tr><th>Порог</th><th>Источник</th><th>Функциональный термин</th><th>Скорректированное p</th><th>Гены пересечения</th></tr></thead><tbody>{significantPathways.slice(0,80).map((p:any,i:number)=><tr key={`${p.term_id}-${p.top_n}-${i}`}><td>Top-{p.top_n}</td><td>{p.source}</td><td>{p.term_name}</td><td>{formatNumber(p.p_value_adjusted,5)}</td><td>{p.intersecting_gene_symbols_json}</td></tr>)}</tbody></table></div> : <div className="empty-state">Для этого гена пока нет значимых pathway/complex-связей в текущем анализе.</div>}
    </section>

    <section className="section next-step-banner">
      <div><div className="eyebrow">ЧТО ПРОВЕРЯТЬ ДАЛЬШЕ В M3.4</div><h3>Можно ли превратить {symbol} в лекарственную мишень?</h3><p className="section-copy">Следующий слой должен добавить лекарственную достижимость, известные лиганды и препараты, 3D-структуры и карманы, нормальную тканевую экспрессию, безопасность и внешние онкологические доказательства.</p></div>
      <Link href="/methodology" className="primary-link">Открыть план фармакологической оценки →</Link>
    </section>
  </>;
}
