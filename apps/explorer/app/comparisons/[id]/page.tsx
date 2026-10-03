import Link from "next/link";
import DependencyScatter from "../../../components/DependencyScatter";
import ResearchTrail from "../../../components/ResearchTrail";
import { apiGet, formatNumber } from "../../../lib/api";

type Gene = Record<string, any>;
type GenePage = { comparison_id:string; total:number; items:Gene[] };
type Comparison = Record<string, any>;

function qcLabel(status:string){
  const s=String(status||"PASS").toUpperCase();
  if(s==="ERROR") return "критическая ошибка";
  if(s==="WARNING") return "есть ограничения";
  return "проверки пройдены";
}

export default async function ComparisonPage({ params }: { params: Promise<{id:string}> }) {
  const {id} = await params;
  const comparisonId = decodeURIComponent(id);
  const [meta, genes, scatter] = await Promise.all([
    apiGet<Comparison>(`/api/comparisons/${encodeURIComponent(comparisonId)}`),
    apiGet<GenePage>(`/api/comparisons/${encodeURIComponent(comparisonId)}/genes?page_size=100&exclude_broad=true&exclude_low_sample=true`),
    apiGet<GenePage>(`/api/comparisons/${encodeURIComponent(comparisonId)}/genes?page_size=2000&exclude_broad=true&exclude_low_sample=true`),
  ]);

  const stableN = scatter.items.filter((g:any)=>String(g.present_all_thresholds).toLowerCase()==="true").length;

  return <>
    <section className="hero">
      <div className="eyebrow hero-eyebrow">ШАГ 2 · {meta.cancer_id} · DepMap {meta.depmap_release}</div>
      <h1>{meta.label}</h1>
      <p><b>Научный вопрос:</b> какие генетические зависимости сильнее выражены в целевой группе клеточных моделей по сравнению с группой сравнения?</p>
    </section>

    <ResearchTrail current={2} />

    <section className="comparison-question">
      <div><span>Целевая группа</span><b>{meta.context_definition || "Молекулярный контекст"}</b><small>{meta.context_models_n ?? "—"} клеточных моделей</small></div>
      <div className="question-arrow">сравниваем с</div>
      <div><span>Группа сравнения</span><b>{meta.comparator_definition || "Контрольная группа"}</b><small>{meta.comparator_models_n ?? "—"} клеточных моделей</small></div>
    </section>

    <section className="kpi-strip">
      <div className="kpi"><strong>{formatNumber(meta.genes_analyzed_n,0)}</strong><span>генов проанализировано</span></div>
      <div className="kpi"><strong>{formatNumber(scatter.total,0)}</strong><span>генов после текущих фильтров пригодности</span></div>
      <div className="kpi"><strong>{stableN}</strong><span>устойчивых генов среди показанных точек</span></div>
      <div className="kpi"><strong>{qcLabel(meta.qc_status)}</strong><span>итог контроля качества</span></div>
    </section>

    <section className="section split">
      <div>
        <div className="section-header">
          <div>
            <div className="eyebrow">КАРТА ЗАВИСИМОСТЕЙ</div>
            <h2>Чем отличаются две группы клеток?</h2>
            <div className="section-copy">Каждая точка — ген. По горизонтали показан медианный <abbr title="Gene Effect — изменение жизнеспособности клеток после CRISPR-выключения гена">Gene Effect</abbr> в группе сравнения, по вертикали — в целевой группе. Точки ниже диагонали сильнее необходимы целевой группе.</div>
          </div>
        </div>
        <DependencyScatter points={scatter.items} />
      </div>
      <aside className="callout">
        <div className="eyebrow">КАК ИНТЕРПРЕТИРОВАТЬ</div>
        <h3>Сильная CRISPR-зависимость — это ещё не готовая лекарственная мишень</h3>
        <p className="section-copy">Здесь мы ищем биологическое отличие между группами. Гены с широкой клеточной зависимостью и записи с ограниченной выборкой исключены из этого визуального слоя, но остаются в исходных данных и контроле качества.</p>
        <div className="definition-list">
          <div><b>Ниже диагонали</b><span>выключение гена сильнее вредит целевой группе</span></div>
          <div><b>Тёмные точки</b><span>устойчивые повторяющиеся гены при Top-50/100/200</span></div>
          <div><b>Gene Effect</b><span>эффект CRISPR-выключения гена на жизнеспособность клетки</span></div>
        </div>
      </aside>
    </section>

    <section className="section">
      <div className="section-header">
        <div>
          <div className="eyebrow">КАНДИДАТЫ ДЛЯ ДАЛЬНЕЙШЕГО РАЗБОРА</div>
          <h2>Наиболее выраженные зависимости</h2>
          <div className="section-copy">Первые 100 генов после исключения широкой клеточной зависимости и ограниченной выборки. Чем отрицательнее Δ Gene Effect, тем сильнее зависимость сдвинута в сторону целевой группы.</div>
        </div>
      </div>
      <div className="table-wrap"><table><thead><tr><th>Ген</th><th><abbr title="Разница медианного Gene Effect между целевой группой и группой сравнения">Δ Gene Effect</abbr></th><th>Целевая группа</th><th>Группа сравнения</th><th><abbr title="Размер эффекта: насколько хорошо разделяются две группы по зависимости от гена">Cliff's δ</abbr></th><th><abbr title="Статистическая значимость после поправки на множественные проверки">q-value / FDR</abbr></th><th>Устойчивость</th></tr></thead><tbody>
      {genes.items.map((g:any) => <tr key={g.gene_symbol}><td><Link href={`/genes/${g.gene_symbol}`} style={{color:"var(--accent)",fontWeight:900}}>{g.gene_symbol}</Link></td><td>{formatNumber(g.delta_gene_effect)}</td><td>{formatNumber(g.context_median_gene_effect)}</td><td>{formatNumber(g.comparator_median_gene_effect)}</td><td>{formatNumber(g.cliffs_delta)}</td><td>{formatNumber(g.q_value,4)}</td><td>{String(g.present_all_thresholds).toLowerCase()==="true"?<span className="badge stable">устойчив при 50/100/200</span>:"—"}</td></tr>)}
      </tbody></table></div>
    </section>

    <section className="section next-step-banner">
      <div><div className="eyebrow">СЛЕДУЮЩИЙ ШАГ</div><h3>Не выбирайте мишень только по этой таблице</h3><p className="section-copy">Перейдите к устойчивым генам: там мы проверяем, сохраняется ли кандидат при разных порогах отбора и повторяется ли в нескольких сравнениях.</p></div>
      <Link href="/genes" className="primary-link">Перейти к устойчивым генам →</Link>
    </section>
  </>;
}
