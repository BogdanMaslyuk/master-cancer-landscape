import Link from "next/link";
import { apiGet } from "../../lib/api";
import NetworkMap from "../../components/NetworkMap";
import ResearchTrail from "../../components/ResearchTrail";

type Graph={nodes:{id:string;type:string;label:string}[];edges:{id:string;source:string;target:string}[]};

export default async function NetworkPage(){
  const graph=await apiGet<Graph>("/api/network?stable_only=true");
  const genes=graph.nodes.filter(x=>x.type==="gene");
  const pathways=graph.nodes.filter(x=>x.type==="pathway");
  return <>
    <div className="section-header" style={{marginBottom:18}}>
      <div>
        <div className="eyebrow">СВЯЗИ МЕЖДУ КАНДИДАТАМИ</div>
        <h1 style={{margin:"4px 0 8px"}}>Карта генов и функциональных модулей</h1>
        <div className="section-copy">Сеть показывает, какие гены формируют устойчивые функциональные сигналы. Она помогает увидеть не отдельные названия, а биологические группы, которые повторяются в анализе.</div>
      </div>
    </div>

    <ResearchTrail current={4} />

    <section className="kpi-strip" style={{marginTop:18}}>
      <div className="kpi"><strong>{genes.length}</strong><span>генов в текущей сети</span></div>
      <div className="kpi"><strong>{pathways.length}</strong><span>устойчивых функциональных модулей</span></div>
      <div className="kpi"><strong>{graph.edges.length}</strong><span>связей «ген → модуль»</span></div>
      <div className="kpi"><strong>3/3</strong><span>порогов отбора пройдено</span></div>
    </section>

    <section className="section split">
      <div><NetworkMap graph={graph}/></div>
      <aside className="callout">
        <div className="eyebrow">КАК ЧИТАТЬ КАРТУ</div>
        <h3>Ген важнее интерпретировать в контексте системы</h3>
        <div className="definition-list">
          <div><b>Круглый узел</b><span>отдельный ген-кандидат</span></div>
          <div><b>Крупный узел</b><span>функциональный путь или белковый комплекс</span></div>
          <div><b>Связь</b><span>ген входит в пересечение, которое формирует соответствующий сигнал обогащения</span></div>
        </div>
        <p className="section-copy" style={{marginTop:12}}>Карта не доказывает прямое физическое взаимодействие всех соединённых сущностей. Она отображает структуру текущего enrichment-анализа.</p>
      </aside>
    </section>

    <section className="section pathway-grid">
      {pathways.map(p=><article className="pathway-card" key={p.id}>
        <div className="eyebrow">ФУНКЦИОНАЛЬНЫЙ МОДУЛЬ</div>
        <h3>{p.label}</h3>
        <div className="gene-cloud">{graph.edges.filter(e=>e.target===p.id).map(e=>graph.nodes.find(n=>n.id===e.source)?.label).filter(Boolean).map(g=><Link href={`/genes/${g}`} className="gene-pill" key={String(g)}>{g}</Link>)}</div>
      </article>)}
    </section>
  </>;
}
