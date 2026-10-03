import { apiGet } from "../../lib/api";
import NetworkMap from "../../components/NetworkMap";

type Graph={nodes:{id:string;type:string;label:string}[];edges:{id:string;source:string;target:string}[]};

export default async function NetworkPage(){
  const graph=await apiGet<Graph>("/api/network?stable_only=true");
  const genes=graph.nodes.filter(x=>x.type==="gene");
  const pathways=graph.nodes.filter(x=>x.type==="pathway");
  return <>
    <div className="section-header" style={{marginBottom:18}}>
      <div>
        <div className="eyebrow">RELATIONSHIP VIEW</div>
        <h1 style={{margin:"4px 0 8px"}}>Gene–Pathway Network</h1>
        <div className="section-copy">Сеть связывает устойчивые функциональные термины с генами, которые формируют enrichment-сигнал. Это визуальный слой над текущими M3.3.1 результатами.</div>
      </div>
    </div>

    <section className="kpi-strip" style={{marginTop:0}}>
      <div className="kpi"><strong>{genes.length}</strong><span>gene nodes</span></div>
      <div className="kpi"><strong>{pathways.length}</strong><span>pathway / complex nodes</span></div>
      <div className="kpi"><strong>{graph.edges.length}</strong><span>gene–term edges</span></div>
      <div className="kpi"><strong>3/3</strong><span>threshold stability filter</span></div>
    </section>

    <section className="section"><NetworkMap graph={graph}/></section>

    <section className="section pathway-grid">
      {pathways.map(p=><article className="pathway-card" key={p.id}>
        <div className="eyebrow">PATHWAY / COMPLEX</div>
        <h3>{p.label}</h3>
        <div className="gene-cloud">{graph.edges.filter(e=>e.target===p.id).map(e=>graph.nodes.find(n=>n.id===e.source)?.label).filter(Boolean).map(g=><span className="gene-pill" key={String(g)}>{g}</span>)}</div>
      </article>)}
    </section>
  </>;
}
