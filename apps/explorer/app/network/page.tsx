import { apiGet } from "../../lib/api";

type Graph={nodes:{id:string;type:string;label:string}[];edges:{id:string;source:string;target:string}[]};
export default async function NetworkPage(){
  const graph=await apiGet<Graph>("/api/network?stable_only=true");
  const genes=graph.nodes.filter(x=>x.type==="gene");
  const pathways=graph.nodes.filter(x=>x.type==="pathway");
  return <><h1>Gene–Pathway Network</h1><p className="muted">Первая серверная версия сети. Интерактивный Cytoscape.js слой будет следующим шагом.</p>
  <section className="grid cards"><div className="card"><div className="label">Gene nodes</div><div className="value">{genes.length}</div></div><div className="card"><div className="label">Pathway nodes</div><div className="value">{pathways.length}</div></div><div className="card"><div className="label">Edges</div><div className="value">{graph.edges.length}</div></div></section>
  <section className="section grid cards">{pathways.map(p=><div className="card" key={p.id}><div className="label">PATHWAY / COMPLEX</div><h3>{p.label}</h3><div className="muted">{graph.edges.filter(e=>e.target===p.id).map(e=>graph.nodes.find(n=>n.id===e.source)?.label).filter(Boolean).join(" · ")}</div></div>)}</section>
  </>;
}
