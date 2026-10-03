"use client";

import { useMemo } from "react";

type Node = { id:string; type:string; label:string };
type Edge = { id:string; source:string; target:string };
type Graph = { nodes:Node[]; edges:Edge[] };

export default function NetworkMap({ graph }:{graph:Graph}){
  const layout = useMemo(()=>{
    const pathways = graph.nodes.filter(n=>n.type==="pathway");
    const genes = graph.nodes.filter(n=>n.type==="gene");
    const width = 1000;
    const height = 520;
    const pos:Record<string,{x:number;y:number}> = {};
    pathways.forEach((p,i)=>{
      const x = pathways.length===1 ? width/2 : 250 + i*(500/Math.max(1,pathways.length-1));
      pos[p.id]={x,y:height/2};
    });
    genes.forEach((g,i)=>{
      const angle = (Math.PI*2*i)/Math.max(genes.length,1) - Math.PI/2;
      pos[g.id]={x:width/2 + Math.cos(angle)*360, y:height/2 + Math.sin(angle)*190};
    });
    return {width,height,pos};
  },[graph]);

  return <div className="network-stage">
    <svg className="network-svg" viewBox={`0 0 ${layout.width} ${layout.height}`} preserveAspectRatio="none">
      {graph.edges.map(e=>{
        const a=layout.pos[e.source], b=layout.pos[e.target];
        if(!a||!b) return null;
        return <line key={e.id} x1={a.x} y1={a.y} x2={b.x} y2={b.y}/>;
      })}
    </svg>
    {graph.nodes.map(n=>{
      const p=layout.pos[n.id]; if(!p) return null;
      return <div key={n.id} className={`network-node ${n.type}`} style={{left:`${(p.x/layout.width)*100}%`,top:`${(p.y/layout.height)*100}%`}}>{n.label}</div>;
    })}
  </div>;
}
