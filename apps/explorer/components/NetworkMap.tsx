"use client";

import { useMemo, useState } from "react";

type Node = { id:string; type:string; label:string };
type Edge = { id:string; source:string; target:string };
type Graph = { nodes:Node[]; edges:Edge[] };

export default function NetworkMap({ graph }:{graph:Graph}){
  const [selected, setSelected] = useState<string | null>(null);

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

  const related = useMemo(()=>{
    if(!selected) return new Set(graph.nodes.map(n=>n.id));
    const ids = new Set<string>([selected]);
    graph.edges.forEach(e=>{
      if(e.source===selected) ids.add(e.target);
      if(e.target===selected) ids.add(e.source);
    });
    return ids;
  },[graph,selected]);

  const selectedNode = selected ? graph.nodes.find(n=>n.id===selected) : null;

  return <div className="network-stage" onClick={()=>selected && setSelected(null)}>
    <div style={{position:"absolute",left:16,top:14,zIndex:4,padding:"7px 10px",border:"1px solid #e2e8f0",borderRadius:10,background:"rgba(255,255,255,.9)",color:"#66778c",fontSize:10.5,fontWeight:750,backdropFilter:"blur(8px)"}}>
      {selectedNode ? `Фокус: ${selectedNode.label} · нажмите на фон, чтобы сбросить` : "Нажмите на ген или модуль, чтобы выделить его связи"}
    </div>
    <svg className="network-svg" viewBox={`0 0 ${layout.width} ${layout.height}`} preserveAspectRatio="none">
      {graph.edges.map(e=>{
        const a=layout.pos[e.source], b=layout.pos[e.target];
        if(!a||!b) return null;
        const active=!selected || e.source===selected || e.target===selected;
        return <line key={e.id} x1={a.x} y1={a.y} x2={b.x} y2={b.y} style={{opacity:active?0.9:0.1,stroke:active&&selected?"#5f86ad":"#c0ccda",strokeWidth:active&&selected?1.8:1.2,transition:"opacity .15s ease, stroke .15s ease"}}/>;
      })}
    </svg>
    {graph.nodes.map(n=>{
      const p=layout.pos[n.id]; if(!p) return null;
      const active=related.has(n.id);
      const chosen=selected===n.id;
      return <div
        key={n.id}
        className={`network-node ${n.type}`}
        role="button"
        tabIndex={0}
        aria-pressed={chosen}
        onClick={(event)=>{event.stopPropagation();setSelected(chosen?null:n.id)}}
        onKeyDown={(event)=>{if(event.key==="Enter"||event.key===" "){event.preventDefault();setSelected(chosen?null:n.id)}}}
        style={{left:`${(p.x/layout.width)*100}%`,top:`${(p.y/layout.height)*100}%`,opacity:active?1:0.22,cursor:"pointer",outline:chosen?"3px solid rgba(31,95,168,.18)":"none",outlineOffset:3,transition:"opacity .15s ease, transform .15s ease"}}
      >{n.label}</div>;
    })}
  </div>;
}
