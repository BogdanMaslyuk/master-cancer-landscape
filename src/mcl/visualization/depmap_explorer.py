from __future__ import annotations

import html
import json
from pathlib import Path

import pandas as pd


def _records(df: pd.DataFrame) -> list[dict]:
    if df.empty:
        return []
    return json.loads(df.to_json(orient="records"))


def _json_for_script(value: object) -> str:
    # Prevent an accidental </script> sequence from terminating the data block.
    return json.dumps(value, ensure_ascii=False, separators=(",", ":")).replace("</", "<\\/")


def build_depmap_explorer_html(
    results: pd.DataFrame,
    cohort: pd.DataFrame,
    gene_effect: pd.DataFrame,
    meta: dict,
    output_path: Path,
    *,
    top_n_heatmap: int = 30,
) -> Path:
    """Build a self-contained, offline HTML explorer for one DepMap context comparison."""
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    if results.empty:
        raise ValueError("Cannot build explorer from an empty genome-wide result table")

    has_comparator = int(meta.get("comparator_models_n", 0)) > 0
    if has_comparator and "delta_gene_effect" in results:
        top = results.dropna(subset=["delta_gene_effect"]).nsmallest(top_n_heatmap, "delta_gene_effect")
    else:
        top = results.dropna(subset=["context_median_gene_effect"]).nsmallest(
            top_n_heatmap, "context_median_gene_effect"
        )
    top_genes = top["gene_symbol"].astype(str).tolist()

    if "analysis_group" in cohort.columns:
        selected_models = cohort[cohort["analysis_group"].isin(["context", "comparator"])].copy()
    else:
        selected_models = cohort.iloc[0:0].copy()
    if selected_models.empty:
        model_ids = list(gene_effect.index.astype(str))
        selected_models = pd.DataFrame({"model_id": model_ids, "cell_line_name": model_ids, "analysis_group": "context"})
    model_ids = [x for x in selected_models["model_id"].astype(str).tolist() if x in gene_effect.index]

    heatmap_rows: list[dict] = []
    for gene in top_genes:
        values = []
        for model_id in model_ids:
            value = gene_effect.at[model_id, gene] if gene in gene_effect.columns else None
            values.append(None if pd.isna(value) else float(value))
        heatmap_rows.append({"gene": gene, "values": values})

    model_payload = []
    cohort_by_id = selected_models.set_index("model_id", drop=False)
    for model_id in model_ids:
        row = cohort_by_id.loc[model_id]
        if isinstance(row, pd.DataFrame):
            row = row.iloc[0]
        model_payload.append(
            {
                "model_id": model_id,
                "name": str(row.get("cell_line_name") or model_id),
                "group": str(row.get("analysis_group") or ""),
                "status": str(row.get("molecular_status") or ""),
            }
        )

    explorer_columns = [
        "gene_symbol", "entrez_gene_id", "context_median_gene_effect",
        "comparator_median_gene_effect", "delta_gene_effect", "cliffs_delta",
        "context_dependency_fraction", "comparator_dependency_fraction",
        "broad_dependency_fraction", "broad_dependency_warning",
        "p_value", "q_value", "fdr_0_05", "low_sample_size",
    ]
    result_records = _records(results[[c for c in explorer_columns if c in results.columns]])
    summary = {
        "release": meta.get("depmap_release", ""),
        "cancer_id": meta.get("cancer_id", ""),
        "comparison": meta.get("comparison_type", ""),
        "context_definition": meta.get("context_definition", ""),
        "comparator_definition": meta.get("comparator_definition", ""),
        "context_models_n": meta.get("context_models_n", 0),
        "comparator_models_n": meta.get("comparator_models_n", 0),
        "genes_analyzed_n": meta.get("genes_analyzed_n", len(results)),
        "fdr_genes_n": int(results.get("fdr_0_05", pd.Series(dtype=bool)).fillna(False).sum()),
        "broad_genes_n": int(results.get("broad_dependency_warning", pd.Series(dtype=bool)).fillna(False).sum()),
    }

    data_js = _json_for_script(
        {
            "summary": summary,
            "results": result_records,
            "models": model_payload,
            "heatmap": heatmap_rows,
        }
    )

    title = f"MCL DepMap Explorer — {summary['cancer_id']} / {summary['comparison']}"
    escaped_title = html.escape(title)
    page = f"""<!doctype html>
<html lang="ru">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>{escaped_title}</title>
<style>
:root{{--bg:#f5f7fb;--panel:#fff;--ink:#0f172a;--muted:#64748b;--line:#e2e8f0;--accent:#2563eb;--danger:#b91c1c;--warn:#d97706;--context:#2563eb;--comp:#7c3aed}}
*{{box-sizing:border-box}} body{{margin:0;font-family:Inter,ui-sans-serif,system-ui,-apple-system,"Segoe UI",sans-serif;background:var(--bg);color:var(--ink)}}
header{{padding:28px 34px;background:#0f172a;color:white}} header h1{{margin:0 0 7px;font-size:25px}} header p{{margin:0;color:#cbd5e1;max-width:1100px;line-height:1.45}}
main{{max-width:1500px;margin:0 auto;padding:24px}} .grid{{display:grid;gap:16px}} .cards{{grid-template-columns:repeat(5,minmax(0,1fr));margin-bottom:16px}}
.card,.panel{{background:var(--panel);border:1px solid var(--line);border-radius:16px;box-shadow:0 1px 2px rgba(15,23,42,.04)}} .card{{padding:16px}} .card .k{{font-size:12px;color:var(--muted);text-transform:uppercase;letter-spacing:.05em}} .card .v{{font-size:25px;font-weight:750;margin-top:5px}}
.panel{{padding:20px;margin-bottom:16px}} .panel h2{{font-size:18px;margin:0 0 8px}} .panel .sub{{color:var(--muted);font-size:13px;margin-bottom:16px}}
.two{{grid-template-columns:minmax(0,1.25fr) minmax(360px,.75fr)}} canvas{{width:100%;height:460px;border:1px solid var(--line);border-radius:12px;background:#fff}}
.legend{{display:flex;gap:18px;align-items:center;font-size:12px;color:var(--muted);margin:8px 0 0}} .dot{{width:9px;height:9px;border-radius:50%;display:inline-block;margin-right:5px}}
.profile{{min-height:460px}} .profile .gene{{font-size:32px;font-weight:800;margin:12px 0 18px}} .metric{{display:grid;grid-template-columns:1fr auto;gap:10px;border-top:1px solid var(--line);padding:10px 0;font-size:13px}} .metric b{{font-variant-numeric:tabular-nums}}
.controls{{display:flex;flex-wrap:wrap;gap:10px;margin:12px 0 14px}} input,select{{border:1px solid var(--line);border-radius:9px;padding:9px 10px;background:white;color:var(--ink)}} label.chk{{display:flex;align-items:center;gap:6px;font-size:13px;color:var(--muted)}}
table{{width:100%;border-collapse:collapse;font-size:12px}} th{{position:sticky;top:0;background:#f8fafc;text-align:left;color:#475569;padding:9px;border-bottom:1px solid var(--line)}} td{{padding:8px 9px;border-bottom:1px solid #edf2f7;font-variant-numeric:tabular-nums}} tbody tr{{cursor:pointer}} tbody tr:hover{{background:#f8fafc}} .tag{{display:inline-block;padding:2px 7px;border-radius:999px;font-size:11px;background:#eef2ff;color:#3730a3}} .tag.warn{{background:#fff7ed;color:#9a3412}}
.tablewrap{{max-height:520px;overflow:auto;border:1px solid var(--line);border-radius:12px}}
.heatmap-wrap{{overflow:auto;border:1px solid var(--line);border-radius:12px;padding:10px;background:#fff}} .heat-row{{display:grid;align-items:center;gap:2px;margin:2px 0}} .heat-gene{{font-size:11px;font-weight:650;overflow:hidden;text-overflow:ellipsis;white-space:nowrap;padding-right:6px}} .heat-cell{{height:17px;border-radius:2px;min-width:17px}} .model-head{{font-size:9px;writing-mode:vertical-rl;transform:rotate(180deg);height:76px;overflow:hidden;color:#64748b;white-space:nowrap}}
.note{{font-size:12px;line-height:1.5;color:#64748b;background:#f8fafc;border:1px solid var(--line);border-radius:12px;padding:12px;margin-top:14px}}
#tooltip{{position:fixed;display:none;pointer-events:none;background:#0f172a;color:white;padding:8px 10px;border-radius:8px;font-size:12px;z-index:30;box-shadow:0 8px 24px rgba(0,0,0,.2)}}
@media(max-width:980px){{.cards{{grid-template-columns:repeat(2,1fr)}}.two{{grid-template-columns:1fr}}main{{padding:12px}}header{{padding:22px 18px}}}}
</style>
</head>
<body>
<header><h1>Genome-wide Dependency Explorer</h1><p><b>{html.escape(str(summary['cancer_id']))}</b> · {html.escape(str(summary['context_definition']))}<br>Контроль: {html.escape(str(summary['comparator_definition'] or 'нет — описательный режим'))}</p></header>
<main>
<section class="grid cards">
<div class="card"><div class="k">DepMap</div><div class="v">{html.escape(str(summary['release']))}</div></div>
<div class="card"><div class="k">Гены</div><div class="v">{summary['genes_analyzed_n']:,}</div></div>
<div class="card"><div class="k">Контекст</div><div class="v">{summary['context_models_n']}</div></div>
<div class="card"><div class="k">Контроль</div><div class="v">{summary['comparator_models_n']}</div></div>
<div class="card"><div class="k">FDR q&lt;0,05</div><div class="v">{summary['fdr_genes_n']}</div></div>
</section>
<section class="grid two">
<div class="panel"><h2>Карта геномных зависимостей</h2><div class="sub">Каждая точка — ген. Ось Y: медианный Gene Effect в изучаемом контексте; X: в контрольной группе. Ниже диагонали — более сильная зависимость в контексте.</div><canvas id="scatter" width="1000" height="520"></canvas><div class="legend"><span><i class="dot" style="background:#2563eb"></i>q&lt;0,05 и Δ&lt;0</span><span><i class="dot" style="background:#d97706"></i>широкая зависимость</span><span><i class="dot" style="background:#94a3b8"></i>остальные</span></div></div>
<div class="panel profile"><h2>Карточка гена</h2><div class="sub">Нажмите на точку или строку таблицы.</div><div id="profile"></div></div>
</section>
<section class="panel"><h2>Тепловая карта наиболее контекстно-селективных генов</h2><div class="sub">Показаны гены с наиболее отрицательным Δ медианного Gene Effect; это прозрачная сортировка, а не общий «рейтинг мишени».</div><div class="heatmap-wrap" id="heatmap"></div><div class="note">Красный оттенок означает более отрицательный Gene Effect (более сильную генетическую зависимость). Тепловая карта показывает отдельные клеточные модели и позволяет увидеть, является ли сигнал однородным или создаётся несколькими линиями.</div></section>
<section class="panel"><h2>Все гены</h2><div class="controls"><input id="search" placeholder="Поиск гена…"><select id="sort"><option value="delta">Δ Gene Effect</option><option value="context">Gene Effect контекста</option><option value="q">q-value</option><option value="broad">Широкая зависимость</option></select><label class="chk"><input type="checkbox" id="fdrOnly"> только q&lt;0,05</label><label class="chk"><input type="checkbox" id="hideBroad"> скрыть широкие зависимости</label></div><div class="tablewrap"><table><thead><tr><th>Ген</th><th>GE контекст</th><th>GE контроль</th><th>Δ</th><th>Cliff</th><th>Dep. контекст</th><th>Dep. контроль</th><th>Broad</th><th>p</th><th>q</th></tr></thead><tbody id="tbody"></tbody></table></div><div class="note">CRISPR-нокаут показывает генетическую зависимость, а не автоматически пригодность белка для лекарственного воздействия. Broad dependency — внутренний скрининговый флаг, не заключение о безопасности.</div></section>
</main><div id="tooltip"></div>
<script>const DATA={data_js};
const fmt=(v,d=3)=>v===null||v===undefined||Number.isNaN(v)?'—':Number(v).toFixed(d); const pct=v=>v===null||v===undefined?'—':(100*Number(v)).toFixed(1)+'%'; const sci=v=>v===null||v===undefined?'—':Number(v).toExponential(2);
const byGene=new Map(DATA.results.map(r=>[r.gene_symbol,r])); const tooltip=document.getElementById('tooltip');
function showProfile(r){{if(!r)return;document.getElementById('profile').innerHTML=`<div class="gene">${{r.gene_symbol}}</div><div class="metric"><span>Median Gene Effect — контекст</span><b>${{fmt(r.context_median_gene_effect)}}</b></div><div class="metric"><span>Median Gene Effect — контроль</span><b>${{fmt(r.comparator_median_gene_effect)}}</b></div><div class="metric"><span>Δ Gene Effect</span><b>${{fmt(r.delta_gene_effect)}}</b></div><div class="metric"><span>Cliff's delta</span><b>${{fmt(r.cliffs_delta)}}</b></div><div class="metric"><span>Зависимые модели — контекст</span><b>${{pct(r.context_dependency_fraction)}}</b></div><div class="metric"><span>Зависимые модели — контроль</span><b>${{pct(r.comparator_dependency_fraction)}}</b></div><div class="metric"><span>Широкая зависимость по DepMap</span><b>${{pct(r.broad_dependency_fraction)}}</b></div><div class="metric"><span>p-value</span><b>${{sci(r.p_value)}}</b></div><div class="metric"><span>q-value</span><b>${{sci(r.q_value)}}</b></div>${{r.broad_dependency_warning?'<div class="note">⚠ Ген широко необходим многим моделям DepMap; это уменьшает специфичность сигнала, но само по себе не является оценкой безопасности.</div>':''}}`;}}
function drawScatter(){{const c=document.getElementById('scatter'),ctx=c.getContext('2d'),pts=DATA.results.filter(r=>r.context_median_gene_effect!==null&&r.comparator_median_gene_effect!==null);ctx.clearRect(0,0,c.width,c.height);if(!pts.length){{ctx.fillStyle='#64748b';ctx.font='16px sans-serif';ctx.fillText('Для этого контекста нет контрольной группы: scatter недоступен.',30,50);return;}} const vals=pts.flatMap(r=>[r.context_median_gene_effect,r.comparator_median_gene_effect]);let lo=Math.min(...vals),hi=Math.max(...vals);lo=Math.min(lo,-.2);hi=Math.max(hi,.2);const pad=55,sx=v=>pad+(v-lo)/(hi-lo)*(c.width-2*pad),sy=v=>c.height-pad-(v-lo)/(hi-lo)*(c.height-2*pad);ctx.strokeStyle='#e2e8f0';ctx.lineWidth=1;for(let i=0;i<=5;i++){{const v=lo+i*(hi-lo)/5;ctx.beginPath();ctx.moveTo(sx(v),pad);ctx.lineTo(sx(v),c.height-pad);ctx.stroke();ctx.beginPath();ctx.moveTo(pad,sy(v));ctx.lineTo(c.width-pad,sy(v));ctx.stroke();ctx.fillStyle='#64748b';ctx.font='11px sans-serif';ctx.fillText(v.toFixed(1),sx(v)-10,c.height-pad+18);ctx.fillText(v.toFixed(1),8,sy(v)+4);}}ctx.strokeStyle='#94a3b8';ctx.setLineDash([5,5]);ctx.beginPath();ctx.moveTo(sx(lo),sy(lo));ctx.lineTo(sx(hi),sy(hi));ctx.stroke();ctx.setLineDash([]);const hit=[];for(const r of pts){{const x=sx(r.comparator_median_gene_effect),y=sy(r.context_median_gene_effect);let color='#94a3b8';if(r.broad_dependency_warning)color='#d97706';if(r.q_value!==null&&r.q_value<.05&&r.delta_gene_effect<0)color='#2563eb';ctx.fillStyle=color;ctx.globalAlpha=.75;ctx.beginPath();ctx.arc(x,y,(r.q_value!==null&&r.q_value<.05)?4:2.5,0,Math.PI*2);ctx.fill();hit.push({{x,y,r}});}}ctx.globalAlpha=1;ctx.fillStyle='#334155';ctx.font='12px sans-serif';ctx.fillText('Median Gene Effect — контроль',c.width/2-75,c.height-12);ctx.save();ctx.translate(14,c.height/2+60);ctx.rotate(-Math.PI/2);ctx.fillText('Median Gene Effect — контекст',0,0);ctx.restore();c.onmousemove=e=>{{const rect=c.getBoundingClientRect(),mx=(e.clientX-rect.left)*c.width/rect.width,my=(e.clientY-rect.top)*c.height/rect.height;let best=null,dist=12;for(const h of hit){{const d=Math.hypot(h.x-mx,h.y-my);if(d<dist){{best=h;dist=d}}}}if(best){{tooltip.style.display='block';tooltip.style.left=(e.clientX+12)+'px';tooltip.style.top=(e.clientY+12)+'px';tooltip.innerHTML=`<b>${{best.r.gene_symbol}}</b><br>контекст: ${{fmt(best.r.context_median_gene_effect)}}<br>контроль: ${{fmt(best.r.comparator_median_gene_effect)}}<br>Δ: ${{fmt(best.r.delta_gene_effect)}}<br>q: ${{sci(best.r.q_value)}}`;}}else tooltip.style.display='none';}};c.onclick=e=>{{const rect=c.getBoundingClientRect(),mx=(e.clientX-rect.left)*c.width/rect.width,my=(e.clientY-rect.top)*c.height/rect.height;let best=null,dist=14;for(const h of hit){{const d=Math.hypot(h.x-mx,h.y-my);if(d<dist){{best=h;dist=d}}}}if(best)showProfile(best.r);}};c.onmouseleave=()=>tooltip.style.display='none';}}
function heatColor(v){{if(v===null||v===undefined)return '#e2e8f0';v=Math.max(-3,Math.min(1,Number(v)));if(v<=0){{const t=Math.min(1,Math.abs(v)/2);const r=Math.round(255-(255-153)*t),g=Math.round(255-(255-27)*t),b=Math.round(255-(255-27)*t);return `rgb(${{r}},${{g}},${{b}})`;}}const t=Math.min(1,v);return `rgb(${{Math.round(255-64*t)}},${{Math.round(255-119*t)}},255)`;}}
function drawHeatmap(){{const el=document.getElementById('heatmap'),cols=DATA.models.length;let h=`<div class="heat-row" style="grid-template-columns:110px repeat(${{cols}},minmax(17px,1fr))"><div></div>`;for(const m of DATA.models)h+=`<div class="model-head" title="${{m.name}} · ${{m.group}}">${{m.name}}</div>`;h+='</div>';for(const row of DATA.heatmap){{h+=`<div class="heat-row" style="grid-template-columns:110px repeat(${{cols}},minmax(17px,1fr))"><div class="heat-gene">${{row.gene}}</div>`;row.values.forEach((v,i)=>{{const m=DATA.models[i];h+=`<div class="heat-cell" style="background:${{heatColor(v)}}" title="${{row.gene}} · ${{m.name}} · ${{m.group}} · Gene Effect ${{fmt(v)}}"></div>`}});h+='</div>';}}el.innerHTML=h;}}
function renderTable(){{let rows=[...DATA.results],q=document.getElementById('search').value.trim().toUpperCase(),fdr=document.getElementById('fdrOnly').checked,hideBroad=document.getElementById('hideBroad').checked,sort=document.getElementById('sort').value;if(q)rows=rows.filter(r=>r.gene_symbol.toUpperCase().includes(q));if(fdr)rows=rows.filter(r=>r.q_value!==null&&r.q_value<.05);if(hideBroad)rows=rows.filter(r=>!r.broad_dependency_warning);const val=(r,k)=>r[k]===null||r[k]===undefined?Infinity:Number(r[k]);if(sort==='delta')rows.sort((a,b)=>val(a,'delta_gene_effect')-val(b,'delta_gene_effect'));if(sort==='context')rows.sort((a,b)=>val(a,'context_median_gene_effect')-val(b,'context_median_gene_effect'));if(sort==='q')rows.sort((a,b)=>val(a,'q_value')-val(b,'q_value'));if(sort==='broad')rows.sort((a,b)=>val(b,'broad_dependency_fraction')-val(a,'broad_dependency_fraction'));rows=rows.slice(0,500);const tb=document.getElementById('tbody');tb.innerHTML=rows.map(r=>`<tr data-gene="${{r.gene_symbol}}"><td><b>${{r.gene_symbol}}</b> ${{r.fdr_0_05?'<span class="tag">FDR</span>':''}} ${{r.broad_dependency_warning?'<span class="tag warn">broad</span>':''}}</td><td>${{fmt(r.context_median_gene_effect)}}</td><td>${{fmt(r.comparator_median_gene_effect)}}</td><td>${{fmt(r.delta_gene_effect)}}</td><td>${{fmt(r.cliffs_delta)}}</td><td>${{pct(r.context_dependency_fraction)}}</td><td>${{pct(r.comparator_dependency_fraction)}}</td><td>${{pct(r.broad_dependency_fraction)}}</td><td>${{sci(r.p_value)}}</td><td>${{sci(r.q_value)}}</td></tr>`).join('');tb.querySelectorAll('tr').forEach(tr=>tr.onclick=()=>showProfile(byGene.get(tr.dataset.gene)));}}
['search','sort','fdrOnly','hideBroad'].forEach(id=>document.getElementById(id).addEventListener(id==='search'?'input':'change',renderTable));drawScatter();drawHeatmap();renderTable();showProfile(DATA.results[0]);
</script></body></html>"""
    output_path.write_text(page, encoding="utf-8")
    return output_path
