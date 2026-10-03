import Link from "next/link";
import { apiGet } from "../../lib/api";

type Gene = Record<string, any>;
export default async function GenesPage(){
  const items = await apiGet<Gene[]>("/api/genes/stable");
  return <><h1>Stable Gene Explorer</h1><p className="muted">Гены, которые сохраняют recurrent-статус при Top-50, Top-100 и Top-200.</p>
    <div className="table-wrap"><table><thead><tr><th>Gene</th><th>Top-50</th><th>Top-100</th><th>Top-200</th><th>Thresholds</th></tr></thead><tbody>
    {items.map((g:any)=><tr key={g.gene_symbol}><td><Link href={`/genes/${g.gene_symbol}`} style={{color:"var(--accent)",fontWeight:800}}>{g.gene_symbol}</Link></td><td>{String(g.recurrent_top50).toLowerCase()==="true"?"●":"—"}</td><td>{String(g.recurrent_top100).toLowerCase()==="true"?"●":"—"}</td><td>{String(g.recurrent_top200).toLowerCase()==="true"?"●":"—"}</td><td>{g.thresholds_present || g.thresholds_n}</td></tr>)}
    </tbody></table></div>
  </>;
}
