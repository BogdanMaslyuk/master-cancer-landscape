import CRISPRModelTable from "../../components/CRISPRModelTable";
import { apiGet } from "../../lib/api";

export const dynamic = "force-dynamic";

type ModelsResponse = {
  total: number;
  status?: string;
  canonical_index_available?: boolean;
  build_command?: string | null;
  note_ru?: string;
  items: Record<string, any>[];
};

export default async function ModelsPage({
  searchParams,
}: {
  searchParams: Promise<{ organ?: string; cancer?: string }>;
}) {
  const params = await searchParams;
  const models = await apiGet<ModelsResponse>("/api/crispr-models?limit=5000");
  const organs = new Set(models.items.map((m) => m.mcl_organ_id).filter(Boolean)).size;
  const cancers = new Set(models.items.map((m) => m.mcl_cancer_id).filter(Boolean)).size;
  const curated = models.items.filter((m) => Number(m.curated_contexts_n || 0) > 0).length;

  return <>
    <section className="page-intro compact-intro">
      <div>
        <div className="eyebrow">ИССЛЕДОВАТЕЛЬ ЭКСПЕРИМЕНТАЛЬНЫХ МОДЕЛЕЙ</div>
        <h1>Клеточные линии с CRISPR-профилем</h1>
        <p>Здесь каждая строка — отдельная клеточная модель DepMap, для которой доступен CRISPR Gene Effect. Найдите модель по названию или последовательно сузьте выбор по органу, опухоли и подтипу.</p>
      </div>
    </section>

    <section className="kpi-strip">
      <div className="kpi"><strong>{models.total}</strong><span>уникальных CRISPR-моделей</span></div>
      <div className="kpi"><strong>{organs}</strong><span>органов / систем</span></div>
      <div className="kpi"><strong>{cancers}</strong><span>групп опухолей</span></div>
      <div className="kpi"><strong>{curated}</strong><span>моделей уже входят в контексты MCL</span></div>
    </section>

    {!models.canonical_index_available && <section className="section callout atlas-callout">
      <div className="eyebrow">ВРЕМЕННЫЙ РЕЖИМ</div>
      <h3>Полный индекс CRISPR-моделей ещё не материализован</h3>
      <p className="section-copy">{models.note_ru} После построения локального индекса эта страница автоматически расширится на все модели из CRISPRGeneEffect.csv.</p>
      {models.build_command && <code>{models.build_command}</code>}
    </section>}

    <section className="section">
      <div className="section-header">
        <div>
          <div className="eyebrow">ПОИСК И ФИЛЬТРАЦИЯ</div>
          <h2>Найдите подходящую экспериментальную модель</h2>
          <div className="section-copy">В отличие от Атласа опухолей, эта страница начинается не с заболевания, а с самой модели. Одна линия показывается один раз; её участие в уже настроенных исследованиях MCL отмечено отдельно.</div>
        </div>
      </div>
      <CRISPRModelTable models={models.items} initialOrgan={params.organ || ""} initialCancer={params.cancer || ""}/>
    </section>
  </>;
}
