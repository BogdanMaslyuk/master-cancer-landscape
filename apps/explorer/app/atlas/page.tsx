import Link from "next/link";
import AtlasIcon from "../../components/AtlasIcon";
import { apiGet } from "../../lib/api";

export const dynamic = "force-dynamic";

type Cancer = {
  id: string;
  name: string;
  models_n: number;
  subtypes_n: number;
  subtypes: { id: string; name: string; models_n: number }[];
  curated_context_ids: string[];
};
type Organ = {
  id: string;
  system_ru: string;
  name_ru: string;
  name_en?: string | null;
  icon: string;
  models_n: number;
  cancers_n: number;
  cancers: Cancer[];
};
type Atlas = {
  status: string;
  canonical_index_available: boolean;
  build_command?: string | null;
  note_ru?: string;
  depmap_release?: string | null;
  models_n: number;
  organs_n: number;
  cancers_n: number;
  subtypes_n: number;
  requires_review_n: number;
  organs: Organ[];
};

function modelHref(organId: string, cancerId: string) {
  return `/models?organ=${encodeURIComponent(organId)}&cancer=${encodeURIComponent(cancerId)}`;
}

export default async function AtlasPage() {
  const atlas = await apiGet<Atlas>("/api/crispr-atlas");
  return <>
    <section className="page-intro atlas-intro">
      <div>
        <div className="eyebrow">CRISPR-АТЛАС ЭКСПЕРИМЕНТАЛЬНЫХ МОДЕЛЕЙ</div>
        <h1>Атлас опухолей</h1>
        <p>Атлас строится снизу вверх из реальных клеточных моделей, для которых доступен CRISPR Gene Effect. Выберите орган, затем опухоль — и перейдите к экспериментальным линиям, на которых можно исследовать генетические зависимости.</p>
      </div>
      <div className="atlas-route-card">
        <span>Путь исследования</span>
        <strong>Орган / система → опухоль → подтип → клеточные модели → CRISPR-зависимости</strong>
      </div>
    </section>

    <section className="kpi-strip atlas-kpis">
      <div className="kpi"><strong>{atlas.organs_n}</strong><span>органов / систем представлены</span></div>
      <div className="kpi"><strong>{atlas.cancers_n}</strong><span>групп опухолей</span></div>
      <div className="kpi"><strong>{atlas.subtypes_n}</strong><span>подтипов опухолей</span></div>
      <div className="kpi"><strong>{atlas.models_n}</strong><span>уникальных CRISPR-моделей</span></div>
    </section>

    {!atlas.canonical_index_available && <section className="section callout atlas-callout">
      <div className="eyebrow">ВРЕМЕННЫЙ РЕЖИМ</div>
      <h3>Полный CRISPR-атлас ещё не материализован</h3>
      <p className="section-copy">{atlas.note_ru} После запуска команды ниже страница автоматически перестроится по всем строкам CRISPRGeneEffect.csv и аннотациям Model.csv.</p>
      {atlas.build_command && <code>{atlas.build_command}</code>}
    </section>}

    <section className="section callout atlas-callout">
      <div className="eyebrow">ПРАВИЛО КЛАССИФИКАЦИИ</div>
      <h3>Опухоль определяется по OncoTree, а не по месту получения образца</h3>
      <p className="section-copy">Если клеточную линию получили из метастаза, например в печени, она остаётся в разделе первичной опухоли. Поля «первичная опухоль / метастаз» и место взятия образца сохраняются в карточке модели, но не меняют её положение в Атласе.</p>
    </section>

    <section className="section">
      <div className="section-header">
        <div>
          <div className="eyebrow">НАВИГАЦИЯ ПО ЗАБОЛЕВАНИЯМ</div>
          <h2>Выберите орган или систему</h2>
          <div className="section-copy">Внутри каждого раздела показаны опухоли, для которых реально существуют CRISPR-профилированные модели в используемом релизе DepMap{atlas.depmap_release ? ` ${atlas.depmap_release}` : ""}. Числа отражают модели, а не частоту заболеваний у пациентов.</div>
        </div>
      </div>

      <div className="organ-grid">
        {atlas.organs.map((organ) => <article className="organ-card" key={organ.id}>
          <div className="organ-card-head">
            <div className={`organ-icon ${organ.icon}`}><AtlasIcon name={organ.icon}/></div>
            <div>
              <div className="eyebrow">{organ.system_ru}</div>
              <h3>{organ.name_ru}</h3>
              {organ.name_en && organ.name_en !== organ.name_ru && <div className="technical-note">OncoTree lineage: {organ.name_en}</div>}
            </div>
          </div>
          <div className="organ-metrics">
            <span><b>{organ.cancers_n}</b> опухолей</span>
            <span><b>{organ.models_n}</b> CRISPR-моделей</span>
          </div>
          <div className="organ-contexts">
            {organ.cancers.slice(0, 12).map((cancer) => <Link href={modelHref(organ.id, cancer.id)} className="organ-context-link" key={cancer.id}>
              <div>
                <strong>{cancer.name}</strong>
                <span>{cancer.subtypes_n ? `${cancer.subtypes_n} подтипов` : "Подтип не указан"}</span>
              </div>
              <div className="context-link-meta">
                <small>{cancer.models_n} клеточных моделей</small>
                {cancer.curated_context_ids.length > 0 && <span className="analysis-state ready">есть анализ MCL</span>}
              </div>
              <b className="context-arrow">→</b>
            </Link>)}
            {organ.cancers.length > 12 && <Link href={`/models?organ=${encodeURIComponent(organ.id)}`} className="organ-context-link">
              <div><strong>Показать все опухоли раздела</strong><span>Ещё {organ.cancers.length - 12}</span></div>
              <b className="context-arrow">→</b>
            </Link>}
          </div>
        </article>)}
      </div>
    </section>

    <section className="section callout atlas-callout">
      <div className="eyebrow">ДВА РАЗНЫХ ИНСТРУМЕНТА</div>
      <h3>Атлас отвечает «где искать», Исследователь моделей — «какую линию изучать»</h3>
      <p className="section-copy">Атлас организует CRISPR-модели по заболеванию. После выбора опухоли вы переходите на страницу «Исследователь моделей», где можно искать конкретные линии, фильтровать их и открывать индивидуальную генетику и функциональные данные.</p>
    </section>
  </>;
}
