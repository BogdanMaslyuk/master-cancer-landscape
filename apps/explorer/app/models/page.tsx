import CellModelTable from "../../components/CellModelTable";
import { apiGet } from "../../lib/api";
import type { ModelsResponse } from "../../lib/generated/api-types";

export default async function ModelsPage(){
  const models=await apiGet<ModelsResponse>("/api/models?limit=5000");
  const unique=new Set(models.items.map((m)=>m.model_id)).size;
  const sequenced=new Set(models.items.filter((m)=>m.sequencing_available).map((m)=>m.model_id)).size;
  const contexts=new Set(models.items.map((m)=>m.cancer_id)).size;

  return <>
    <section className="page-intro compact-intro">
      <div>
        <div className="eyebrow">КАТАЛОГ МОДЕЛЕЙ</div>
        <h1>Клеточные линии</h1>
        <p>Все клеточные модели, которые сейчас участвуют в настроенных опухолевых контекстах MCL. Используйте каталог, чтобы понять, какие экспериментальные системы доступны до интерпретации генетических зависимостей.</p>
      </div>
    </section>

    <section className="kpi-strip">
      <div className="kpi"><strong>{unique}</strong><span>уникальных DepMap-моделей</span></div>
      <div className="kpi"><strong>{contexts}</strong><span>опухолевых контекста</span></div>
      <div className="kpi"><strong>{sequenced}</strong><span>моделей с мутационным профилированием</span></div>
      <div className="kpi"><strong>{models.total}</strong><span>записей «модель × контекст»</span></div>
    </section>

    <section className="section">
      <div className="section-header"><div><div className="eyebrow">ПОИСК И ФИЛЬТРАЦИЯ</div><h2>Найдите клеточную модель</h2><div className="section-copy">Одна линия может встречаться в разных аналитических контекстах. Поэтому каталог показывает не только название, но и заболевание, молекулярный контекст и роль в конкретной выборке.</div></div></div>
      <CellModelTable models={models.items} showContext/>
    </section>
  </>;
}
