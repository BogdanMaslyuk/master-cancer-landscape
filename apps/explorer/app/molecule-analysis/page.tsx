"use client";

import { FormEvent, useState } from "react";
import { useRouter } from "next/navigation";
import { apiBase } from "../../lib/api";
import styles from "./molecule-analysis.module.css";

export default function MoleculeAnalysisPage() {
  const router = useRouter();
  const [smiles, setSmiles] = useState("");
  const [name, setName] = useState("");
  const [busy, setBusy] = useState(false);
  const [error, setError] = useState("");

  async function submit(event: FormEvent) {
    event.preventDefault();
    setError("");
    if (!smiles.trim()) {
      setError("Вставьте SMILES молекулы.");
      return;
    }
    setBusy(true);
    try {
      const response = await fetch(`${apiBase}/api/molecules/analyze`, {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ smiles: smiles.trim(), name: name.trim() || null }),
      });
      if (!response.ok) throw new Error(await response.text());
      const payload = await response.json();
      router.push(`/molecule-analysis/${encodeURIComponent(payload.analysis_id)}`);
    } catch (err) {
      setError(err instanceof Error ? err.message : "Не удалось выполнить анализ.");
    } finally {
      setBusy(false);
    }
  }

  async function loadSmilesFile(file: File | null) {
    if (!file) return;
    const text = await file.text();
    const first = text.split(/\r?\n/).map((line) => line.trim()).find(Boolean);
    if (first) setSmiles(first.split(/\s+/)[0]);
  }

  return (
    <>
      <section className={styles.hero}>
        <div>
          <div className="eyebrow">MCL vNEXT · НОВАЯ МОЛЕКУЛА</div>
          <h1>От SMILES к клеточной модели</h1>
          <p>
            Загрузите новую малую молекулу. MCL стандартизирует структуру, сравнит её с экспериментальными
            лигандами, предложит проверяемые мишени и подберёт наиболее информативные опухолевые модели и контроли.
          </p>
        </div>
        <div className={styles.route}>
          <span>SMILES</span><b>→</b><span>мишени</span><b>→</b><span>опухоли</span><b>→</b><span>культуры</span><b>→</b><span>эксперимент</span>
        </div>
      </section>

      <section className={styles.panel}>
        <div className={styles.panelHead}>
          <div>
            <div className="eyebrow">АНАЛИЗ ОДНОЙ СТРУКТУРЫ</div>
            <h2>Загрузите молекулу</h2>
          </div>
          <div className={styles.status}>MVP · локальный анализ</div>
        </div>
        <form onSubmit={submit} className={styles.form}>
          <label>
            <span>Название · необязательно</span>
            <input value={name} onChange={(e) => setName(e.target.value)} placeholder="Например, соединение 12" />
          </label>
          <label>
            <span>SMILES</span>
            <textarea
              value={smiles}
              onChange={(e) => setSmiles(e.target.value)}
              placeholder="O=C1CC(C)(C)CC(/C1=N\\NC2=NC=CC=C2)=O"
              rows={5}
              spellCheck={false}
            />
          </label>
          <div className={styles.actions}>
            <label className={styles.fileButton}>
              Загрузить .smi / .txt
              <input type="file" accept=".smi,.smiles,.txt" onChange={(e) => loadSmilesFile(e.target.files?.[0] || null)} />
            </label>
            <button className={styles.primary} type="submit" disabled={busy}>{busy ? "Анализируем…" : "Проанализировать молекулу"}</button>
          </div>
          {error && <div className={styles.error}>{error}</div>}
        </form>
      </section>

      <section className={styles.grid}>
        <article><strong>1</strong><h3>Структура</h3><p>RDKit Cleanup → parent structure → canonical SMILES → InChIKey и 2D-формула.</p></article>
        <article><strong>2</strong><h3>Мишени</h3><p>Morgan r2/2048 + Tanimoto против локального Target Ligand Space ChEMBL/BindingDB.</p></article>
        <article><strong>3</strong><h3>Клеточные модели</h3><p>Candidate v2, DepMap CRISPR, PRISM, опухолевый контекст и реестр нашей лаборатории.</p></article>
        <article><strong>4</strong><h3>Эксперимент</h3><p>Положительные модели, отрицательный опухолевый контроль и общий неопухолевый контроль.</p></article>
      </section>

      <div className={styles.guardrail}>
        <b>Научная граница.</b> Система выбирает максимально информативный следующий эксперимент. Структурное сходство, CRISPR и фармакология известных веществ не доказывают активность загруженной молекулы.
      </div>
    </>
  );
}
