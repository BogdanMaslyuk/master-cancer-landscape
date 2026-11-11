"use client";

import { useState } from "react";

type Props = { smiles: string | null; compoundId: string };

export function CopySmiles({ smiles, compoundId }: Props) {
  const [feedback, setFeedback] = useState("");
  async function copy() {
    if (!smiles) return;
    try {
      await navigator.clipboard.writeText(smiles);
      setFeedback("SMILES скопирован");
    } catch {
      setFeedback("Не удалось скопировать — выделите SMILES вручную");
    }
  }
  return (
    <span style={{ display: "inline-flex", alignItems: "center", flexWrap: "wrap", gap: 10 }}>
      <button type="button" onClick={copy} disabled={!smiles} aria-label={"Скопировать SMILES " + compoundId}
        style={{ minHeight: 44, border: "1px solid #b8cbe1", borderRadius: 10, background: "#edf5fd", padding: "0 14px", fontWeight: 700, color: "#1c518c", cursor: smiles ? "pointer" : "not-allowed" }}>
        Скопировать SMILES
      </button>
      <span role="status" aria-live="polite" style={{ color: "#50667a", fontSize: 12 }}>{feedback}</span>
    </span>
  );
}

export function StructurePreview({ compoundId, src }: { compoundId: string; src: string }) {
  const [failed, setFailed] = useState(false);
  if (failed) {
    return (
      <div role="img" aria-label={"Двумерная структура " + compoundId + " недоступна"}
        style={{ padding: 24, textAlign: "center", color: "#52667b", lineHeight: 1.6 }}>
        Изображение структуры недоступно. Проверьте SMILES ниже и работу локального API.
      </div>
    );
  }
  // API renders the structural SVG from the exact registered compound SMILES.
  // eslint-disable-next-line @next/next/no-img-element
  return <img src={src} alt={"Двумерная химическая структура " + compoundId}
    onError={() => setFailed(true)}
    style={{ display: "block", width: "100%", height: "100%", maxHeight: 320, objectFit: "contain" }} />;
}
