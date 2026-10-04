"use client";

import { FormEvent, useEffect, useMemo, useState } from "react";
import { useRouter } from "next/navigation";
import styles from "./GeneSearch.module.css";

type Suggestion = {
  gene_symbol: string;
  best_context_label?: string | null;
  best_delta_gene_effect?: number | null;
  best_model_gene_effect?: number | null;
  present_all_thresholds?: boolean | null;
};

const API_BASE = process.env.NEXT_PUBLIC_MCL_API_URL || "http://127.0.0.1:8000";

function cleanSymbol(value: string) {
  return value.trim().toUpperCase().replace(/\s+/g, "");
}

function smallNumber(value: unknown) {
  const n = Number(value);
  return Number.isFinite(n) ? n.toFixed(2) : "—";
}

export default function GeneSearch({ compact = false }: { compact?: boolean }) {
  const router = useRouter();
  const [value, setValue] = useState("");
  const [suggestions, setSuggestions] = useState<Suggestion[]>([]);
  const [focused, setFocused] = useState(false);
  const symbol = useMemo(() => cleanSymbol(value), [value]);

  useEffect(() => {
    if (symbol.length < 1) {
      setSuggestions([]);
      return;
    }
    const controller = new AbortController();
    const timer = window.setTimeout(async () => {
      try {
        const response = await fetch(`${API_BASE}/api/genes/suggest?q=${encodeURIComponent(symbol)}&limit=8`, {
          signal: controller.signal,
          cache: "no-store",
        });
        if (!response.ok) return;
        const payload = (await response.json()) as Suggestion[];
        setSuggestions(Array.isArray(payload) ? payload : []);
      } catch {
        // Autocomplete is optional; direct submit still works if the API is warming up.
      }
    }, 140);
    return () => {
      window.clearTimeout(timer);
      controller.abort();
    };
  }, [symbol]);

  function openGene(next: string) {
    const cleaned = cleanSymbol(next);
    if (!cleaned) return;
    setFocused(false);
    setSuggestions([]);
    router.push(`/genes/${encodeURIComponent(cleaned)}`);
  }

  function submit(event: FormEvent) {
    event.preventDefault();
    openGene(value);
  }

  return (
    <div className={`${styles.root} ${compact ? styles.compact : ""}`}>
      <form className={styles.form} onSubmit={submit} role="search">
        <span className={styles.icon} aria-hidden="true">⌕</span>
        <input
          value={value}
          onChange={(event) => setValue(event.target.value)}
          onFocus={() => setFocused(true)}
          onBlur={() => window.setTimeout(() => setFocused(false), 120)}
          placeholder="Найти ген: AHR, KRAS, DNM1L…"
          aria-label="Найти ген или белок"
          autoComplete="off"
          spellCheck={false}
        />
        {value && <button type="button" className={styles.clear} onClick={() => setValue("")} aria-label="Очистить поиск">×</button>}
        <button type="submit" className={styles.submit}>Найти</button>
      </form>

      {focused && suggestions.length > 0 && (
        <div className={styles.menu} role="listbox">
          {suggestions.map((item) => (
            <button
              type="button"
              className={styles.option}
              key={item.gene_symbol}
              onMouseDown={(event) => event.preventDefault()}
              onClick={() => openGene(item.gene_symbol)}
            >
              <span className={styles.symbol}>{item.gene_symbol}</span>
              <span className={styles.meta}>
                {item.present_all_thresholds ? <b>устойчивый</b> : <span>{item.best_context_label || "геномный набор DepMap"}</span>}
                <small>ΔGE {smallNumber(item.best_delta_gene_effect)} · min GE {smallNumber(item.best_model_gene_effect)}</small>
              </span>
            </button>
          ))}
        </div>
      )}
    </div>
  );
}
