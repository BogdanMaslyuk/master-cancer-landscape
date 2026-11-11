"use client";

import Link from "next/link";
import { useEffect } from "react";
import styles from "./profile-v2.module.css";

export default function PyzCardError({ error, reset }: { error: Error & { digest?: string }; reset: () => void }) {
  useEffect(() => {
    // Error details remain in the local browser console; avoid displaying stack traces.
    console.error("PYZ card failed:", error);
  }, [error]);
  return <main className={styles.page}>
    <div className={styles.identity} role="alert">
      <span className={styles.eyebrow}>НЕ УДАЛОСЬ ЗАГРУЗИТЬ КАРТОЧКУ</span>
      <h1 style={{ fontSize: 30, margin: "10px 0" }}>Ошибка данных или подключения к API</h1>
      <p className={styles.description}>Исходные научные данные не изменены. Проверьте, что сервер API запущен и его версия соответствует интерфейсу.</p>
      <div className={styles.actions}>
        <button className={styles.primaryAction} type="button" onClick={() => reset()}>Повторить загрузку</button>
        <Link href="/pyz" className={styles.secondaryAction}>Вернуться к молекулам</Link>
      </div>
    </div>
  </main>;
}
