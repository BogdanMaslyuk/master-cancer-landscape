import styles from "./profile-v2.module.css";

export default function LoadingPyzCard() {
  return <main className={styles.page} aria-busy="true" aria-label="Загрузка карточки PYZ">
    <div className={styles.pageHeader}><div><span className={styles.eyebrow}>MCL · НАУЧНАЯ КАРТОЧКА 2.0</span>
      <h1>Загружаем молекулу…</h1><p>Читаем зарегистрированную структуру и научные данные.</p></div></div>
    <div className={styles.heroGrid}>
      <div className={styles.structure} style={{ minHeight: 310, background: "#f2f6fa" }}/>
      <div className={styles.identity}><p>Подготавливаем структуру, свойства и ограничения…</p></div>
    </div>
  </main>;
}
