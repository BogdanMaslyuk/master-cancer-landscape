import Link from "next/link";
import styles from "./profile-v2.module.css";

export default function NotFoundPyzCard() {
  return <main className={styles.page}>
    <div className={styles.identity}>
      <span className={styles.eyebrow}>404 · НЕИЗВЕСТНАЯ МОЛЕКУЛА</span>
      <h1 style={{ fontSize: 30, margin: "10px 0" }}>Соединение отсутствует в реестре PYZ</h1>
      <p className={styles.description}>Номер должен соответствовать зарегистрированному соединению PYZ-001…PYZ-065. Мы не подставляем вымышленные результаты.</p>
      <Link className={styles.primaryAction} href="/pyz">Открыть каталог PYZ</Link>
    </div>
  </main>;
}
