# MCL TLS V3 / P0.8 — запуск в Windows PowerShell 5.1

**Программа не скачивает научные базы самостоятельно.** По умолчанию делает безопасную проверку P0/P0.6/P0.7, формирует контрольную панель из шести белков и показывает, каких локальных источников ещё нет. Работающие матрица 65×13, данные PYZ и настройки сайта не меняются.

## Шаг 1. Получить исходники P0.8

```powershell
$ErrorActionPreference = 'Stop'
$Root = "$env:USERPROFILE\MCL_Project"
$Repo = "$Root\repositories\master-cancer-landscape"
$Review = "$Root\repositories\mcl-tls-v3-p08-review"
$Branch = 'feature/tls-v3-p08-local-experimental-sources'
git -C $Repo fetch origin "refs/heads/${Branch}:refs/remotes/origin/${Branch}"
if ($LASTEXITCODE -ne 0) { throw 'Git fetch failed' }
if (-not (Test-Path -LiteralPath $Review)) {
    git -C $Repo worktree add --detach $Review "refs/remotes/origin/$Branch"
} else {
    $Change = @(git -C $Review status --porcelain)
    if ($LASTEXITCODE -ne 0 -or $Change.Count -ne 0) { throw 'Local review worktree modifications' }
    git -C $Review switch --detach "refs/remotes/origin/$Branch"
}
if ($LASTEXITCODE -ne 0) { throw 'Worktree switch failed' }
```

## Шаг 2. Сначала предварительная проверка — без ChEMBL/BindingDB

```powershell
powershell.exe -NoProfile -ExecutionPolicy Bypass `
    -File "$Review\scripts\run-target-space-v3-p08.ps1" -Root $Root
if ($LASTEXITCODE -ne 0) { throw 'P0.8 preflight failed' }
```

Нужно получить `INPUT_REQUIRED`, пока полноценные выгрузки не указаны.

## Шаг 3. С полными локальными источниками (после загрузки вручную)

Загрузить версионированную SQLite из [официального ChEMBL](https://chembl.gitbook.io/chembl-interface-documentation/downloads) и полную `BindingDB_All_<YYYYMM>_tsv.zip` с [официальной страницы BindingDB](https://www.bindingdb.org/rwd/bind/chemsearch/marvin/Download.jsp). Вставить **реальные пути и реальные метки выпуска** в следующую команду. Папки и имена в образце — только примеры, **не** утверждение о наличии файлов на вашем компьютере.

```powershell
$ChEMBL = 'D:\ScientificData\chembl_36\chembl_36.db'         # заменить на реальный путь
$BindingDB = 'D:\ScientificData\BindingDB_All_202611_tsv.zip' # заменить на реальный путь

powershell.exe -NoProfile -ExecutionPolicy Bypass `
    -File "$Review\scripts\run-target-space-v3-p08.ps1" `
    -Root $Root `
    -ChEMBLSqlite $ChEMBL -ChEMBLRelease 'chembl_36' `
    -BindingDBTsv $BindingDB -BindingDBRelease '202611' `
    -Panel focus6 -MaxRowsPerTarget 1000
if ($LASTEXITCODE -ne 0) { throw 'P0.8 source audit failed' }
```

Результат: `MCL_Project\reports\TLS_V3_P08_<timestamp>_<commit>.zip`.

- Если отсутствует только одна база, можно передать только другой источник и получить частичный отчёт со статусом `INPUT_REQUIRED`.
- Не распаковывать многогигабайтный BindingDB ZIP ради P0.8: приложение читает TSV из ZIP потоком.
- Не запускать второй раз для исправления результата без анализа QC. Ни числа ligand hits, ни количество `REVIEW` строк не равны достоверным научным открытиям.
- При желании проверить все 40 генов исходного пилота можно указать `-Panel sample40`, но до ревизии первой шестёрки это преждевременно.

Пришлите ZIP нового отчёта. До научной приёмки слияние PR и обновление рабочей матрицы запрещены.