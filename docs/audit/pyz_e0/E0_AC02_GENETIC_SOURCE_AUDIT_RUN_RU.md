# Э0 / AC-02 — запуск научного аудита из MCL_Project

**Не требуется** включать PS1/PS2, открывать браузер, скачивать DepMap, запускать `build-explorer.ps1` или менять исходники. Достаточно уже существующих локальных папок и Python из `.venv`.

Откройте PowerShell и выполните **один блок**:

```powershell
$ErrorActionPreference = 'Stop'
$Root = "$env:USERPROFILE\MCL_Project"
$Repo = "$Root\repositories\master-cancer-landscape"
$Review = "$Root\repositories\mcl-pyz-e0-ac02-review"
$Branch = 'feature/pyz-e0-ac02-genetic-source-audit'

git -C $Repo fetch origin "refs/heads/${Branch}:refs/remotes/origin/${Branch}"
if ($LASTEXITCODE -ne 0) { throw 'Не удалось скачать ветку AC-02' }

if (-not (Test-Path $Review)) {
    git -C $Repo worktree add --detach $Review "refs/remotes/origin/$Branch"
} else {
    $Changed = @(git -C $Review status --porcelain)
    if ($LASTEXITCODE -ne 0 -or $Changed.Count -gt 0) {
        throw 'Есть локальные изменения в проверяемой копии; не перезаписываем их'
    }
    git -C $Review switch --detach "refs/remotes/origin/$Branch"
}
if ($LASTEXITCODE -ne 0) { throw 'Не удалось подготовить копию AC-02' }

powershell.exe -NoProfile -ExecutionPolicy Bypass `
    -File "$Review\scripts\run-pyz-e0-ac02.ps1" -Root $Root
if ($LASTEXITCODE -ne 0) { throw 'Проверка AC-02 завершилась ошибкой запуска' }
```

Результат появится в `C:\Users\Bogdan\MCL_Project\reports\MCL_E0_AC02_GENETIC_<дата_время>.zip`. Если анализатор выявит научный блокер, он **всё равно создаст отчёт**, а затем предупредит о проблеме. Передайте ZIP для разбора; это не основание автоматически перестраивать данные.

Чтение `gene_effect.npy`, `expression.npy`, `copy_number.npy` идёт в режиме memory-mapped read-only. Parquet используется только для проверки осей и небольшой выборки. Все исходники сохраняются в первоначальном виде. **Э0 не принимается автоматически.**