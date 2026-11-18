# MCL TLS V3 / P0 — запуск локальной инвентаризации

**Не изменять научные файлы.** Никакой интернет-загрузки, пересборки старой матрицы 65 × 13 или слияния PR.

Команда Windows PowerShell 5.1 (вставить целиком):

```powershell
$ErrorActionPreference = 'Stop'
$Root = "$env:USERPROFILE\MCL_Project"
$Repo = "$Root\repositories\master-cancer-landscape"
$Review = "$Root\repositories\mcl-tls-v3-p0-review"
$Branch = 'feature/tls-v3-p0-target-universe-inventory'

git -C $Repo fetch origin "refs/heads/${Branch}:refs/remotes/origin/${Branch}"
if ($LASTEXITCODE -ne 0) { throw 'Git fetch failed' }
if (-not (Test-Path -LiteralPath $Review)) {
    git -C $Repo worktree add --detach $Review "refs/remotes/origin/$Branch"
} else {
    $Changes = @(git -C $Review status --porcelain)
    if ($LASTEXITCODE -ne 0 -or $Changes.Count -gt 0) {
        throw 'Review worktree has local changes'
    }
    git -C $Review switch --detach "refs/remotes/origin/$Branch"
}
if ($LASTEXITCODE -ne 0) { throw 'Git worktree update failed' }

powershell.exe -NoProfile -ExecutionPolicy Bypass `
    -File "$Review\scripts\run-target-space-v3-p0.ps1" -Root $Root
if ($LASTEXITCODE -ne 0) { throw 'TLS V3 P0 failed' }
```

Отчёт появится в `MCL_Project/reports/TLS_V3_P0_<timestamp>_<commit>/` и ZIP рядом. Прикрепить ZIP для независимой научной оценки.

Тесты (не заменяют запуск на настоящих данных):

```powershell
& "$Root\repositories\master-cancer-landscape\.venv\Scripts\python.exe" -m pytest -q "$Review\tests\test_target_space_v3_p0.py"
```

**Известные ограничения:** P0 считает локальное покрытие. До проверки полной базы прямых лигандов и полноценного реестра UniProt итоговое число T2/T3/T4 — нижняя оценка, а не полный потенциал исследования. Полное число для P1 пока не утверждается.