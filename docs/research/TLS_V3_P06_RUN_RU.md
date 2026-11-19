# Как запускать TLS V3/P0.6 (Windows PowerShell 5.1)

Проект `C:\Users\Bogdan\MCL_Project`, головной Git-репозиторий `repositories\master-cancer-landscape`. Изменения подключаются отдельной копией `repositories\mcl-tls-v3-p06-review` и **не** сливаются с рабочей веткой.

```powershell
$ErrorActionPreference='Stop'
$Root="$env:USERPROFILE\MCL_Project"
$Repo="$Root\repositories\master-cancer-landscape"
$Review="$Root\repositories\mcl-tls-v3-p06-review"
$Branch='feature/tls-v3-p06-bounded-external-pilot'

git -C $Repo fetch origin "refs/heads/${Branch}:refs/remotes/origin/${Branch}"
if ($LASTEXITCODE -ne 0) { throw 'Git fetch failed' }
if (-not (Test-Path $Review)) {
    git -C $Repo worktree add --detach $Review "refs/remotes/origin/$Branch"
} else {
    $Changes=@(git -C $Review status --porcelain)
    if ($LASTEXITCODE -ne 0 -or $Changes.Count -gt 0) { throw 'Review worktree has changes' }
    git -C $Review switch --detach "refs/remotes/origin/$Branch"
}
if ($LASTEXITCODE -ne 0) { throw 'Worktree failed' }

# Один запуск: автоматически пересчитает исправленный P0 без сети,
# затем построит воспроизводимую выборку P0.6 без сети.
powershell.exe -NoProfile -ExecutionPolicy Bypass `
    -File "$Review\scripts\run-target-space-v3-p06.ps1" -Root $Root
if ($LASTEXITCODE -ne 0) { throw 'P0.6 offline failed' }
```

Для **разрешённого внешнего пилота** по UniProt + ChEMBL + BindingDB (по умолчанию 40 мишеней) выполнить эту же команду с явным параметром `-EnableExternalPilot`. Данный шаг делает сетевые обращения с вашего компьютера, может столкнуться с квотами/недоступностью источников и сохраняет все ошибки в отчёте.

Объекты вывода: `TLS_V3_P0_*.zip` и `TLS_V3_P06_*.zip` в `MCL_Project\reports\`; при ошибке `TLS_V3_P06_DIAGNOSTIC_*.log`. Пришлите оба ZIP для научной проверки. Никакого автоматического одобрения новых мишеней нет.