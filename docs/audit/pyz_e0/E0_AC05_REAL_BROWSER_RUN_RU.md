# MCL / Э0-AC-05 — единый запуск

Открыть PowerShell 5.1, убедиться, что сервисы PS1 (API 8100) и PS2 (сайт 3100) уже запущены из центрального `MCL_Project`. Одна команда для существующего корня и отдельной ветки:

```powershell
$ErrorActionPreference = 'Stop'
$Root = "$env:USERPROFILE\MCL_Project"
$Repo = "$Root\repositories\master-cancer-landscape"
$Review = "$Root\repositories\mcl-pyz-e0-ac05-review"
$Branch = 'feature/pyz-e0-ac05-browser-verification'
git -C $Repo fetch origin "refs/heads/${Branch}:refs/remotes/origin/${Branch}"
if ($LASTEXITCODE -ne 0) { throw 'Fetch failed' }
if (-not (Test-Path $Review)) {
    git -C $Repo worktree add --detach $Review "refs/remotes/origin/$Branch"
} else {
    $Changed = @(git -C $Review status --porcelain)
    if ($LASTEXITCODE -ne 0 -or $Changed.Count -ne 0) { throw 'Worktree has local changes; stopped safely' }
    git -C $Review switch --detach "refs/remotes/origin/$Branch"
}
if ($LASTEXITCODE -ne 0) { throw 'Worktree update failed' }
powershell.exe -NoProfile -ExecutionPolicy Bypass -File "$Review\scripts\run-pyz-e0-ac05-browser.ps1" -Root $Root
if ($LASTEXITCODE -ne 0) { throw 'AC-05 browser audit failed; see report and console error' }
```

Если Playwright ещё не установлен, программа установит только библиотеку `playwright-core` в `reports/_tools/ac05`, не скачивая Chrome/Chromium. Для запуска нужен системный Edge или Chrome; если ни один не открывается в автоматическом режиме, отчёт должен содержать `NOT_RUN_BROWSER_UNAVAILABLE`, а не выдавать ложный успех.

После проверки прикрепить полученный `MCL_E0_AC05_BROWSER_*.zip`. Ни Э0, ни AC-05 не утверждаются автоматически. Никакие scientific PYZ файлы не перезаписываются.