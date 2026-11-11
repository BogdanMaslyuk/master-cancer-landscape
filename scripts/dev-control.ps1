param(
    [ValidateSet('menu', 'start', 'restart', 'stop', 'status', 'open')]
    [string]$Action = 'menu',
    [ValidateSet('all', 'api', 'web')]
    [string]$Service = 'all'
)

$ErrorActionPreference = 'Stop'
$Root = [System.IO.Path]::GetFullPath((Join-Path $PSScriptRoot '..'))
$Worker = Join-Path $PSScriptRoot 'dev-service.ps1'
$StateDir = Join-Path $Root '.mcl-dev'
$Settings = @{
    api = @{ Port = 8100; Label = 'PS1 API' }
    web = @{ Port = 3100; Label = 'PS2 SITE' }
}

function Get-RecordPath([string]$Name) {
    return Join-Path $StateDir ($Name + '.json')
}

function Read-Record([string]$Name) {
    $Path = Get-RecordPath $Name
    if (-not (Test-Path -LiteralPath $Path -PathType Leaf)) { return $null }
    try { return (Get-Content -LiteralPath $Path -Raw -Encoding UTF8 | ConvertFrom-Json) }
    catch {
        Write-Host "Cannot read $Path; refusing to control an unverified process." -ForegroundColor Red
        return $null
    }
}

function Find-ManagedProcess([string]$Name) {
    $Record = Read-Record $Name
    if ($null -eq $Record) { return $null }
    $ProcessId = 0
    if (-not [int]::TryParse([string]$Record.process_id, [ref]$ProcessId) -or $ProcessId -le 0) { return $null }
    $Proc = Get-CimInstance -ClassName Win32_Process -Filter "ProcessId = $ProcessId" -ErrorAction SilentlyContinue
    if ($null -eq $Proc) { return $null }
    $Command = [string]$Proc.CommandLine
    if ($Command.IndexOf($Worker, [System.StringComparison]::OrdinalIgnoreCase) -lt 0) { return $null }
    if ($Command -notmatch ('(?i)-Service\s+' + [regex]::Escape($Name) + '(\s|$)')) { return $null }
    $Started = $Proc.CreationDate.ToUniversalTime().ToString('o')
    if ($Started -ne [string]$Record.created_utc) { return $null }
    return $Proc
}

function Port-Listener([int]$Port) {
    $Rows = @(Get-NetTCPConnection -LocalPort $Port -State Listen -ErrorAction SilentlyContinue)
    if ($Rows.Count -eq 0) { return $null }
    return $Rows[0]
}

function Start-Managed([string]$Name) {
    $Info = $Settings[$Name]
    if ($null -ne (Find-ManagedProcess $Name)) {
        Write-Host "$($Info.Label): managed window already exists. Use 'restart' to relaunch." -ForegroundColor Yellow
        return
    }
    $Listener = Port-Listener $Info.Port
    if ($null -ne $Listener) {
        Write-Host "$($Info.Label): port $($Info.Port) is occupied by PID $($Listener.OwningProcess)." -ForegroundColor Red
        Write-Host 'This process was NOT started by this manager, so it will NOT be stopped.' -ForegroundColor Yellow
        Write-Host 'Stop the old manual PS window with Ctrl+C, then retry.' -ForegroundColor Yellow
        return
    }
    if (-not (Test-Path -LiteralPath $Worker -PathType Leaf)) { throw "Missing worker: $Worker" }
    New-Item -ItemType Directory -Path $StateDir -Force | Out-Null
    $Arguments = '-NoProfile -NoExit -ExecutionPolicy Bypass -File "{0}" -Service {1}' -f $Worker, $Name
    $Started = Start-Process -FilePath 'powershell.exe' -WorkingDirectory $Root -ArgumentList $Arguments -PassThru
    $Proc = $null
    for ($i = 0; $i -lt 15; $i++) {
        Start-Sleep -Milliseconds 150
        $Proc = Get-CimInstance -ClassName Win32_Process -Filter "ProcessId = $($Started.Id)" -ErrorAction SilentlyContinue
        if ($null -ne $Proc -and $Proc.CommandLine) { break }
    }
    if ($null -eq $Proc -or -not $Proc.CommandLine) {
        Write-Host 'Window was opened but could not be registered safely.' -ForegroundColor Red
        Write-Host 'Close the newly opened window manually and try again.' -ForegroundColor Yellow
        return
    }
    $Creation = $Proc.CreationDate.ToUniversalTime().ToString('o')
    $Record = [ordered]@{
        service = $Name
        process_id = [int]$Started.Id
        created_utc = $Creation
        working_root = $Root
        port = [int]$Info.Port
    }
    $Record | ConvertTo-Json | Set-Content -LiteralPath (Get-RecordPath $Name) -Encoding UTF8
    Write-Host "$($Info.Label): opened separate PowerShell window (PID $($Started.Id), port $($Info.Port))." -ForegroundColor Green
}

function Stop-Managed([string]$Name) {
    $Info = $Settings[$Name]
    $Proc = Find-ManagedProcess $Name
    if ($null -eq $Proc) {
        if (Test-Path -LiteralPath (Get-RecordPath $Name)) {
            Write-Host "$($Info.Label): stale/unknown registration; no process will be killed." -ForegroundColor Yellow
            Remove-Item -LiteralPath (Get-RecordPath $Name) -Force
        }
        else { Write-Host "$($Info.Label): not managed, nothing to stop." -ForegroundColor DarkGray }
        return
    }
    $ProcessId = [int]$Proc.ProcessId
    Write-Host "$($Info.Label): stopping verified window PID $ProcessId and its own process tree..." -ForegroundColor Cyan
    & taskkill.exe /PID $ProcessId /T /F | Out-Host
    if ($LASTEXITCODE -ne 0) { throw "Could not stop $Name: taskkill returned $LASTEXITCODE" }
    Remove-Item -LiteralPath (Get-RecordPath $Name) -Force -ErrorAction SilentlyContinue
}

function Show-Status {
    Write-Host ''
    Write-Host '=== MCL PYZ 2.0 | STATUS ===' -ForegroundColor Cyan
    Write-Host "Working folder: $Root" -ForegroundColor Gray
    foreach ($Name in @('api','web')) {
        $Info = $Settings[$Name]
        $Proc = Find-ManagedProcess $Name
        $Listener = Port-Listener $Info.Port
        if ($null -ne $Proc -and $null -ne $Listener) {
            Write-Host "$($Info.Label) : window PID $($Proc.ProcessId) | port $($Info.Port) LISTENING" -ForegroundColor Green
        }
        elseif ($null -ne $Proc) {
            Write-Host "$($Info.Label) : window open, port not listening (starting or failed). Check its logs." -ForegroundColor Yellow
        }
        elseif ($null -ne $Listener) {
            Write-Host "$($Info.Label) : port $($Info.Port) occupied by UNMANAGED PID $($Listener.OwningProcess)" -ForegroundColor Yellow
        }
        else {
            Write-Host "$($Info.Label) : stopped" -ForegroundColor DarkGray
        }
    }
    Write-Host 'PYZ: http://127.0.0.1:3100/pyz/PYZ-035' -ForegroundColor Cyan
}

function Chosen-Services([string]$Name) {
    if ($Name -eq 'all') { return @('api','web') }
    return @($Name)
}

function Invoke-Action([string]$Command, [string]$Scope) {
    $Names = @(Chosen-Services $Scope)
    switch ($Command) {
        'start' { foreach ($Name in $Names) { Start-Managed $Name } }
        'stop' { foreach ($Name in $Names) { Stop-Managed $Name } }
        'restart' {
            foreach ($Name in $Names) { Stop-Managed $Name }
            Start-Sleep -Milliseconds 500
            foreach ($Name in $Names) { Start-Managed $Name }
        }
        'status' { Show-Status }
        'open' { Start-Process 'http://127.0.0.1:3100/pyz/PYZ-035' }
    }
}

if ($Action -ne 'menu') {
    Invoke-Action $Action $Service
    return
}

Write-Host 'MCL PYZ 2.0 | LOCAL CONTROL' -ForegroundColor Cyan
Write-Host 'Only manager-owned windows may be stopped. Old MCL is protected.' -ForegroundColor Gray
while ($true) {
    Write-Host ''
    Write-Host '[1] Start PS1 + PS2       [2] RESTART BOTH' -ForegroundColor White
    Write-Host '[3] Restart PS1 (API)    [4] Restart PS2 (site)' -ForegroundColor White
    Write-Host '[5] Status               [6] Open PYZ-035' -ForegroundColor White
    Write-Host '[7] Stop both            [0] Exit control (keep servers running)' -ForegroundColor White
    $Choice = (Read-Host 'Select').Trim()
    try {
        switch ($Choice) {
            '1' { Invoke-Action 'start' 'all' }
            '2' { Invoke-Action 'restart' 'all' }
            '3' { Invoke-Action 'restart' 'api' }
            '4' { Invoke-Action 'restart' 'web' }
            '5' { Invoke-Action 'status' 'all' }
            '6' { Invoke-Action 'open' 'all' }
            '7' { Invoke-Action 'stop' 'all' }
            '0' { return }
            default { Write-Host 'Unknown option.' -ForegroundColor Yellow }
        }
    }
    catch {
        Write-Host "Error: $($_.Exception.Message)" -ForegroundColor Red
    }
}
