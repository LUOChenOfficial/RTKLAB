param([string]$Baseline = 'e16b799')
$ErrorActionPreference = 'Stop'
$repo = (Resolve-Path (Join-Path $PSScriptRoot '..\..')).Path
$taskFolder = Join-Path $repo 'tmp\ransac_validation'
New-Item -ItemType Directory -Force -Path $taskFolder | Out-Null
$vcvars = 'C:\Program Files\Microsoft Visual Studio\18\Community\VC\Auxiliary\Build\vcvars32.bat'
$compilerEnv = & cmd.exe /c ('call "' + $vcvars + '" >nul && set')
if ($LASTEXITCODE) { throw 'vcvars32 failed' }
foreach ($line in $compilerEnv) {
    if ($line -match '^([^=]+)=(.*)$') {
        [Environment]::SetEnvironmentVariable($matches[1], $matches[2], 'Process')
    }
}
$originalCL = $env:CL
Push-Location $taskFolder
try {
    $archive = Join-Path $taskFolder 'baseline.zip'
    & git -C $repo archive --format=zip "--output=$archive" $Baseline code
    if ($LASTEXITCODE) { throw 'baseline archive failed' }
    $baselineFolder = Join-Path $taskFolder 'baseline_source'
    Expand-Archive -LiteralPath $archive -DestinationPath $baselineFolder -Force
    foreach ($variant in @('baseline', 'off', 'on', 'slope')) {
        $source = if ($variant -eq 'baseline') { $baselineFolder } else { $repo }
        $project = Join-Path $source 'code\rtklibdemo5\rtklibdemo5.vcxproj'
        $buildFolder = Join-Path $taskFolder $variant
        New-Item -ItemType Directory -Force -Path $buildFolder | Out-Null
        $enabled = if ($variant -in @('on','slope')) { 1 } else { 0 }
        $combined = if ($variant -eq 'slope') { 1 } else { 0 }
        $env:CL = "/DENABLE_RTK_ONLINE_MC=0 /DENABLE_RTK_RANSAC=$enabled /DENABLE_RTK_RANSAC_SLOPE=$combined /DENABLE_RTK_RANSAC_ARAIM=0"
        & msbuild $project /t:Rebuild /p:Configuration=Debug /p:Platform=Win32 "/p:IntDir=$buildFolder\" "/p:OutDir=$buildFolder\" /verbosity:quiet /nologo *> (Join-Path $buildFolder 'build.log')
        if ($LASTEXITCODE) { Get-Content (Join-Path $buildFolder 'build.log') -Tail 20; throw "$variant build failed" }
        $objects = @(Get-ChildItem -LiteralPath $buildFolder -Filter '*.obj' | Where-Object { $_.Name -notin @('main.obj','smoke.obj','kernel.obj','slope_test.obj','callbacks.obj') } | ForEach-Object FullName)
        $include = Join-Path $source 'code\rtklibdemo5\src'
        & cl /nologo /utf-8 /MDd /D_DEBUG "/I$include" "/Fo$buildFolder\callbacks.obj" /c (Join-Path $PSScriptRoot 'test_callbacks.c')
        if ($LASTEXITCODE) { throw 'callbacks compile failed' }
        $objects += Join-Path $buildFolder 'callbacks.obj'
        $smoke = Join-Path $buildFolder 'smoke.exe'
        & cl /nologo /utf-8 /MDd /D_DEBUG "/I$include" "/Fo$buildFolder\smoke.obj" "/Fe$smoke" (Join-Path $PSScriptRoot 'ransac_smoke.c') @objects /link ws2_32.lib winmm.lib
        if ($LASTEXITCODE) { throw "$variant smoke link failed" }
        & $smoke (Join-Path $buildFolder 'static.pos') 2> (Join-Path $buildFolder 'run.log')
        if ($LASTEXITCODE) { throw "$variant smoke run failed" }
        if ($variant -eq 'on') {
            $kernel = Join-Path $buildFolder 'kernel.exe'
            & cl /nologo /utf-8 /MDd /D_DEBUG "/I$include" "/Fo$buildFolder\kernel.obj" "/Fe$kernel" (Join-Path $PSScriptRoot 'ransac_test.c') @objects /link ws2_32.lib winmm.lib
            if ($LASTEXITCODE) { throw 'kernel link failed' }
            & $kernel
            if ($LASTEXITCODE) { throw 'kernel assertions failed' }
        }
        if ($variant -eq 'slope') {
            $slopeTest = Join-Path $buildFolder 'slope_test.exe'
            & cl /nologo /utf-8 /MDd /D_DEBUG "/I$include" "/Fo$buildFolder\slope_test.obj" "/Fe$slopeTest" (Join-Path $PSScriptRoot 'slope_test.c') @objects /link ws2_32.lib winmm.lib
            if ($LASTEXITCODE) { throw 'slope test link failed' }
            & $slopeTest
            if ($LASTEXITCODE) { throw 'slope assertions failed' }
        }
    }
    $baseRows = @(Get-Content (Join-Path $taskFolder 'baseline\static.pos') | Where-Object { $_ -and $_ -notmatch '^%' })
    $offRows = @(Get-Content (Join-Path $taskFolder 'off\static.pos') | Where-Object { $_ -and $_ -notmatch '^%' })
    if (!$baseRows.Count -or (Compare-Object $baseRows $offRows -SyncWindow 0)) { throw 'disabled regression differs or is empty' }
    $onRows = @(Get-Content (Join-Path $taskFolder 'on\static.pos') | Where-Object { $_ -and $_ -notmatch '^%' })
    if (!$onRows.Count) { throw 'enabled smoke output is empty' }
    $records = @(Get-Content (Join-Path $taskFolder 'on\static.rtk_debug') | Where-Object { $_ -match '^RANSAC_SD ' })
    if (!$records.Count) { throw 'no RANSAC records' }
    $rejected = @($records | Where-Object { ($_ -split '\s+')[8] -eq '0' })
    $mask = @{}
    $checked = 0
    foreach ($line in Get-Content (Join-Path $taskFolder 'on\static.rtk_debug')) {
        $a = $line -split '\s+'
        if ($a[0] -eq 'RANSAC_SD') { $mask["$($a[1]) $($a[2]) $($a[3]) $($a[4])"] = $a[8] }
        if ($a[0] -eq 'AMB_FIX') {
            foreach ($satId in @($a[7],$a[8])) {
                $key = "$($a[1]) $($a[2]) $satId $($a[9])"
                $checked++
                if ($mask.ContainsKey($key) -and $mask[$key] -eq '0') { throw "rejected carrier used in AR: $key" }
            }
        }
    }
    "PASS: $checked fixed-ambiguity endpoints contain no rejected carriers"
    "PASS: disabled matches baseline ($($baseRows.Count) epochs); enabled $($onRows.Count) epochs, $($records.Count) SD records, $($rejected.Count) rejected"
    $slopeRows = @(Get-Content (Join-Path $taskFolder 'slope\static.pos') | Where-Object { $_ -and $_ -notmatch '^%' })
    if (!$slopeRows.Count -or @($slopeRows | Where-Object { ($_ -split '\s+')[5] -ne '2' }).Count) { throw 'slope mode must output float estimates' }
    $slopeRecords = @(Get-Content (Join-Path $taskFolder 'slope\static.rtk_debug') | Where-Object { $_ -match '^SLOPE_KF ' })
    $available = @($slopeRecords | Where-Object { ($_ -split '\s+')[3] -eq '1' })
    if (!$available.Count) { throw 'no available slope PL in smoke test' }
    "PASS: slope $($slopeRows.Count) float epochs, $($available.Count)/$($slopeRecords.Count) available PL"
}
finally {
    $env:CL = $originalCL
    Pop-Location
}
