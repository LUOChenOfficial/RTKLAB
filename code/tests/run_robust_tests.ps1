$ErrorActionPreference='Stop'
$repo=(Resolve-Path (Join-Path $PSScriptRoot '..\..')).Path
$out=Join-Path $repo 'tmp\robust_kernel_tests'
New-Item -ItemType Directory -Force -Path $out | Out-Null
& 'C:\Program Files\Microsoft Visual Studio\18\Community\MSBuild\Current\Bin\MSBuild.exe' "$repo\code\rtklibdemo5\rtklibdemo5.vcxproj" /t:Build /p:Configuration=Debug /p:Platform=Win32 /verbosity:quiet /nologo
if($LASTEXITCODE){throw 'Debug build failed'}
$compilerEnv=& cmd.exe /c 'call "C:\Program Files\Microsoft Visual Studio\18\Community\VC\Auxiliary\Build\vcvars32.bat" >nul && set'
if($LASTEXITCODE){throw 'compiler environment failed'}
foreach($line in $compilerEnv){if($line -match '^([^=]+)=(.*)$'){[Environment]::SetEnvironmentVariable($matches[1],$matches[2],'Process')}}
$objects=@(Get-ChildItem "$repo\code\rtklibdemo5\rtklibdemo5\Debug" -Filter '*.obj' | Where-Object {$_.Name -ne 'main.obj'} | ForEach-Object FullName)
& cl /nologo /utf-8 /MDd /D_DEBUG "/I$repo\code\rtklibdemo5\src" "/Fo$out\" "/Fe$out\robust_test.exe" "$PSScriptRoot\robust_test.c" "$PSScriptRoot\test_callbacks.c" @objects /link /DEBUG ws2_32.lib winmm.lib
if($LASTEXITCODE){throw 'kernel test build failed'}
& "$out\robust_test.exe"
if($LASTEXITCODE){throw 'kernel assertions failed'}
