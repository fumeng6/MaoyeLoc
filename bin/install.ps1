$ErrorActionPreference = 'Stop'
$bundleRoot = $PSScriptRoot
$binDirectory = Join-Path $bundleRoot 'bin'
if (!(Test-Path -LiteralPath (Join-Path $binDirectory 'myloc.exe') -PathType Leaf)) {
    throw 'myloc.exe was not found. Extract the complete bundle before installing.'
}
if (!(Test-Path -LiteralPath (Join-Path $bundleRoot 'myloc.cmd') -PathType Leaf)) {
    throw 'myloc.cmd was not found. Extract the complete bundle before installing.'
}
$userPath = [Environment]::GetEnvironmentVariable('Path', 'User')
$entries = @()
$rootListed = $false
foreach ($entry in ($userPath -split ';')) {
    if (!$entry) { continue }
    $expanded = [Environment]::ExpandEnvironmentVariables($entry.Trim().Trim('"')).TrimEnd('\')
    # Replace this bundle's earlier bin entry; leave other installations alone.
    if ($expanded -ieq $binDirectory.TrimEnd('\')) { continue }
    if ($expanded -ieq $bundleRoot.TrimEnd('\')) { $rootListed = $true }
    $entries += $entry
}
if (!$rootListed) { $entries += $bundleRoot }
[Environment]::SetEnvironmentVariable('Path', ($entries -join ';'), 'User')
Write-Host 'MaoyeLoc: added the myloc command directory to your user PATH.'
Write-Host 'Open a new terminal and run: myloc --version'
Write-Host 'Keep the complete bundle at this location; rerun install.cmd if you move it.'
