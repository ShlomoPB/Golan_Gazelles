param(
  [string]$SourceMarkdown = "mountain_gazelle_manuscript_final_verification.md",
  [string]$OutputDocx = "Mountain_Gazelle_Manuscript_Final_Verification.docx"
)

Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

$workspace = (Get-Location).Path
$sourcePath = Join-Path $workspace $SourceMarkdown
$outputPath = Join-Path $workspace $OutputDocx

if (-not (Test-Path -LiteralPath $sourcePath)) {
  throw "Source markdown not found: $sourcePath"
}

Add-Type -AssemblyName System.IO.Compression

function Escape-XmlText {
  param([AllowNull()][string]$Text)
  if ($null -eq $Text) { return "" }
  return [System.Security.SecurityElement]::Escape($Text)
}

function Strip-MarkdownInline {
  param([string]$Text)
  $out = $Text -replace '\[([^\]]+)\]\(([^)]+)\)', '$1 ($2)'
  $out = $out -replace '\*\*([^*]+)\*\*', '$1'
  $out = $out -replace '\*([^*]+)\*', '$1'
  $out = $out -replace '`([^`]+)`', '$1'
  return $out
}

function New-Paragraph {
  param(
    [string]$Text,
    [ValidateSet("Normal", "Title", "Subtitle", "Heading1", "Heading2", "Bullet")]
    [string]$Style = "Normal"
  )

  $pPr = '<w:pPr><w:spacing w:after="120"/><w:jc w:val="both"/></w:pPr>'
  $rPr = '<w:rPr><w:rFonts w:ascii="Aptos" w:hAnsi="Aptos"/><w:sz w:val="21"/></w:rPr>'
  $prefix = ""

  switch ($Style) {
    "Title" {
      $pPr = '<w:pPr><w:spacing w:after="180"/><w:jc w:val="center"/></w:pPr>'
      $rPr = '<w:rPr><w:rFonts w:ascii="Aptos Display" w:hAnsi="Aptos Display"/><w:b/><w:color w:val="173F35"/><w:sz w:val="34"/></w:rPr>'
    }
    "Subtitle" {
      $pPr = '<w:pPr><w:spacing w:after="280"/><w:jc w:val="center"/></w:pPr>'
      $rPr = '<w:rPr><w:rFonts w:ascii="Aptos" w:hAnsi="Aptos"/><w:i/><w:color w:val="555555"/><w:sz w:val="20"/></w:rPr>'
    }
    "Heading1" {
      $pPr = '<w:pPr><w:keepNext/><w:spacing w:before="280" w:after="100"/></w:pPr>'
      $rPr = '<w:rPr><w:rFonts w:ascii="Aptos Display" w:hAnsi="Aptos Display"/><w:b/><w:color w:val="173F35"/><w:sz w:val="27"/></w:rPr>'
    }
    "Heading2" {
      $pPr = '<w:pPr><w:keepNext/><w:spacing w:before="220" w:after="80"/></w:pPr>'
      $rPr = '<w:rPr><w:rFonts w:ascii="Aptos Display" w:hAnsi="Aptos Display"/><w:b/><w:sz w:val="23"/></w:rPr>'
    }
    "Bullet" {
      $pPr = '<w:pPr><w:ind w:left="360" w:hanging="180"/><w:spacing w:after="70"/></w:pPr>'
      $prefix = [char]0x2022 + " "
    }
  }

  $clean = Escape-XmlText ($prefix + (Strip-MarkdownInline $Text))
  return "<w:p>$pPr<w:r>$rPr<w:t xml:space=`"preserve`">$clean</w:t></w:r></w:p>"
}

function Split-MarkdownTableRow {
  param([string]$Line)
  return @(($Line.Trim().Trim('|') -split '\|') | ForEach-Object { (Strip-MarkdownInline $_.Trim()) })
}

function New-WordTable {
  param(
    [string[]]$Headers,
    [object[]]$Rows
  )

  $xml = '<w:tbl><w:tblPr><w:tblW w:w="0" w:type="auto"/><w:tblLayout w:type="autofit"/><w:tblBorders><w:top w:val="single" w:sz="6" w:color="7C948C"/><w:left w:val="nil"/><w:bottom w:val="single" w:sz="6" w:color="7C948C"/><w:right w:val="nil"/><w:insideH w:val="single" w:sz="3" w:color="D5DDD9"/><w:insideV w:val="nil"/></w:tblBorders><w:tblCellMar><w:top w:w="70" w:type="dxa"/><w:left w:w="80" w:type="dxa"/><w:bottom w:w="70" w:type="dxa"/><w:right w:w="80" w:type="dxa"/></w:tblCellMar></w:tblPr>'
  $xml += '<w:tr><w:trPr><w:tblHeader/></w:trPr>'
  foreach ($header in $Headers) {
    $xml += '<w:tc><w:tcPr><w:shd w:fill="E8EFEC"/></w:tcPr><w:p><w:pPr><w:keepNext/></w:pPr><w:r><w:rPr><w:rFonts w:ascii="Aptos" w:hAnsi="Aptos"/><w:b/><w:sz w:val="18"/></w:rPr><w:t>' + (Escape-XmlText $header) + '</w:t></w:r></w:p></w:tc>'
  }
  $xml += '</w:tr>'

  foreach ($row in $Rows) {
    $xml += '<w:tr>'
    for ($i = 0; $i -lt $Headers.Count; $i++) {
      $value = if ($i -lt $row.Count) { [string]$row[$i] } else { "" }
      $xml += '<w:tc><w:p><w:r><w:rPr><w:rFonts w:ascii="Aptos" w:hAnsi="Aptos"/><w:sz w:val="18"/></w:rPr><w:t xml:space="preserve">' + (Escape-XmlText $value) + '</w:t></w:r></w:p></w:tc>'
    }
    $xml += '</w:tr>'
  }
  $xml += '</w:tbl><w:p><w:pPr><w:spacing w:after="80"/></w:pPr></w:p>'
  return $xml
}

$body = New-Object System.Collections.Generic.List[string]
$lines = Get-Content -LiteralPath $sourcePath
$paragraph = New-Object System.Collections.Generic.List[string]

function Flush-Paragraph {
  if ($paragraph.Count -gt 0) {
    $text = ($paragraph -join " ").Trim()
    if ($text.Length -gt 0) { $body.Add((New-Paragraph $text "Normal")) }
    $paragraph.Clear()
  }
}

$i = 0
while ($i -lt $lines.Count) {
  $line = $lines[$i]
  $trim = $line.Trim()

  if ($trim.Length -eq 0) {
    Flush-Paragraph
    $i++
    continue
  }

  if ($trim.StartsWith('|') -and $i + 1 -lt $lines.Count -and $lines[$i + 1].Trim() -match '^\|[\s:|-]+\|$') {
    Flush-Paragraph
    $headers = Split-MarkdownTableRow $trim
    $rows = New-Object System.Collections.Generic.List[object]
    $i += 2
    while ($i -lt $lines.Count -and $lines[$i].Trim().StartsWith('|')) {
      $rows.Add(@(Split-MarkdownTableRow $lines[$i]))
      $i++
    }
    $body.Add((New-WordTable -Headers $headers -Rows $rows.ToArray()))
    continue
  }

  if ($trim -match '^#\s+(.+)$') {
    Flush-Paragraph
    $body.Add((New-Paragraph $matches[1] "Title"))
    $i++
    if ($i -lt $lines.Count -and $lines[$i].Trim().Length -eq 0) { $i++ }
    if ($i -lt $lines.Count -and -not $lines[$i].Trim().StartsWith('#')) {
      $body.Add((New-Paragraph $lines[$i].Trim() "Subtitle"))
      $i++
    }
    continue
  }
  if ($trim -match '^##\s+(.+)$') {
    Flush-Paragraph
    $body.Add((New-Paragraph $matches[1] "Heading1"))
    $i++
    continue
  }
  if ($trim -match '^###\s+(.+)$') {
    Flush-Paragraph
    $body.Add((New-Paragraph $matches[1] "Heading2"))
    $i++
    continue
  }
  if ($trim -match '^[-*]\s+(.+)$') {
    Flush-Paragraph
    $body.Add((New-Paragraph $matches[1] "Bullet"))
    $i++
    continue
  }

  $paragraph.Add($trim)
  $i++
}
Flush-Paragraph

$documentXml = @"
<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<w:document xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main">
  <w:body>
    $($body -join "`n")
    <w:sectPr>
      <w:pgSz w:w="12240" w:h="15840"/>
      <w:pgMar w:top="1080" w:right="900" w:bottom="1080" w:left="900" w:header="720" w:footer="720"/>
    </w:sectPr>
  </w:body>
</w:document>
"@

$stylesXml = @"
<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<w:styles xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main">
  <w:style w:type="paragraph" w:default="1" w:styleId="Normal">
    <w:name w:val="Normal"/>
    <w:rPr><w:rFonts w:ascii="Aptos" w:hAnsi="Aptos"/><w:sz w:val="21"/></w:rPr>
  </w:style>
</w:styles>
"@

$contentTypesXml = @"
<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">
  <Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>
  <Default Extension="xml" ContentType="application/xml"/>
  <Override PartName="/word/document.xml" ContentType="application/vnd.openxmlformats-officedocument.wordprocessingml.document.main+xml"/>
  <Override PartName="/word/styles.xml" ContentType="application/vnd.openxmlformats-officedocument.wordprocessingml.styles+xml"/>
  <Override PartName="/docProps/core.xml" ContentType="application/vnd.openxmlformats-package.core-properties+xml"/>
  <Override PartName="/docProps/app.xml" ContentType="application/vnd.openxmlformats-officedocument.extended-properties+xml"/>
</Types>
"@

$packageRelsXml = @"
<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">
  <Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="word/document.xml"/>
  <Relationship Id="rId2" Type="http://schemas.openxmlformats.org/package/2006/relationships/metadata/core-properties" Target="docProps/core.xml"/>
  <Relationship Id="rId3" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/extended-properties" Target="docProps/app.xml"/>
</Relationships>
"@

$documentRelsXml = @"
<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">
  <Relationship Id="rIdStyles" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/styles" Target="styles.xml"/>
</Relationships>
"@

$created = (Get-Date).ToUniversalTime().ToString("s") + "Z"
$coreXml = @"
<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<cp:coreProperties xmlns:cp="http://schemas.openxmlformats.org/package/2006/metadata/core-properties" xmlns:dc="http://purl.org/dc/elements/1.1/" xmlns:dcterms="http://purl.org/dc/terms/" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">
  <dc:title>Mountain Gazelle Manuscript Final Verification Report</dc:title>
  <dc:creator>Mountain gazelle publication analysis workflow</dc:creator>
  <cp:lastModifiedBy>Mountain gazelle publication analysis workflow</cp:lastModifiedBy>
  <dcterms:created xsi:type="dcterms:W3CDTF">$created</dcterms:created>
  <dcterms:modified xsi:type="dcterms:W3CDTF">$created</dcterms:modified>
</cp:coreProperties>
"@

$appXml = @"
<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Properties xmlns="http://schemas.openxmlformats.org/officeDocument/2006/extended-properties">
  <Application>Mountain gazelle OOXML report builder</Application>
</Properties>
"@

$buildRoot = Join-Path $env:TEMP ("gazelle_verification_docx_" + [guid]::NewGuid().ToString("N"))
$null = New-Item -ItemType Directory -Path $buildRoot
$null = New-Item -ItemType Directory -Path (Join-Path $buildRoot "_rels")
$null = New-Item -ItemType Directory -Path (Join-Path $buildRoot "word")
$null = New-Item -ItemType Directory -Path (Join-Path $buildRoot "word\_rels")
$null = New-Item -ItemType Directory -Path (Join-Path $buildRoot "docProps")

try {
  Set-Content -LiteralPath (Join-Path $buildRoot "[Content_Types].xml") -Value $contentTypesXml -Encoding UTF8
  Set-Content -LiteralPath (Join-Path $buildRoot "_rels\.rels") -Value $packageRelsXml -Encoding UTF8
  Set-Content -LiteralPath (Join-Path $buildRoot "word\document.xml") -Value $documentXml -Encoding UTF8
  Set-Content -LiteralPath (Join-Path $buildRoot "word\styles.xml") -Value $stylesXml -Encoding UTF8
  Set-Content -LiteralPath (Join-Path $buildRoot "word\_rels\document.xml.rels") -Value $documentRelsXml -Encoding UTF8
  Set-Content -LiteralPath (Join-Path $buildRoot "docProps\core.xml") -Value $coreXml -Encoding UTF8
  Set-Content -LiteralPath (Join-Path $buildRoot "docProps\app.xml") -Value $appXml -Encoding UTF8

  if (Test-Path -LiteralPath $outputPath) { Remove-Item -LiteralPath $outputPath -Force }
  $stream = [System.IO.File]::Open($outputPath, [System.IO.FileMode]::Create)
  try {
    $archive = [System.IO.Compression.ZipArchive]::new($stream, [System.IO.Compression.ZipArchiveMode]::Create)
    try {
      foreach ($file in Get-ChildItem -LiteralPath $buildRoot -Recurse -File) {
        $entryName = $file.FullName.Substring($buildRoot.Length + 1).Replace('\', '/')
        $entry = $archive.CreateEntry($entryName, [System.IO.Compression.CompressionLevel]::Optimal)
        $entryStream = $entry.Open()
        $sourceStream = [System.IO.File]::OpenRead($file.FullName)
        try { $sourceStream.CopyTo($entryStream) } finally { $sourceStream.Dispose(); $entryStream.Dispose() }
      }
    } finally { $archive.Dispose() }
  } finally { $stream.Dispose() }
} finally {
  if (Test-Path -LiteralPath $buildRoot) { Remove-Item -LiteralPath $buildRoot -Recurse -Force }
}

Write-Host "Wrote $outputPath"
